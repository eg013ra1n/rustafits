# Analysis Pipeline Redesign

Date: 2026-03-08

## Goal

Redesign the analysis pipeline to be robust, precise, and fast with a single unified codepath. Eliminate algorithm-selection config knobs — every remaining option tunes a threshold, not a codepath.

## Use Cases

1. **Subframe scoring** — FWHM, eccentricity, SNR for keep/reject decisions
2. **Plate solving** — precise sub-pixel centroids for catalog matching
3. **Optical diagnostics** — per-star PSF data for spatial quality analysis (tilt, curvature, coma)

## Pipeline Architecture

```
Input → Luminance → Background (mesh+MRS)
  → Pass 1 (discover PSF model from bright stars)
  → Pass 2 (full detection + fixed-beta Moffat)
  → Metrics (per-star + frame-level)
  → AnalysisResult
```

Single unified pipeline, no mode switching.

### Eliminated Codepaths

- Global background → always mesh-grid, auto-tuned cell size
- MAD noise → always MRS wavelet
- Moments-only mode → always fit (fallback chain)
- Gaussian-only mode → Moffat primary

### Remaining Config

```rust
ImageAnalyzer::new()
    .with_detection_sigma(f32)       // default: 5.0
    .with_max_stars(usize)           // cap returned stars
    .with_min_star_area(usize)       // default: 5
    .with_max_star_area(usize)       // default: 2000
    .with_saturation_fraction(f32)   // default: 0.95
    .with_mrs_layers(usize)          // default: 1
    .with_trail_threshold(f32)       // default: 0.5
    .without_debayer()               // skip green interpolation
    .with_thread_pool(Arc<ThreadPool>)
```

No "choose your algorithm" knobs. Every option tunes a threshold or a cap.

## Stage 1: Background & Noise Estimation

One path, always:

1. **Auto-tune cell size** — `max(16, longest_axis / 32)`, ensuring ~1000+ pixels per cell
2. **Per-cell sigma-clipped statistics** — 3 iterations, 3-sigma MAD clipping. Reject cells with >30% clipped.
3. **Fill invalid cells** — nearest-neighbor search
4. **3x3 median filter** on grid — suppress residual star contamination
5. **Bicubic Catmull-Rom interpolation** to full resolution
6. **MRS wavelet noise** — B3-spline, configurable layers (default 1). Per-cell noise interpolated bilinearly.
7. **Source mask iteration** — after Pass 1 detection, mask stars (r = 2.5 x FWHM), re-estimate background. One mandatory iteration.

Global background/noise derived as aggregates of the mesh (median of cell values), not a separate codepath.

### Dilated a trous wavelet (new, for layers 2+)

Same B3 kernel `[1/16, 1/4, 3/8, 1/4, 1/16]` with spacing `2^(layer-1)` between taps. Separable passes with SIMD. Layer 1 already exists; layers 2+ insert zeros between kernel taps (dilated convolution).

## Stage 2: Detection & PSF Model Discovery

### Pass 1 — Discovery

1. DAOFIND matched filter with sigma=3.0 px kernel
2. CCL + deblending (Voronoi for pass 1 — only need rough centroids)
3. Select calibration stars — brightest ~50-100, filtered by:
   - Not saturated
   - Not blended (single peak in CCL blob)
   - Not on border
   - Eccentricity < 0.5
4. Free-beta Moffat fit on each calibration star
5. Derive field-wide PSF model: sigma-clipped median of beta and FWHM
6. Source mask → trigger background re-estimation (Stage 1, step 7)

### Pass 2 — Full Detection

1. DAOFIND with measured FWHM as kernel size
2. CCL + deblending with iterative PSF subtraction:
   - Single-peak blobs → fit directly
   - Multi-peak blobs → sort by brightness, fit brightest, subtract model, fit next on residual
3. Fixed-beta Moffat (field median beta) on all stars
4. Fallback chain per star: fixed-beta Moffat → Gaussian → windowed moments (flagged unreliable)
5. Detection depth controlled by `detection_sigma`

### Iterative PSF Subtraction (new, replaces Voronoi for Pass 2)

1. Sort blended peaks by brightness (descending)
2. Fit Moffat on brightest peak using full stamp data
3. Evaluate fitted model on pixel grid, subtract from stamp
4. Fit next peak on residual
5. Repeat (rarely more than 2-3 stars per blend)

Requires exposing `evaluate_moffat_2d(params, x, y) -> f32` from fitting.rs.

## Stage 3: Fitting Engine

Levenberg-Marquardt with analytical Jacobian, f64 arithmetic.

### Convergence Parameters

| Parameter | Value |
|---|---|
| Max iterations | 50 |
| Convergence tolerance | 1e-6 |
| sigma/alpha minimum | 0.5 px |
| sigma/alpha maximum | 5x initial |
| Beta range | 1.5 - 10.0 |
| Amplitude check | A > 3 x local_noise |
| Centroid drift | < 2 px from initial |

### Lambda Strategy

Keep existing Nielsen gain ratio adaptation. Add Cholesky failure retry: multiply lambda by 10, retry up to 3 times before declaring non-convergence.

### Convergence Check

Dual criterion:
- Parameter norm: `||delta|| / ||params|| < 1e-6`
- Residual: `|delta_chi2 / chi2| < 1e-4`

Stop when either is satisfied.

### Fallback Chain

```
Free-beta Moffat (8 params)
  → Fixed-beta Moffat (7 params, field median beta)
    → Gaussian (7 params)
      → Windowed moments (flagged unreliable)
```

## Stage 4: Metrics & Output

### Per-Star Metrics

- Centroid (x, y) — sub-pixel, from fit
- FWHM_x, FWHM_y — from alpha and beta
- FWHM — geometric mean
- Eccentricity — from fit axes
- Position angle theta — from fit rotation
- HFR — half-flux radius (aperture-based)
- Aperture SNR
- Beta — from fit (or field median if fixed)
- Fit method — enum: FreeMoffat | FixedMoffat | Gaussian | Moments
- Converged — bool

### Frame-Level Metrics

- Background (median of mesh)
- Noise (MRS sigma)
- Star count
- Median FWHM (sigma-clipped)
- Median eccentricity (sigma-clipped)
- Background SNR — mean_signal / noise
- SNR Weight — median(star_flux)^2 / (noise^2 * background)
- PSF Signal — median(star_peaks) / noise
- Trail detection — Rayleigh R^2 test (advisory flag)

### Removed Metrics

- SNR dB — redundant with background SNR

## Implementation Summary

### New Work

1. **Dilated a trous wavelet** — MRS layers 2+ in convolution.rs
2. **Iterative PSF subtraction** — in detection.rs, evaluate + subtract fitted Moffat model
3. **Auto-tuned mesh cell size** — in background.rs
4. **Calibration star selection** — in mod.rs, pass 1 filtering + field-wide beta/FWHM derivation
5. **FitMethod enum** — in types.rs

### Refactoring

6. **Merge measurement functions** — metrics.rs: one `measure_star()` with internal fallback chain
7. **Simplify orchestration** — mod.rs: remove config branching, single linear flow
8. **Fitting convergence** — fitting.rs: dual criterion, Cholesky retry, tightened bounds, centroid drift check

### Unchanged

- convolution.rs — DAOFIND kernel + separable convolution (layer 1 B3 spline)
- snr.rs — keep aperture SNR, SNR weight, PSF signal. Remove only `compute_snr_db()`
- render.rs — debug visualization
- annotate.rs — consumes AnalysisResult

### Config Methods Removed

| Removed | Replacement |
|---|---|
| `.without_gaussian_fit()` | Always Moffat (Gaussian is fallback only) |
| `.with_background_mesh(cell_size)` | Always mesh, auto-tuned |
| `.with_iterative_background(n)` | Always one iteration |
| `.with_moffat_beta(f32)` | Beta auto-derived from pass 1 |
| `.with_max_distortion(f32)` | Internal to calibration star selection |
| `.with_moffat_fit()` / `.without_moffat_fit()` | Always Moffat |
