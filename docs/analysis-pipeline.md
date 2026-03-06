# Analysis Pipeline Overview

The analysis module performs automatic star detection, PSF measurement, and image quality
assessment on astronomical images. It outputs per-star metrics and image-wide statistics
comparable to PixInsight's SubframeSelector.

All detected stars are measured and contribute to statistics — no artificial cap on
measurement count. Like PixInsight, statistics reflect the full stellar population.

## Data Flow

```
FITS / XISF File
       |
       v
+------------------+
|  Read & Decode   |  u16 big-endian (FITS) or f32 little-endian (XISF)
+------------------+
       |
       v
+------------------+
| u16 -> f32       |  SIMD-accelerated, [0..65535] range
+------------------+
       |
       v
+-------------------------------+
| To Mono Luminance             |  Three mutually exclusive paths:
|                               |
| Mono: pass through            |  Already single-channel, nothing to do
|                               |
| OSC (Bayer): Green Interp     |  Native-resolution green-channel interpolation
|   green CFA pixels keep       |  R/B pixels get distance-weighted average of
|   original values              |  green neighbors (3x3). Builds green_mask
|                               |  so fitter only uses true green pixels
|                               |
| RGB (3-ch): Extract Lum       |  L = 0.2126R + 0.7152G + 0.0722B
|   (pre-debayered input via    |  Converts 3-channel planar to single-channel
|    analyze_data only)          |
+-------------------------------+
       |
       v
+-------------------------------+
| Background Estimation         |  Global (sigma-clipped) or Mesh-grid (spatially varying)
| -> background level (ADU)     |  Optional iterative source-masked re-estimation
| -> noise sigma (ADU)          |  (with_iterative_background + with_background_mesh)
| -> optional bg/noise maps     |
+-------------------------------+
       |
       v
+-------------------------------+
| Two-Pass Adaptive Detection   |
|                               |
| Pass 1: Separable Gaussian    |  SIMD-accelerated convolution (AVX2/SSE2/NEON)
|   convolution with FWHM=3.0   |  Zero-sum DAOFIND-style matched filter kernel
|   + CCL + peak deblending     |  Connected Component Labeling with Union-Find
|                               |  Voronoi deblending for multi-peak components
| -> estimate FWHM from top-20  |  Stamp-based half-max radial profile
|                               |
| Pass 2 (if FWHM differs >30%)|  Re-detect with measured FWHM kernel
|   convolution with measured   |  Adapts to actual seeing conditions
|   FWHM + CCL + deblending    |
+-------------------------------+
       |
       | (optional: iterative background loop — mask sources, re-estimate, re-detect)
       v
+-------------------------------+
| Trail Detection               |  Rayleigh test on position angles (2 theta)
| (image-level, advisory)       |  Dual-path: strong R^2 OR eccentricity-gated
| -> trail_r_squared            |  Advisory only — never rejects, always continues
| -> possibly_trailed flag      |  Computed on ALL detected candidates
+-------------------------------+
       |
       v
+-------------------------------+
| PSF Measurement (ALL stars)   |  Extract stamp around each detected star
| -> FWHM (x, y, geometric)    |  2D Gaussian LM fit (default) or windowed moments
| -> Eccentricity               |  Optional: 2D Moffat LM fit (8-param, reports beta)
| -> HFR, theta                 |  Moments-based theta (always computed)
| -> Moffat beta (optional)     |  OSC: only green CFA pixels fed to fitter
|                               |  Stars that fail fitting are excluded
+-------------------------------+
       |
       v
+-------------------------------+
| Per-Star SNR (ALL measured)   |  Aperture photometry with local sky annulus
| -> aperture flux              |  CCD noise equation: SNR = F / sqrt(F + n*sigma^2)
| -> local background           |  Aperture radius = 2.5 * median_fwhm
+-------------------------------+
       |
       v
+-------------------------------+
| Statistics (full population)  |  Computed from ALL measured stars, not a subset
| -> median FWHM                |
| -> median eccentricity        |  stars_detected = total with valid measurements
| -> median SNR, HFR            |
| -> median Moffat beta         |
+-------------------------------+
       |
       v
+-------------------------------+
| Image-Wide Metrics            |
| -> SNR dB                     |  20 * log10(mean / noise)  -- PI SNRViews
| -> SNR Weight                 |  (MeanDev / noise)^2       -- PI SubframeSelector
| -> PSF Signal                 |  median(star_peaks) / noise
+-------------------------------+
       |
       v
+-------------------------------+
| Late Cap (display only)       |  measured.truncate(max_stars)
| -> stars: Vec<StarMetrics>    |  Only limits the returned per-star vector
|                               |  Statistics already computed from full population
+-------------------------------+
       |
       v
+-------------------------------+
| AnalysisResult                |  Summary medians + capped per-star vector
+-------------------------------+
```

## Detection Filters

Stars are filtered during detection (before measurement) by:

| Filter | Threshold | Rejects |
|--------|-----------|---------|
| Matched filter SNR | detection_sigma * noise * sqrt(kernel_energy) | Faint sources below threshold |
| Non-maximum suppression | Spatial radius = kernel half-width | Duplicate peaks |
| Min area | min_star_area (default 5 px) | Hot pixels, noise spikes |
| Max area | max_star_area (default 2000 px) | Extended objects (galaxies, nebulae, comae) |
| Saturation | saturation_fraction * 65535 (default 0.95) | Saturated stars |
| Border | Touches image edge | Truncated PSFs |
| Aspect ratio | Bounding box > 8:1 | Cosmic rays (aspect ~20-100), satellite trails (~50+) |

No eccentricity filter is applied — all stars with valid measurements contribute
to statistics, matching PixInsight's approach.

## Module Map

| Module           | File              | Purpose                                      |
|------------------|-------------------|----------------------------------------------|
| Background       | `background.rs`   | Global, mesh-grid, and iterative background   |
| Convolution      | `convolution.rs`  | Separable Gaussian convolution (SIMD dispatch) |
| Detection        | `detection.rs`    | DAOFIND matched filter + CCL + deblending     |
| Fitting          | `fitting.rs`      | LM 2D Gaussian (7-param) and Moffat (8-param) |
| Metrics          | `metrics.rs`      | Per-star FWHM, eccentricity, HFR measurement |
| SNR              | `snr.rs`          | Per-star and image-wide SNR computations      |
| Orchestration    | `mod.rs`          | Builder API, trail detection, pipeline wiring |

## Trail Detection

Detects directional coherence across all detected stars using circular statistics
on position angles. Computes an R^2 statistic and advisory `possibly_trailed` flag.
The pipeline **always continues** to PSF measurement regardless of the flag — the
caller decides whether to reject based on `trail_r_squared` and `possibly_trailed`.

Two detection paths cover different regimes:

| Path | Condition | Catches |
|------|-----------|---------|
| A — Strong coherence | R^2 > threshold AND p < 0.01 | Trails with consistent theta (low-ecc knots) |
| B — Eccentricity-gated | median ecc > 0.6 AND p < 0.05 | Undersampled trails (elongated blobs) |

The Path A threshold defaults to 0.5 and is configurable via `with_trail_threshold()`.

See [Trail Detection](trail-rejection.md) for the full algorithm and rationale.

## PixInsight Comparison

| Metric              | PixInsight Tool            | Our Equivalent        | mono.fits | osc.fits  |
|---------------------|----------------------------|-----------------------|-----------|-----------|
| FWHM                | FWHMEccentricity           | `median_fwhm`         | ~2.16 px  | ~2.66 px  |
| FWHM (PI)           | FWHMEccentricity           | —                     | ~2.16 px  | ~2.73 px  |
| Eccentricity        | FWHMEccentricity           | `median_eccentricity` | ~0.49     | ~0.44     |
| Stars detected      | SubframeSelector           | `stars_detected`      | ~4,485    | ~712      |
| Image SNR           | SNRViews (dB)              | `snr_db`              | ~27.7 dB  | —         |
| SNR Weight          | SubframeSelector           | `snr_weight`          | ~1.1      | —         |
| PSF Signal          | SubframeSelector           | `psf_signal`          | ~14       | —         |

FWHM accuracy: mono ~1% vs PixInsight, OSC ~3% vs PixInsight.

See individual algorithm documents for full details:

- [Background Estimation](background.md)
- [Star Detection](detection.md)
- [Trail Rejection (Rayleigh Test)](trail-rejection.md)
- [Gaussian & Moffat Fitting](fitting.md)
- [PSF Metrics (FWHM, Eccentricity, HFR)](metrics.md)
- [SNR Computations](snr.md)
- [Star Annotation Overlay](annotation.md)
