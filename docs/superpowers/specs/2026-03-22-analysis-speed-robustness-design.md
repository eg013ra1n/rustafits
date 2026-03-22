# Analysis Pipeline Speed & Robustness Optimization

**Date:** 2026-03-22
**Status:** Approved
**Goal:** Reduce per-frame analysis time from ~1.25s to ~0.5-0.6s on M4 while improving star selection robustness (reject nebula knots, edge stars, blended stars).

## Problem Statement

The star analysis pipeline has become slow (~1.25s/frame on M4, 2.5 minutes for 120 frames) due to accumulated complexity in the two-pass calibration pipeline. PSF fitting dominates at ~50% of wall time, measuring up to 2000 stars per frame when only ~100-200 well-measured stars are needed for stable aggregate statistics.

Additionally, the pipeline lacks morphological filters for non-stellar sources. Nebula knots, edge stars, and blended stars pass detection filters and either waste fitting time or corrupt FWHM/eccentricity statistics.

Athenaeum (the primary consumer) uses only the aggregate statistics (median FWHM, eccentricity, SNR) for frame ranking — it does not need deep per-star catalogs.

## Constraints

- Preserve PixInsight-comparison accuracy: FWHM R^2 >= 0.995, Ecc R^2 >= 0.943
- All final reported metrics come from PSF fits, never from moments
- measure_cap remains configurable for users who need deeper catalogs
- No changes to the public API surface (AnalysisResult struct, ImageAnalyzer builder)

## Design

### 1. Detection-Stage Filters (detection.rs)

Three new filters added during component validation in `detect_stars()`, after CCL and before the final star list is returned.

#### 1a. DAOFIND Sharpness

Computed from the convolved image (already allocated during detection):

```
sharpness = (conv[peak] - mean(conv[8_neighbors])) / conv[peak]
```

Discriminates:
- Stars: 0.3-0.9 (flux concentrated at peak relative to surroundings)
- Nebula knots: < 0.2 (diffuse; peak barely above neighbors in convolved space)
- Cosmic rays: > 0.9 (all flux in one pixel; convolution peak much higher than neighbors)

**Thresholds:** `sharplo = 0.2`, `sharphi = 0.9`. Hardcoded initially; may expose as DetectionParams fields later if tuning is needed.

**Cost:** 9 pixel lookups per candidate from the convolved image. Effectively zero.

**Implementation:** In the main detection loop (after `process_component` returns a star, before pushing to `stars` vec — around line 252 of detection.rs). Use the NMS peak pixel coordinates from `peak_map[label]` (not the centroid) to look up the 8 neighbors in the `conv` buffer, which is in scope at this point. Compute mean of neighbors, compute sharpness ratio, reject if outside bounds.

#### 1b. Concentration Index

Flux ratio at two radii within the detection stamp, computed from raw background-subtracted pixels:

```
CI = flux_within(1.0 * field_sigma) / flux_within(3.0 * field_sigma)
```

Where `field_sigma = field_fwhm / 2.3548`. For a Gaussian PSF, CI is approximately 0.39 regardless of FWHM. Deviations indicate non-stellar morphology:
- Nebula knots (diffuse): CI < 0.25
- Cosmic rays (sharp): CI > 0.7

**Threshold:** Initially use a fixed range `[0.15, 0.75]` which is intentionally generous — it covers the PSF variation from Gaussian (CI ~0.39) to Moffat with beta=2 (CI ~0.25-0.30). The wide range avoids false rejections; plan to tighten empirically after testing on the M42 dataset. The lower bound primarily catches nebula knots; the upper bound catches cosmic rays/hot pixels.

**Cost:** O(stamp_area) per candidate — traverses the same pixels already accessed for centroid/moments computation.

**Dependency:** Requires `field_fwhm` from Pass 1 calibration. Applied only in Pass 2 detection. Pass 1 runs without this filter (field_fwhm not yet known).

**Interface:** Add an `Option<f32>` parameter `field_fwhm` to `detect_stars()`. When `Some(fwhm)`, enable concentration index and edge margin filters. When `None` (Pass 1), skip them. This avoids confusion with the existing `fwhm` parameter which controls the matched filter kernel size (may differ from the calibrated PSF FWHM when field_fwhm is within 30% of 3.0).

**Implementation:** In the main detection loop (after `process_component` returns, before pushing to `stars`), construct a square stamp of radius `ceil(3 * field_sigma)` centered on the centroid. Compute two flux sums within circular apertures at radii `field_sigma` and `3*field_sigma`, using `data - background` (global) or `data - bg_map[pixel]` (adaptive) for background subtraction — matching whichever background source is used for detection. Reject if ratio is outside bounds.

#### 1c. Edge Margin

Replace the current 1-pixel border check (`touches_border`) with a centroid-based FWHM-proportional margin:

```
margin = max(2.0 * field_fwhm, 8.0)
```

Reject stars whose **centroid** is within `margin` of any image edge. This is a behavioral change from the current check, which tests whether any CCL **pixel** touches the image border. The centroid-based check is more appropriate because it ensures the full PSF core and fitting stamp are within bounds. The existing bounds check in `metrics.rs` (line 86-92) provides a secondary safety net for any edge stars that still reach measurement.

**Rationale:** The current check only rejects stars with pixels literally on the border. Stars 2-3 pixels from the edge produce incomplete stamps and unreliable fits. A margin of 2x FWHM ensures the full PSF core is measurable.

**Dependency:** Same as concentration index — requires `field_fwhm` via the `Option<f32>` parameter, applied in Pass 2 only. Pass 1 retains the original pixel-level border check.

**Impact:** Removes 1-3% additional stars near edges. These stars currently waste fitting time and produce high fit-residuals that dilute statistics.

### 2. Moments-Based Pre-Screening (metrics.rs)

A fast screening gate inserted into `measure_single_star()` before the Moffat/Gaussian/Moments fallback chain. Uses the stamp pixels already extracted for fitting.

#### The Gate

For each candidate passed to measurement:

1. Extract stamp (existing code)
2. Compute image moments: centroid (2 iterations) + second-order moments (Ixx, Iyy, Ixy)
3. Derive: moment_fwhm = 2.355 * sqrt(sigma_x * sigma_y), moment_ecc from eigenvalues, peak_sharpness = peak / (flux / area)
4. Apply reject criteria:
   - `moment_fwhm` outside `[0.7, 2.5] * field_fwhm` → reject
   - `moment_ecc > 0.8` → reject
   - `peak_sharpness` outside `[0.3, 5.0]` → reject
5. If rejected: return `None` (star excluded from results and statistics)
6. If accepted: proceed to Moffat/Gaussian fitting as normal

#### Rationale for Thresholds

- **FWHM [0.7, 2.5]x:** Below 0.7x = sub-PSF artifact (hot pixel cluster, cosmic ray). Above 2.5x = extended object (nebula knot, galaxy, badly blended group). Wide enough to accommodate natural PSF variation across the field (typically +/- 30%).
- **Eccentricity 0.8:** Matches the existing post-measurement filter. Catches elongated artifacts before wasting LM iterations. Bypassed when `possibly_trailed` is true (same logic as existing stat filter). **Edge case:** On frames with non-coherent trails (wind shake — high eccentricity but no angle coherence), the Rayleigh test may not flag `possibly_trailed`, so the 0.8 gate stays active. This is acceptable because most wind-shake stars have eccentricity 0.5-0.7 (below the gate). For extreme wind shake (ecc > 0.8 for most stars), add a safety fallback: if more than 50% of candidates fail the eccentricity gate, disable it and re-run measurement for the frame.
- **Peak/mean sharpness [0.3, 5.0]:** Reinforces detection-stage sharpness using raw pixels. Very diffuse sources have low peak-to-mean; single-pixel spikes have high.

#### Interface Change

`measure_stars()` receives `field_fwhm: f32` and `possibly_trailed: bool` as new parameters. Both are already available in `run_analysis()` at the call site.

#### What Happens to Rejected Stars

Rejected candidates are not measured, not included in the final star list, and not included in statistics. They don't fall through to the moments fallback — they are discarded. This is the primary speed win.

#### Cost

Moments: ~2N FLOPs per candidate (N = stamp pixel count).
LM fitting: ~1000N FLOPs per candidate (10 iterations average, 8 params).
**Ratio: ~500x cheaper per rejection.**

Expected rejection rate: 20-40% of candidates, concentrated on frames with nebulosity or other non-stellar features.

### 3. Reduced Measure Cap with Spatial Sampling (mod.rs)

#### Lower Default measure_cap

From 2000 to **500**.

**Statistical justification:** For a sigma-clipped weighted median with MAD-based outlier rejection, ~100-200 well-measured stars give 2% precision on median FWHM (1/sqrt(n) convergence). With 500 candidates and ~70% LM success rate after moments screening, ~350 measured stars reach statistics — well above the threshold. PixInsight recommends 500-2000; DES Y6 operates with median 164 per CCD.

**The measure_cap remains configurable.** Users needing deeper catalogs can still set it higher via `ImageAnalyzer::with_measure_cap()`.

#### Spatial Distribution Guard

Brightness-based selection risks spatial clustering of measured stars near the field center (where optical quality is best), systematically underestimating field-average FWHM.

**Solution:** Divide the image into a 4x4 grid (16 cells). Instead of taking the global brightest N stars, round-robin select from each cell by flux rank:

1. Bucket all detected stars into grid cells by centroid position
2. Sort each cell by flux descending (already globally sorted; bucketing preserves order). Add a `debug_assert!` that the detected list is flux-sorted before bucketing to guard against future sort-order changes in `detect_stars()`.
3. Round-robin: take one star from each non-empty cell, cycling until measure_cap is reached or all cells exhausted
4. Apply moments screening and LM fitting on the spatially-balanced set

**Cost:** O(N) bucketing + O(16) round-robin overhead. Negligible.

**Fallback:** If a cell has no stars (e.g., corner of a heavily vignetted frame), its share is redistributed to other cells naturally by the round-robin.

#### Net Speedup on Fitting

Current: up to 2000 LM fits × ~0.3ms each = ~600ms
After: ~500 candidates → ~350 LM fits × ~0.3ms = ~105ms
**~5x reduction in fitting time.**

## Performance Budget

Note: MRS wavelet convolution is already parallelized via `rayon::par_chunks_mut` in `b3_spline_smooth()` and `b3_spline_smooth_dilated()`. No further parallelization needed there.

| Component | Before | After | Savings |
|---|---|---|---|
| PSF fitting (LM) | ~625ms | ~105ms | ~520ms |
| Detection (2 passes) | ~312ms | ~300ms | ~12ms (sharpness/CI add negligible cost) |
| MRS wavelet + background | ~313ms | ~313ms | — (already parallelized) |
| **Total** | **~1250ms** | **~718ms** | **~532ms** |

**With moments pre-screening rejecting 20-40% before fitting:**
Fitting drops further to ~75-85ms. **Total ~0.5-0.6s per frame.**

For 120 frames: **~60-72 seconds** (down from 2.5 minutes). **~2-2.5x overall speedup.**

Conservative estimate — the actual speedup could be higher on frames with significant nebulosity (more candidates rejected by sharpness/CI/moments).

## Files Modified

| File | Changes |
|---|---|
| `src/analysis/detection.rs` | Add sharpness, concentration index, edge margin filters in component validation. Accept `field_fwhm` parameter for Pass 2. |
| `src/analysis/metrics.rs` | Add moments pre-screening gate in `measure_single_star()`. Accept `field_fwhm` and `possibly_trailed` parameters. |
| `src/analysis/mod.rs` | Lower default measure_cap to 500. Add spatial grid selection before measurement. Pass field_fwhm/possibly_trailed to measure_stars(). |
| `src/analysis/background.rs` | No changes needed (MRS convolution already parallelized). |

## Testing Strategy

1. **Accuracy regression:** Run on all 4 test files (cocoon.fits, mono.fits, osc.fits, test.xisf). Compare median_fwhm, median_eccentricity, median_hfr, median_snr against baseline values. Accept < 2% deviation.
2. **PixInsight comparison:** Re-run the 77-frame M42 comparison. Verify FWHM R^2 >= 0.99 and Ecc R^2 >= 0.93 (allowing small regression from current 0.995/0.943).
3. **Performance:** Time 120-frame batch before and after. Target: < 75 seconds on M4 (down from 150s).
4. **Edge cases:**
   - Frame with < 50 detected stars (sparse field): verify statistics still computed
   - Frame with heavy nebulosity (M42 core): verify nebula knots rejected
   - Trailed frame: verify moments eccentricity gate is bypassed
   - OSC (Bayer) frame: verify sharpness/CI work on debayered data
5. **Visual inspection:** Generate annotated outputs for all 4 test files. Verify rejected stars are indeed junk (no good stars lost).
6. **Unit tests for new filters:** Add unit tests in `detection.rs` and `metrics.rs` using synthetic data:
   - Sharpness: synthetic Gaussian star (expect ~0.5), synthetic nebula blob (expect < 0.2), synthetic cosmic ray (expect > 0.9)
   - Concentration index: Gaussian star (CI ~0.39), diffuse blob (CI < 0.2), sharp spike (CI > 0.7)
   - Edge margin: stars at various distances from the border
   - Moments screening: stars with known FWHM relative to field_fwhm, verify correct accept/reject at threshold boundaries
   - Reuse existing `make_star_field` pattern in detection.rs tests for synthetic data generation.

## Risks

1. **Over-filtering on unusual frames:** Tight sharpness/CI/moments thresholds might reject real stars on severely defocused or heavily trailed frames. Mitigation: thresholds are relative to field_fwhm (adapts to conditions); eccentricity gate bypassed on trailed frames.
2. **Spatial grid on sparse fields:** If < 500 stars detected, all go to measurement regardless. Grid only affects distribution within the cap.
3. **Accuracy regression from fewer measured stars:** Statistical theory and DES Y6 experience show 200+ stars is sufficient. Fit-residual weighting and sigma-clipped median provide additional robustness.
