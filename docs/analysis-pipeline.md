# Analysis Pipeline Overview

The analysis module performs automatic star detection, PSF measurement, and image quality
assessment on astronomical images. It outputs per-star metrics and image-wide statistics
suitable for subframe evaluation and stacking weight computation.

A configurable measure cap (default 2000) limits how many stars undergo PSF fitting.
Stars are sorted by flux (brightest first) before capping. Statistics reflect the
measured population, while `stars_detected` reports the raw detection count.

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
| Background Estimation         |  Mesh-grid (spatially varying)
| -> background level (ADU)     |  Cell size = auto_cell_size(w,h) = max(16, max(w,h)/32)
| -> noise sigma (ADU)          |  Bicubic interpolation between cell medians
| -> bg/noise maps              |
|                               |
| Noise: MRS wavelet (default)  |  B3-spline à trous wavelet, 4-layer significance masking
|   (configurable layers via    |  Isolates noise from nebulosity/gradients
|    with_mrs_layers)            |
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
| Pass 2 (if FWHM differs >30%)|  Re-detect with calibration FWHM kernel
|   convolution with calibration|  Uses field_fwhm from calibration pass
|   FWHM + CCL + deblending    |
+-------------------------------+
       |
       v
+-------------------------------+
| Calibration Pass              |  Free-beta Moffat on top 100 bright stars
|                               |  (eccentricity < 0.5, area >= 5 px)
| -> field_beta                 |  Sigma-clipped median of fitted beta values
| -> field_fwhm                 |  Sigma-clipped median of fitted FWHM values
+-------------------------------+
       |
       v
+-------------------------------+
| Source-Mask Background        |  Mask star pixels (r = 2.5 * field_fwhm)
| Re-estimation                 |  Re-run mesh-grid background on masked image
| -> refined background map     |  Removes stellar flux contamination
| -> refined noise map          |
+-------------------------------+
       |
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
| Measure Cap                   |  Slice detected stars to measure_cap (default 2000)
| stars_detected = raw count    |  Stars already sorted by flux (brightest first)
| 0 = measure all               |  Brightest stars have highest-SNR fits
+-------------------------------+
       |
       v
+-------------------------------+
| PSF Measurement (capped)      |  Extract stamp around each star to measure
| -> FWHM (x, y, geometric)    |  Fixed-beta Moffat (pass 2): beta = field_beta
| -> Eccentricity               |  Fallback chain: Moffat -> Gaussian -> Moments
| -> HFR, theta                 |  FitMethod tracked per star (Moffat/Gaussian/Moments)
| -> Moffat beta                |  Moments-based theta (always computed)
|                               |  OSC: only green CFA pixels fed to fitter
|                               |  Stars that fail all fitting are excluded
|                               |  LM: configurable max_iter/conv_tol/max_rejects
+-------------------------------+
       |
       v
+-------------------------------+
| Per-Star SNR (all measured)   |  Aperture photometry with local sky annulus
| -> aperture flux              |  CCD noise equation: SNR = F / sqrt(F + n*sigma^2)
| -> local background           |  Aperture radius = 2.5 * median_fwhm
+-------------------------------+
       |
       v
+-------------------------------+
| Statistics (measured stars)   |  Computed from all measured stars (after cap)
| -> median FWHM                |
| -> median eccentricity        |  stars_detected = raw detection count (before cap)
| -> median SNR, HFR            |
| -> median Moffat beta         |
+-------------------------------+
       |
       v
+-------------------------------+
| Image-Wide Metrics            |
| -> SNR Weight                 |  (MeanDev / noise)^2  -- subframe ranking
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

All measured stars (after measure cap) contribute to statistics.

## Module Map

| Module           | File              | Purpose                                      |
|------------------|-------------------|----------------------------------------------|
| Background       | `background.rs`   | Mesh-grid background, MRS wavelet noise        |
| Convolution      | `convolution.rs`  | Separable Gaussian convolution (SIMD), B3-spline smoothing |
| Detection        | `detection.rs`    | DAOFIND matched filter + CCL + deblending     |
| Fitting          | `fitting.rs`      | Two-pass calibration, Moffat->Gaussian->Moments fallback |
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

## Typical Results

| Metric              | Field               | mono.fits | osc.fits  |
|---------------------|---------------------|-----------|-----------|
| FWHM                | `median_fwhm`       | ~2.16 px  | ~2.66 px  |
| Eccentricity        | `median_eccentricity` | ~0.49   | ~0.44     |
| Stars detected      | `stars_detected`    | ~4,485    | ~712      |
| SNR Weight          | `snr_weight`        | ~1.1      | —         |
| PSF Signal          | `psf_signal`        | ~14       | —         |
| Frame SNR           | `frame_snr`         | —         | —         |

FWHM accuracy: within ~0.3% of professional tools.

## Configurable Parameters

The pipeline exposes a small set of tuning knobs via the builder API. Most algorithm
choices (mesh-grid background, MRS noise, two-pass calibration, Moffat fallback chain)
are now fixed defaults.

| Parameter | Default | Builder Method | Effect |
|-----------|---------|----------------|--------|
| Detection sigma | 5.0 | `with_detection_sigma(f32)` | Matched-filter SNR threshold |
| Min star area | 5 px | `with_min_star_area(u32)` | Reject sources smaller than N pixels |
| Max star area | 2000 px | `with_max_star_area(u32)` | Reject extended objects larger than N pixels |
| Saturation fraction | 0.95 | `with_saturation_fraction(f32)` | Fraction of 65535 above which stars are rejected |
| Max stars | 5000 | `with_max_stars(usize)` | Late cap on returned per-star vector (statistics use all) |
| Measure cap | 2000 | `with_measure_cap(usize)` | Max stars to PSF-fit (0 = all). Dense fields benefit most. |
| Fit max iter | 25 | `with_fit_max_iter(usize)` | LM max iterations for measurement pass (calibration always 50) |
| Fit tolerance | 1e-4 | `with_fit_tolerance(f64)` | LM convergence tolerance for measurement pass (calibration always 1e-6) |
| Fit max rejects | 5 | `with_fit_max_rejects(usize)` | LM consecutive reject bailout |
| MRS noise layers | 4 | `with_mrs_layers(usize)` | Wavelet layers for iterative significance masking |
| Trail threshold | 0.5 | `with_trail_threshold(f32)` | R^2 threshold for trail advisory flag |
| Debayer | auto | `with_debayer(bool)` | Force debayer on/off for OSC data |
| Thread pool | None | `with_thread_pool(Arc<ThreadPool>)` | Optional shared Rayon thread pool |

See individual algorithm documents for full details:

- [Background Estimation](background.md)
- [Star Detection](detection.md)
- [Trail Rejection (Rayleigh Test)](trail-rejection.md)
- [Gaussian & Moffat Fitting](fitting.md)
- [PSF Metrics (FWHM, Eccentricity, HFR)](metrics.md)
- [SNR Computations](snr.md)
- [Star Annotation Overlay](annotation.md)
