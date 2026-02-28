# Analysis Pipeline Overview

The analysis module performs automatic star detection, PSF measurement, and image quality
assessment on astronomical images. It outputs per-star metrics and image-wide statistics
comparable to PixInsight.

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
+------------------+
| Debayer (if OSC) |  Super-pixel: 2x2 Bayer -> 1 RGB pixel (half resolution)
+------------------+
       |
       v
+------------------+
| Extract Lum      |  L = 0.2126R + 0.7152G + 0.0722B  (ITU Rec.709)
| (if RGB)         |  Mono images skip this step
+------------------+
       |
       v
+----------------------------+
| Background Estimation      |  Global (sigma-clipped) or Mesh-grid (spatially varying)
| -> background level (ADU)  |
| -> noise sigma (ADU)       |
| -> optional background map |
+----------------------------+
       |
       v
+----------------------------+
| Star Detection             |  DAOFIND matched filter + Connected Component Labeling
| -> candidate list          |  Stamp-based theta & eccentricity per star
| -> sorted by flux          |
+----------------------------+
       |
       v
+----------------------------+
| Trail Rejection            |  Rayleigh test on position angles (2 theta)
| (Phase 1: image-level)     |  Dual-path: strong R^2 OR eccentricity-gated
| -> zero result if trailing |
+----------------------------+
       |
       v
+----------------------------+
| PSF Measurement (per star) |  Extract stamp around each star
| -> FWHM (x, y, geometric) |  2D Gaussian fit (default) or windowed moments
| -> Eccentricity            |  Moments-based theta (always)
| -> HFR, theta              |
+----------------------------+
       |
       v
+----------------------------+
| Eccentricity Filter        |  Reject stars with ecc > max_eccentricity (default 0.5)
| (Phase 2: per-star)        |  Removes remaining cosmic rays, satellite trails
+----------------------------+
       |
       v
+----------------------------+
| Per-Star SNR               |  Aperture photometry with local sky annulus
| -> aperture flux           |  CCD noise equation: SNR = F / sqrt(F + n*sigma^2)
| -> local background        |
+----------------------------+
       |
       v
+----------------------------+
| Image-Wide Metrics         |
| -> SNR dB                  |  20 * log10(mean / noise)  -- comparable to PI SNRViews
| -> SNR Weight              |  (MeanDev / noise)^2       -- comparable to PI SubframeSelector
| -> PSF Signal              |  median(star_peaks) / noise
+----------------------------+
       |
       v
+----------------------------+
| AnalysisResult             |  Summary medians + per-star vector
+----------------------------+
```

## Module Map

| Module           | File              | Purpose                                      |
|------------------|-------------------|----------------------------------------------|
| Background       | `background.rs`   | Global and mesh-grid background estimation    |
| Detection        | `detection.rs`    | DAOFIND matched filter + connected components |
| Fitting          | `fitting.rs`      | Levenberg-Marquardt 1D and 2D Gaussian fits   |
| Metrics          | `metrics.rs`      | Per-star FWHM, eccentricity, HFR measurement |
| SNR              | `snr.rs`          | Per-star and image-wide SNR computations      |
| Orchestration    | `mod.rs`          | Builder API, trail rejection, pipeline wiring |

## Trail Rejection

The pipeline has a two-phase strategy for rejecting trailed images:

**Phase 1 — Image-level Rayleigh test** (in `mod.rs`, before PSF measurement):

Detects directional coherence across all detected stars using circular statistics
on position angles. If the position angles of all detected blobs point the same
direction, the image is trailed and the pipeline returns a zero result (no stars).

Two rejection paths cover different regimes:

| Path | Condition | Catches |
|------|-----------|---------|
| A — Strong coherence | R^2 > 0.3 AND p < 0.01 | Oversampled trails (low-ecc knots, consistent theta) |
| B — Eccentricity-gated | median ecc > 0.6 AND p < 0.05 | Undersampled trails (elongated blobs) |

See [Trail Rejection](trail-rejection.md) for the full algorithm and rationale.

**Phase 2 — Per-star eccentricity filter** (after PSF measurement):

Stars with eccentricity above `max_eccentricity` (default 0.5) are removed from the
final result. This catches individual trailed stars, cosmic rays, and other elongated
artifacts that survived Phase 1.

## PixInsight Comparison

| Metric              | PixInsight Tool            | Our Equivalent        | mono.fits |
|---------------------|----------------------------|-----------------------|-----------|
| FWHM                | FWHMEccentricity           | `median_fwhm`         | ~2.16 px  |
| Eccentricity        | FWHMEccentricity           | `median_eccentricity` | ~0.35     |
| Image SNR           | SNRViews (dB)              | `snr_db`              | ~27.7 dB  |
| SNR Weight          | SubframeSelector           | `snr_weight`          | ~1.1      |
| PSF Signal          | SubframeSelector           | `psf_signal`          | ~383      |

See individual algorithm documents for full details:

- [Background Estimation](background.md)
- [Star Detection](detection.md)
- [Trail Rejection (Rayleigh Test)](trail-rejection.md)
- [Gaussian Fitting](fitting.md)
- [PSF Metrics (FWHM, Eccentricity, HFR)](metrics.md)
- [SNR Computations](snr.md)
