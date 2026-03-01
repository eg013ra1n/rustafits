# Analysis API Usage Guide

## Quick Start

```rust
use rustafits::ImageAnalyzer;

let result = ImageAnalyzer::new()
    .analyze("light_001.fits")?;

println!("Stars detected: {}", result.stars_detected);
println!("Median FWHM:    {:.2} px", result.median_fwhm);
println!("Median HFR:     {:.2} px", result.median_hfr);
println!("Median SNR:     {:.1}", result.median_snr);
println!("Eccentricity:   {:.3}", result.median_eccentricity);
println!("SNR (dB):       {:.1}", result.snr_db);
println!("SNR weight:     {:.1}", result.snr_weight);
println!("PSF signal:     {:.1}", result.psf_signal);
```

`analyze()` accepts any path to a FITS or XISF file. It handles byte-swap,
u16-to-f32 conversion, and debayering (for OSC/Bayer images) automatically.

## Configuration

`ImageAnalyzer` uses a builder pattern. All settings have sensible defaults;
override only what you need.

```rust
let analyzer = ImageAnalyzer::new()
    .with_detection_sigma(3.0)
    .with_min_star_area(8)
    .with_max_stars(500)
    .without_gaussian_fit();

let result = analyzer.analyze("light_001.fits")?;
```

### Builder Methods

| Method | Default | Description |
|--------|---------|-------------|
| `with_detection_sigma(f32)` | 5.0 | Detection threshold in sigma above background. Lower values (2-3) find fainter stars but increase false positives. Clamped to minimum 1.0. |
| `with_min_star_area(usize)` | 5 | Minimum connected-component area in pixels. Filters hot pixels and noise. Clamped to minimum 1. |
| `with_max_star_area(usize)` | 2000 | Maximum connected-component area in pixels. Filters extended objects like galaxies and nebulae. |
| `with_saturation_fraction(f32)` | 0.95 | Reject stars with peak above this fraction of 65535. Clamped to 0.5-1.0. |
| `with_max_stars(usize)` | 200 | Keep only the brightest N stars. Clamped to minimum 1. |
| `without_gaussian_fit()` | Gaussian fit enabled | Use fast windowed-moments instead of iterative Gaussian fitting for FWHM/eccentricity. Less accurate but faster. |
| `with_background_mesh(usize)` | Global background | Enable mesh-grid background estimation with the given cell size in pixels. Handles uneven backgrounds and gradients. Cell size clamped to minimum 16. |
| `without_debayer()` | Debayer enabled | Skip green-channel interpolation for OSC images. By default, OSC images get native-resolution green interpolation + green-pixel-only fitting. This flag disables both, running analysis directly on the raw Bayer mosaic (faster but less accurate). |
| `with_trail_threshold(f32)` | 0.5 | R² threshold for trail detection (Path A). Lower = more aggressive. The raw `trail_r_squared` is always reported regardless of this setting. Clamped to 0.0-1.0. |
| `with_thread_pool(Arc<ThreadPool>)` | Global rayon pool | Route all parallel work through a custom rayon thread pool. |

### When to Adjust

- **Crowded fields** — lower `detection_sigma` to 3.0, raise `max_stars` to 500+.
- **Sparse fields / short exposures** — lower `detection_sigma` to 3.0-4.0 to pick up faint stars.
- **Hot pixel problems** — raise `min_star_area` to 8-10.
- **Nebula regions** — lower `max_star_area` to avoid detecting bright nebula knots as stars.
- **Speed over accuracy** — use `without_gaussian_fit()` and `without_debayer()`.
- **Strong gradients** — enable `with_background_mesh(64)` (typical cell sizes: 32-128).
- **OSC images** — green-channel interpolation and green-pixel-only fitting are applied automatically. No configuration needed. FWHM accuracy is within ~5% of PixInsight.

## Analyzing Raw Data

If you already have pixel data in memory, use `analyze_data()` to skip file I/O:

```rust
use rustafits::ImageAnalyzer;

// Planar f32 data in ADU range [0, 65535]
// For 3-channel: layout is RRRGGGBBB (channel planes stored contiguously)
let data: Vec<f32> = load_my_data();
let width = 4096;
let height = 2048;
let channels = 1; // 1 = mono, 3 = RGB

let result = ImageAnalyzer::new()
    .analyze_data(&data, width, height, channels)?;
```

`analyze_data()` skips file reading and debayering — the data must already be
in planar f32 format with values in the 0-65535 ADU range. Channel count must
be 1 (mono) or 3 (RGB in planar RRRGGGBBB layout).

## Reading Results

### `AnalysisResult`

| Field | Type | Description |
|-------|------|-------------|
| `width` | `usize` | Image width after debayer (if applicable). |
| `height` | `usize` | Image height after debayer (if applicable). |
| `source_channels` | `usize` | 1 = mono, 3 = color (after debayer). |
| `background` | `f32` | Global background level in ADU. |
| `noise` | `f32` | Background noise sigma in ADU. |
| `detection_threshold` | `f32` | Actual detection threshold used (ADU above background). |
| `stars_detected` | `usize` | Total stars found before the `max_stars` cap. |
| `stars` | `Vec<StarMetrics>` | Per-star metrics, sorted by flux descending, capped at `max_stars`. |
| `median_fwhm` | `f32` | Median FWHM across all measured stars (pixels). Good seeing: 2-3 px. |
| `median_eccentricity` | `f32` | Median eccentricity. 0 = perfectly round. Values > 0.5 suggest tracking issues. |
| `median_snr` | `f32` | Median per-star SNR. |
| `median_hfr` | `f32` | Median half-flux radius (pixels). Similar to FWHM but more robust to non-Gaussian profiles. |
| `snr_db` | `f32` | Image-wide SNR in decibels: `20 * log10(mean_signal / noise)`. Comparable to PixInsight SNRViews. |
| `snr_weight` | `f32` | PixInsight-style SNR weight: `(MeanDev / noise)^2`. Useful for weighting frames in stacking. |
| `psf_signal` | `f32` | PSF signal strength: `median(star_peaks) / noise`. Higher = better signal. |
| `trail_r_squared` | `f32` | Rayleigh R² statistic for directional coherence of star position angles. 0.0 = uniform (no trail), 1.0 = all aligned (strong trail). |
| `possibly_trailed` | `bool` | True if the image is likely trailed (R² above threshold or eccentricity-gated Rayleigh fires). |

### `StarMetrics`

| Field | Type | Description |
|-------|------|-------------|
| `x`, `y` | `f32` | Subpixel centroid position. |
| `peak` | `f32` | Background-subtracted peak value (ADU). |
| `flux` | `f32` | Total background-subtracted flux (ADU). |
| `fwhm_x` | `f32` | FWHM along major axis (pixels). |
| `fwhm_y` | `f32` | FWHM along minor axis (pixels). |
| `fwhm` | `f32` | Geometric mean FWHM: `sqrt(fwhm_x * fwhm_y)`. |
| `eccentricity` | `f32` | 0 = round, approaching 1 = elongated. |
| `snr` | `f32` | Per-star aperture photometry SNR. |
| `hfr` | `f32` | Half-flux radius (pixels). |
| `theta` | `f32` | PSF position angle in radians, counter-clockwise from +X axis. Orientation of the major axis (`fwhm_x` direction). With `theta`, `fwhm_x`, and `fwhm_y` you can reconstruct the full 2D Gaussian ellipse: `sigma_x = fwhm_x / 2.3548`, `sigma_y = fwhm_y / 2.3548`. |

## Thread Pool Usage

By default, parallel work uses the global rayon thread pool. To control
concurrency (e.g., to limit CPU usage or run multiple analyses concurrently),
pass a custom pool:

```rust
use std::sync::Arc;
use rustafits::{ImageAnalyzer, ThreadPoolBuilder};

let pool = Arc::new(
    ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("Failed to create thread pool")
);

let analyzer = ImageAnalyzer::new()
    .with_thread_pool(pool.clone());

let result = analyzer.analyze("light_001.fits")?;
```

`ThreadPool` and `ThreadPoolBuilder` are re-exported from `rayon` by
`rustafits`, so you don't need to add `rayon` as a direct dependency.

## Star Annotation Overlay

Analysis results can be visualized as color-coded ellipses on converted images.
See the [Annotation Documentation](annotation.md) for the full API reference,
including three tiers of integration (raw geometry, RGBA overlay layer, burn-in)
and configurable color thresholds.

Quick example:

```rust
use rustafits::{ImageConverter, ImageAnalyzer, annotate_image, AnnotationConfig};

let mut image = ImageConverter::new().process("light.fits")?;
let result = ImageAnalyzer::new().analyze("light.fits")?;

annotate_image(&mut image, &result, &AnnotationConfig::default());
ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```
