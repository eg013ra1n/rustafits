# Athenaeum Integration Guide

## Overview

rustafits provides three pipelines that Athenaeum can use independently or together:

1. **Image Conversion** — FITS/XISF to RGB/RGBA pixel data (`ImageConverter`)
2. **Image Analysis** — star detection, FWHM, eccentricity, SNR (`ImageAnalyzer`)
3. **Star Annotation** — color-coded ellipse overlays on converted images (`annotate_image`, `compute_annotations`, `create_annotation_layer`)

All three share a common thread pool via `with_thread_pool()`.

---

## 1. Image Conversion

```rust
use astroimage::{ImageConverter, ProcessedImage};

let image: ProcessedImage = ImageConverter::new()
    .with_downscale(2)
    .with_rgba_output()        // RGBA for Canvas/ImageData
    .with_thread_pool(pool.clone())
    .process("light.fits")?;

// image.data      — Vec<u8>, interleaved RGBA bytes
// image.width     — pixel width
// image.height    — pixel height
// image.channels  — 4 (RGBA) or 3 (RGB)
// image.is_color  — true if debayered/RGB
// image.flip_vertical — true if vertically flipped during processing
```

### Saving a ProcessedImage to disk

```rust
ImageConverter::save_processed(&image, "output.jpg", 95)?;
ImageConverter::save_processed(&image, "output.png", 0)?;  // quality ignored for PNG
```

### Builder Methods

| Method | Default | Description |
|--------|---------|-------------|
| `with_downscale(n)` | 1 | Downscale by factor n. Bayer images: debayer counts as 2x. |
| `with_quality(q)` | 95 | JPEG quality 1-100. |
| `without_debayer()` | debayer on | Skip Bayer debayering. |
| `with_preview_mode()` | off | 2x2 binning for fast previews. |
| `with_rgba_output()` | RGB | Output RGBA instead of RGB. |
| `with_thread_pool(pool)` | global | Route parallel work to a custom rayon pool. |

---

## 2. Image Analysis

```rust
use astroimage::ImageAnalyzer;

let result = ImageAnalyzer::new()
    .with_max_stars(500)
    .with_thread_pool(pool.clone())
    .analyze("light.fits")?;

// Defaults: two-pass Moffat calibration, mesh-grid background (auto cell size),
// MAD noise, detection sigma 5.0, measure cap 500
```

### Key `AnalysisResult` Fields

| Field | Type | Description |
|-------|------|-------------|
| `width`, `height` | `usize` | Image dimensions (native resolution). |
| `stars` | `Vec<StarMetrics>` | Per-star metrics, sorted by flux descending. |
| `stars_detected` | `usize` | Total stars detected (before measure cap). |
| `median_fwhm` | `f32` | Median FWHM (pixels). Good seeing: 2-3 px. |
| `median_eccentricity` | `f32` | Median eccentricity. 0 = round, >0.5 = elongated. |
| `median_snr` | `f32` | Median per-star SNR. |
| `median_hfr` | `f32` | Median half-flux radius (pixels). |
| `measured_fwhm_kernel` | `f32` | FWHM used for final matched filter kernel. |
| `snr_weight` | `f32` | SNR weight for stacking: median(star_flux)² / (noise² × background). Star-based, gradient-immune. |
| `psf_signal` | `f32` | PSF signal strength: median(star_peaks) / noise. |
| `frame_snr` | `f32` | Per-frame SNR: background / noise. For stacking: `stacked_snr = sqrt(sum(frame_snr_i²))`. |
| `trail_r_squared` | `f32` | Rayleigh R² for directional coherence. 0 = no trail, 1 = strong trail. |
| `median_beta` | `Option<f32>` | Median Moffat β. Always `Some(...)` with default settings (Moffat is always-on). |
| `possibly_trailed` | `bool` | Advisory trail flag (Rayleigh angle coherence test). |

### Key `StarMetrics` Fields

| Field | Type | Description |
|-------|------|-------------|
| `x`, `y` | `f32` | Subpixel centroid position (native resolution). |
| `fwhm_x`, `fwhm_y` | `f32` | FWHM along major/minor axes (pixels). |
| `fwhm` | `f32` | Geometric mean FWHM. |
| `eccentricity` | `f32` | 0 = round, approaching 1 = elongated. |
| `theta` | `f32` | PSF position angle (radians, CCW from +X). |
| `snr` | `f32` | Per-star aperture photometry SNR. |
| `hfr` | `f32` | Half-flux radius (pixels). |
| `peak` | `f32` | Background-subtracted peak (ADU). |
| `flux` | `f32` | Total background-subtracted flux (ADU). |
| `beta` | `Option<f32>` | Moffat β shape parameter. `None` if Gaussian fit. |
| `fit_method` | `FitMethod` | Which fitting method produced this measurement. |
| `fit_residual` | `f32` | Normalized fit quality (lower = better). 1.0 for moments fallback. |

### Analyzer Builder Methods

| Method | Default | Description |
|--------|---------|-------------|
| `with_detection_sigma(f32)` | 5.0 | Detection threshold in sigma. |
| `with_min_star_area(usize)` | 5 | Min star area in pixels (filters hot pixels). |
| `with_max_star_area(usize)` | 2000 | Max star area in pixels (filters galaxies/nebulae). |
| `with_saturation_fraction(f32)` | 0.95 | Reject stars above this fraction of 65535. |
| `with_max_stars(usize)` | 200 | Keep only the brightest N stars in returned result. |
| `with_mrs_layers(usize)` | 0 | MRS wavelet noise layers. 0 = fast MAD noise. 4 = MRS wavelet. |
| `without_debayer()` | debayer on | Skip green-channel interpolation for OSC. |
| `with_trail_threshold(f32)` | 0.5 | R² threshold for Rayleigh trail detection. |
| `with_measure_cap(usize)` | 500 | Max stars to PSF-fit. 0 = measure all. |
| `with_fit_max_iter(usize)` | 25 | LM max iterations for measurement pass (calibration always 50). |
| `with_fit_tolerance(f64)` | 1e-4 | LM convergence tolerance for measurement pass (calibration always 1e-6). |
| `with_fit_max_rejects(usize)` | 5 | LM consecutive reject bailout. |
| `with_thread_pool(pool)` | global | Route parallel work to a custom rayon pool. |

### Recommended Settings

The defaults are tuned for general-purpose analysis. Mesh-grid background estimation
is automatic (cell size derived from image dimensions), two-pass Moffat calibration
handles gradients, vignetting, and nebulosity out of the box.

```rust
// Recommended production configuration
let result = ImageAnalyzer::new()
    .with_max_stars(500)             // Return top 500 for annotation
    .with_thread_pool(pool.clone())
    .analyze("light.fits")?;
```

This uses the defaults: two-pass Moffat calibration, mesh-grid background (auto cell
size), MAD noise, detection sigma 5.0, measure cap 500,
fit-residual-weighted statistics. These defaults are well-tested across all image
types (mono, OSC, dense fields, nebulae). Pipeline runs ~250-750ms per 26 MP frame.

**Measure cap:** The `measure_cap` (default 500) limits how many stars undergo
PSF fitting. Stars are sorted by flux (brightest first) before capping. Statistics
(FWHM, eccentricity, HFR) are computed from the measured population, while
`stars_detected` reports the raw detection count. Set to 0 to measure all.

### Accuracy-Tested Settings

The default pipeline produces accurate metrics. Based on file-by-file comparison
of 77 M42 frames (300s + 30s exposures), no special configuration is required:

```rust
let result = ImageAnalyzer::new()
    .with_max_stars(500)
    .with_thread_pool(pool.clone())
    .analyze("light.fits")?;
```

#### Expected accuracy (file-by-file, 77 M42 frames)

| Metric | Accuracy | Notes |
|--------|----------|-------|
| **FWHM** | R²=0.978, slope 0.985 | Near-perfect proportional match |
| **Eccentricity** | R²=0.969, slope 1.121 | Excellent correlation |
| **Noise** | Consistent | MAD noise (default); MRS wavelet optional |
| **Star count** | Conservative | Favors quality over quantity |

#### OSC (Bayer) handling

OSC images are green-interpolated before detection and PSF fitting: non-green
(R/B) pixels are replaced with a distance-weighted average of neighboring green
values. All pixels in the interpolated image are used for PSF fitting — no
green mask. This produces accurate FWHM values for Bayer-pattern data.

#### Eccentricity methodology note

Eccentricity is computed from the **Moffat/Gaussian fit** axis ratio:
`e = sqrt(1 - (b/a)^2)`. Frame-level eccentricity uses a sigma-clipped
weighted median with fit-residual weighting. R²=0.969 indicates excellent
relative ranking agreement.

#### Fit-residual weighting

Statistics (FWHM, eccentricity, HFR) are computed using **sigma-clipped weighted
medians** where each star's weight is `1 / (1 + fit_residual)`. Well-fit stars
(low residual) contribute more than poorly-fit or moments-fallback stars. This
continuous quality weighting is the primary reason for the improved R² values vs v0.7.1.

#### When to use each setting

| Setting | Default | Use when |
|---------|---------|----------|
| `with_mrs_layers(n)` | 0 (MAD) | Default uses fast MAD noise. Set to 4 for MRS wavelet (more robust against nebulosity, ~200ms slower). |
| `with_measure_cap(n)` | 500 | Default is fast and accurate. Increase for very dense fields. Set 0 for no limit. |
| `with_fit_max_iter(n)` | 25 | Increase to 50 for maximum accuracy at the cost of speed. |
| `with_fit_tolerance(tol)` | 1e-4 | Decrease to 1e-6 for tighter convergence (slower). |
| `with_trail_threshold(t)` | 0.5 | Lower to 0.3 for more aggressive trail detection (more false positives). Raise to 0.7 for fewer flags. |

For general-purpose use, the defaults give the best overall results:
accurate FWHM, stable eccentricity, fast pipeline (~250-750ms per frame).

### Analyzing pre-loaded data

```rust
// Planar f32, [0..65535] ADU range. 3-channel: RRRGGGBBB layout.
let result = ImageAnalyzer::new()
    .analyze_data(&pixel_data, width, height, channels)?;
```

---

## 3. Star Annotation Overlay

Three tiers of API — pick the level of control Athenaeum needs:

### Tier 1: `compute_annotations()` — Raw Geometry

Returns `Vec<StarAnnotation>` with transformed coordinates and ellipse parameters.
Use this when Athenaeum renders with its own drawing system (Canvas2D, SwiftUI, etc.).

```rust
use astroimage::{compute_annotations, AnnotationConfig};

let annotations = compute_annotations(
    &result,
    image.width,
    image.height,
    image.flip_vertical,
    &AnnotationConfig::default(),
);

for ann in &annotations {
    // ann.x, ann.y         — centroid in output image coordinates
    // ann.semi_major       — ellipse semi-major axis (output px)
    // ann.semi_minor       — ellipse semi-minor axis (output px)
    // ann.theta            — rotation angle (radians, CCW from +X)
    // ann.eccentricity     — original eccentricity value
    // ann.fwhm             — original FWHM (analysis pixels)
    // ann.color            — [R, G, B] based on color scheme
    my_canvas.draw_ellipse(ann.x, ann.y, ann.semi_major, ann.semi_minor, ann.theta, ann.color);
}
```

### Tier 2: `create_annotation_layer()` — RGBA Overlay

Returns a transparent RGBA buffer (same dimensions as output). Composite over the
base image to toggle annotations on/off without re-rendering.

```rust
use astroimage::create_annotation_layer;

let layer = create_annotation_layer(
    &result,
    image.width,
    image.height,
    image.flip_vertical,
    &config,
);
// layer: Vec<u8>, RGBA, length = width * height * 4
// Transparent where no annotations, opaque colored pixels on ellipses
```

### Tier 3: `annotate_image()` — Burn-In

Draws directly onto a `ProcessedImage`. Simplest path — modifies image in place.

```rust
use astroimage::annotate_image;

annotate_image(&mut image, &result, &AnnotationConfig::default());
ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```

### `AnnotationConfig`

All fields are public — override only what you need:

```rust
use astroimage::{AnnotationConfig, ColorScheme};

let config = AnnotationConfig {
    color_scheme: ColorScheme::Eccentricity,
    show_direction_tick: true,
    min_radius: 6.0,
    max_radius: 60.0,
    line_width: 2,
    ecc_good: 0.5,    // ≤ 0.5 → green (good)
    ecc_warn: 0.6,    // 0.51-0.6 → yellow (warning), > 0.6 → red (problem)
    fwhm_good: 1.3,   // FWHM ratio thresholds (star FWHM / median FWHM)
    fwhm_warn: 2.0,
};
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `color_scheme` | `ColorScheme` | `Eccentricity` | `Eccentricity`, `Fwhm`, or `Uniform`. |
| `show_direction_tick` | `bool` | `true` | Draw ticks along elongation axis (when ecc > 0.15). |
| `min_radius` | `f32` | `6.0` | Min ellipse semi-axis (output px). |
| `max_radius` | `f32` | `60.0` | Max ellipse semi-axis (output px). |
| `line_width` | `u8` | `2` | `1` = 1px, `2` = 3px cross, `3` = 5px diamond. |
| `ecc_good` | `f32` | `0.5` | Eccentricity ≤ this → green. |
| `ecc_warn` | `f32` | `0.6` | Eccentricity between good and warn → yellow. Above → red. |
| `fwhm_good` | `f32` | `1.3` | FWHM ratio below this → green. |
| `fwhm_warn` | `f32` | `2.0` | FWHM ratio between good and warn → yellow. Above → red. |

### Coordinate Transform

Star positions from `AnalysisResult` are in native resolution. The output image
may be smaller due to debayer (2x), downscale, and/or preview binning. The
annotation system handles this automatically:

```
scale_x = output_width  / result.width
scale_y = output_height / result.height
x_out = star.x * scale_x
y_out = star.y * scale_y   (or flipped if flip_vertical)
```

| Scenario | Analysis dims | Output dims | Scale |
|----------|---------------|-------------|-------|
| Mono, no downscale | 4096 x 3072 | 4096 x 3072 | 1.0 |
| OSC, debayer only | 4096 x 3072 | 2048 x 1536 | 0.5 |
| OSC, debayer + 2x | 4096 x 3072 | 1024 x 768 | 0.25 |
| Mono, preview bin | 4096 x 3072 | 2048 x 1536 | 0.5 |

### What Gets Drawn

For each detected star:
- **Rotated ellipse** using `fwhm_x`, `fwhm_y`, `theta`. Semi-axes scaled 2.5x for visibility, clamped to `[min_radius, max_radius]`.
- **Color** based on scheme: green (good) → yellow (warning) → red (problem).
- **Direction ticks** extending from ellipse edge along major axis. Length proportional to eccentricity. Shows drift direction, coma pattern, or curvature.

---

## 4. Trail Detection (Rayleigh Test)

Trail detection is **advisory** — the pipeline always computes full metrics and
reports trail status via `trail_r_squared` and `possibly_trailed`.

Trail detection uses the Rayleigh test on PSF-fit stars (angle coherence).
Requires at least 20 measured stars and FWHM ≥ 2.0 px. Includes optical
aberration suppression: if the elongation pattern is radial (coma/field curvature)
or eccentricity increases with distance from center (tilt), the trail flag is
suppressed.

### `trail_r_squared`

| Value | Meaning |
|-------|---------|
| 0.0 | Uniform angles (no trail), < 20 stars, or FWHM < 2.0 px |
| 0.01-0.15 | Typical undersampled round stars (grid-induced) |
| 0.15-0.40 | Oversampled with field-angle effects (coma, curvature) |
| 0.50-0.70 | Suspicious — possible mild trailing |
| 0.70-1.00 | Strong coherent trail |

### `possibly_trailed`

Set to `true` when the Rayleigh test fires:

- **Path A**: R̄² > threshold AND p < 0.01 (catches coherent trails: RA drift)
- **Path B**: R̄² > 0.05 AND median detection-stage ecc > 0.6 AND p < 0.05 (catches undersampled trails)
- Requires ≥20 detected stars. The threshold defaults to 0.5, configurable via `with_trail_threshold()`.

### Effect on Statistics

When `possibly_trailed = true`, the ecc ≤ 0.8 filter is **bypassed** for FWHM,
eccentricity, and HFR statistics. This ensures trailed frames report accurate (high)
eccentricity values rather than having the signal suppressed by the filter.

---

## 5. Thread Pool

Share a single pool across all pipelines to prevent oversubscription:

```rust
use std::sync::Arc;
use astroimage::{ImageConverter, ImageAnalyzer, ThreadPoolBuilder};

let pool = Arc::new(
    ThreadPoolBuilder::new()
        .num_threads(num_cpus::get())
        .build()?
);

let converter = ImageConverter::new()
    .with_thread_pool(pool.clone());

let analyzer = ImageAnalyzer::new()
    .with_thread_pool(pool.clone());
```

---

## 6. Batch Analysis Performance

### How It Works

`analyze()` is CPU-bound — each call uses rayon internally for parallel PSF fitting,
background estimation, and wavelet transforms. Each call has serial gaps between stages
(background → detection → fitting → stats) where pool threads briefly idle.

On a typical 26-megapixel frame (release build):

| Field type | Stars | Time per frame |
|---|---|---|
| Dense (cocoon) | ~115K | ~4.6s |
| Sparse (mono) | ~5.7K | ~1.4s |
| OSC | ~4.6K | ~1.2s |

### Key Requirement: Shared Thread Pool

All concurrent `analyze()` calls **must** share a single rayon pool via
`with_thread_pool()`. With a shared pool, multiple frames submit work to the same N
threads — rayon's work-stealing handles queuing, no oversubscription occurs. Without a
shared pool, each frame may contend on the global pool without `pool.install()` routing,
losing control over thread allocation.

```rust
use std::sync::Arc;
use astroimage::{ImageAnalyzer, ThreadPoolBuilder};

let pool = Arc::new(
    ThreadPoolBuilder::new()
        .num_threads(num_cpus::get())
        .build()?
);

let analyzer = ImageAnalyzer::new()
    .with_thread_pool(pool.clone());
```

### Recommended Concurrency: buffer_unordered(2)

With a shared pool, moderate concurrency **helps** — while frame A is in a serial
section between stages, frame B's parallel work keeps the pool threads busy.

```rust
// GOOD: 2 concurrent frames, shared pool fills serial gaps
let results = stream::iter(files)
    .map(|path| {
        let analyzer = analyzer.clone();
        tokio::task::spawn_blocking(move || analyzer.analyze(&path))
    })
    .buffer_unordered(2)
    .collect::<Vec<_>>()
    .await;
```

**Why not higher?** Each in-flight frame needs several full-image buffers (~200MB for
26MP). At concurrency 4, that's ~800MB competing for L2/L3 cache, causing thrashing
that offsets the scheduling benefit. Concurrency 2 is the sweet spot: good pool
utilization without excessive memory/cache pressure.

| Concurrency | CPU utilization | Memory per frame | Cache behavior |
|---|---|---|---|
| `buffered(1)` | Threads idle during serial gaps | ~200MB | Clean |
| `buffer_unordered(2)` | Good — gaps filled by second frame | ~400MB | Mild pressure |
| `buffer_unordered(4)` | Marginally better than 2 | ~800MB | Thrashing |

### Migration from buffer_unordered(4)

Change concurrency from 4 to 2 and ensure a shared thread pool:

```rust
// Before: no shared pool, high concurrency
let analyzer = ImageAnalyzer::new();
// ... buffer_unordered(4)

// After: shared pool, concurrency 2
let pool = Arc::new(ThreadPoolBuilder::new().num_threads(num_cpus::get()).build()?);
let analyzer = ImageAnalyzer::new()
    .with_thread_pool(pool.clone());
// ... buffer_unordered(2)
```

### Speed Tuning for Batch

For batch frame ranking where speed matters more than maximum accuracy:

```rust
let analyzer = ImageAnalyzer::new()
    .with_max_stars(50)       // fewer star records returned (stats still use all measured)
    .with_thread_pool(pool.clone());
```

| Optimization | Speedup | Trade-off |
|---|---|---|
| `with_max_stars(50)` | Negligible | Fewer star records in output (medians unchanged) |
| Shared pool + concurrency 2 | ~15-25% vs no pool + concurrency 4 | Requires pool setup |

### Expected Batch Throughput (Release Build)

For 100 frames of 26-megapixel data on an 8-core machine:

| Approach | Dense field | Sparse field |
|---|---|---|
| No shared pool, `buffer_unordered(4)` | ~500-600s | ~180-200s |
| Shared pool, `buffer_unordered(2)` | ~380-420s | ~110-130s |
| Shared pool, `buffer_unordered(2)`, defaults (MAD noise) | ~330-370s | ~95-115s |

---

## 7. Full Integration Example

Complete workflow: convert, analyze, and annotate with toggleable overlay.

```rust
use std::sync::Arc;
use astroimage::{
    ImageConverter, ImageAnalyzer,
    compute_annotations, create_annotation_layer, annotate_image,
    AnnotationConfig, ColorScheme,
    ThreadPoolBuilder,
};

let pool = Arc::new(ThreadPoolBuilder::new().num_threads(8).build()?);

// Configure annotation thresholds (expose in Athenaeum preferences)
let config = AnnotationConfig {
    color_scheme: ColorScheme::Eccentricity,
    ecc_good: 0.5,
    ecc_warn: 0.6,
    line_width: 2,
    ..AnnotationConfig::default()
};

// ── Per-frame processing ──

let image = ImageConverter::new()
    .with_rgba_output()
    .with_thread_pool(pool.clone())
    .process("light_001.fits")?;

// Defaults handle everything: two-pass Moffat calibration,
// mesh-grid background (auto cell size), MAD noise (default),
// measure cap 500, fit-residual-weighted statistics, Rayleigh trail detection
let result = ImageAnalyzer::new()
    .with_max_stars(500)             // Return top 500 stars for annotation
    .with_thread_pool(pool.clone())
    .analyze("light_001.fits")?;

// Display metrics in subframe inspector
println!("Stars: {} (detected: {})  FWHM: {:.2}  Ecc: {:.3}  HFR: {:.2}  R²: {:.3}",
    result.stars.len(),
    result.stars_detected,
    result.median_fwhm,
    result.median_eccentricity,
    result.median_hfr,
    result.trail_r_squared,
);
if result.possibly_trailed {
    println!("  WARNING: possibly trailed (R²={:.3})", result.trail_r_squared);
}
if let Some(beta) = result.median_beta {
    println!("  Moffat beta: {:.2}", beta);
}

// Option A: Get raw geometry for native UI rendering
let annotations = compute_annotations(
    &result, image.width, image.height, image.flip_vertical, &config,
);
// → render with SwiftUI/Canvas2D, toggle visibility in UI

// Option B: Get RGBA overlay for layer compositing
let overlay = create_annotation_layer(
    &result, image.width, image.height, image.flip_vertical, &config,
);
// → composite over base image when user toggles "show stars"

// Option C: Burn-in for export/thumbnail
let mut export_image = ImageConverter::new()
    .with_thread_pool(pool.clone())
    .process("light_001.fits")?;
annotate_image(&mut export_image, &result, &config);
ImageConverter::save_processed(&export_image, "annotated.jpg", 95)?;
```

### Accuracy-tested variant

The defaults already produce accurate metrics. No special configuration
is needed. To increase wavelet layers for nebula-rich fields:

```rust
let result = ImageAnalyzer::new()
    .with_mrs_layers(4)              // MRS wavelet for nebula-rich fields
    .with_max_stars(500)
    .with_thread_pool(pool.clone())
    .analyze("light_001.fits")?;

// FWHM: R²=0.978, slope 0.985 (fit-residual-weighted)
// Eccentricity: R²=0.969, slope 1.121
// Trail detection: Rayleigh angle coherence test
```

---

## 8. Integration Recommendations

### Subframe inspector columns

Core metrics — always available:

| Column | Field | Format | Notes |
|--------|-------|--------|-------|
| FWHM | `median_fwhm` | `{:.2}` | Median FWHM in pixels. Good seeing: 2-3 px |
| Eccentricity | `median_eccentricity` | `{:.3}` | Median eccentricity. 0=round, >0.5=elongated |
| HFR | `median_hfr` | `{:.2}` | Half-flux radius in pixels |
| Stars | `stars.len()` | integer | Stars returned (capped at `max_stars`) |
| Stars Detected | `stars_detected` | integer | Total stars detected (before measure cap) |
| Moffat β | `median_beta` | `{:.2}` | Moffat shape param. Typical: 2-5. Always present with defaults |
| Background | `background` | `{:.0}` | Global background level (ADU) |
| Noise | `noise` | `{:.1}` | Background noise sigma (ADU) |

Additional quality indicators:

| Column | Field | Format | Notes |
|--------|-------|--------|-------|
| Trail R² | `trail_r_squared` | `{:.3}` | 0 = no trail, >0.5 = suspicious, >0.7 = strong |
| Trailed | `possibly_trailed` | icon/badge | Warning, not auto-reject |
| PSF Signal | `psf_signal` | `{:.0}` | median(star_peaks) / noise |
| SNR Weight | `snr_weight` | `{:.1}` | SNR weight for stacking: median(star_flux)² / (noise² × bg). Star-based, gradient-immune |
| Frame SNR | `frame_snr` | `{:.1}` | background / noise, for stacking prediction |

### Annotation toggle in image viewer

- Use Tier 1 (`compute_annotations`) or Tier 2 (`create_annotation_layer`) for
  toggleable overlays. Do not re-render the base image when toggling.
- Cache `Vec<StarAnnotation>` per frame — recompute only when the frame or
  config changes.
- Consider exposing `ecc_good`/`ecc_warn` thresholds as user preferences
  (different optical systems have different baselines).

### Annotation color scheme selector

Let users switch between `Eccentricity`, `Fwhm`, and `Uniform` in the viewer.
Each reveals different quality issues:

- **Eccentricity** — tracking, tilt, coma, curvature
- **FWHM** — focus, field curvature, seeing variation
- **Uniform** — just show star positions

### Trail detection

- Display `possibly_trailed` as a warning icon, not auto-reject
- Allow sorting/filtering by `trail_r_squared`
- Consider exposing trail threshold in preferences (default 0.5, range 0.3-0.8)
- Rayleigh test catches coherent trailing (RA drift, satellite trails)
- When `possibly_trailed = true`, eccentricity values are accurate (not suppressed by filter)

---

## 9. Public Exports

All types and functions are re-exported from the crate root:

```rust
use astroimage::{
    // Conversion
    ImageConverter,
    ProcessedImage,

    // Analysis
    ImageAnalyzer,
    AnalysisResult,
    StarMetrics,
    AnalysisConfig,
    FitMethod,

    // Annotation
    annotate_image,
    compute_annotations,
    create_annotation_layer,
    AnnotationConfig,
    ColorScheme,
    StarAnnotation,

    // Thread pool (re-exported from rayon)
    ThreadPool,
    ThreadPoolBuilder,
};
```

---

## 10. Migrating from v0.6.x to v0.7.0

This section covers every change needed in the Athenaeum codebase to upgrade from
rustafits v0.6.4 to v0.7.0.

### 9.0 Release Summary — What Changed

**Breaking changes:**
- `AnalysisResult.snr_db` field removed
- 9 `ImageAnalyzer` builder methods removed: `with_max_measure`, `without_gaussian_fit`,
  `with_background_mesh`, `with_moffat_fit` / `without_moffat_fit`,
  `with_iterative_background`, `with_mrs_noise`, `with_moffat_beta`, `with_max_distortion`
- `with_mrs_noise(n)` renamed to `with_mrs_layers(n)`

**New defaults (no configuration needed):**
- Two-pass Moffat calibration (free-beta pass 1 → fixed-beta pass 2)
- Mesh-grid background with auto cell size: `max(16, max(w,h)/32)`
- MAD noise estimation (default; MRS wavelet optional via with_mrs_layers)
- Moffat → Gaussian → Moments fallback chain

**New fields:**
- `AnalysisResult.measured_fwhm_kernel: f32` — FWHM used for final matched-filter kernel
- `StarMetrics.fit_method: FitMethod` — which PSF model produced the measurement
  (`FreeMoffat`, `FixedMoffat`, `Gaussian`, `Moments`)

**Behavioral changes:**
- `median_beta` is now always populated (`Some(...)`) — Moffat fitting can't be disabled
- Default noise estimation is MAD; MRS wavelet available via with_mrs_layers(4)
- Default background changed from global to mesh-grid (auto cell size)
- Accuracy improved: FWHM -0.3% median (was +1.1%), ecc R²=0.984

**Remaining builder methods (unchanged from v0.7.0):**
- `with_detection_sigma(f32)` — default 5.0
- `with_min_star_area(usize)` — default 5
- `with_max_star_area(usize)` — default 2000
- `with_saturation_fraction(f32)` — default 0.95
- `with_max_stars(usize)` — default 200
- `with_mrs_layers(usize)` — default 0 (MAD noise)
- `with_trail_threshold(f32)` — default 0.5
- `without_debayer()` — skip OSC green interpolation
- `with_thread_pool(Arc<ThreadPool>)` — custom rayon pool

**New in v0.7.3:**
- `with_measure_cap(usize)` — default 500, max stars to PSF-fit (0 = all)
- `with_fit_max_iter(usize)` — default 25, LM iterations for measurement pass
- `with_fit_tolerance(f64)` — default 1e-4, LM convergence tolerance
- `with_fit_max_rejects(usize)` — default 5, LM consecutive reject bailout

### 9.1 Dependency Version Bump

| File | Change |
|------|--------|
| `crates/athenaeum-core/Cargo.toml` | `rustafits = "0.6.4"` → `rustafits = "0.7"` |
| `crates/athenaeum-tauri/Cargo.toml` | `rustafits = "0.6.4"` → `rustafits = "0.7"` |

### 9.2 Removed Builder Methods

**File:** `crates/athenaeum-core/src/analysis/analyzer.rs` — `build_analyzer()`

9 method calls must be removed from `build_analyzer()`:

| Line | Removed Call | Action |
|------|-------------|--------|
| 19 | `.with_max_measure(config.max_stars)` | Delete — redundant with `with_max_stars` |
| 23-25 | `if !config.use_gaussian_fit { analyzer.without_gaussian_fit() }` | Delete block — Gaussian is always fallback now |
| 27-29 | `if let Some(mesh_size) = config.background_mesh_size { analyzer.with_background_mesh(...) }` | Delete block — mesh-grid is always-on with auto cell size |
| 31-35 | `if config.use_moffat_fit { ... with_moffat_fit() } else { ... without_moffat_fit() }` | Delete block — Moffat is always-on (two-pass calibration) |
| 37-39 | `if config.iterative_background > 0 { analyzer.with_iterative_background(...) }` | Delete block — source-mask re-estimation is always-on |
| 41-43 | `if config.mrs_noise > 0 { analyzer.with_mrs_noise(...) }` | Replace with `.with_mrs_layers(config.mrs_layers as usize)` |
| 45-47 | `if let Some(beta) = config.moffat_beta { analyzer.with_moffat_beta(beta) }` | Delete block — beta auto-derived from calibration pass |
| 49-51 | `if let Some(max_dist) = config.max_distortion { analyzer.with_max_distortion(...) }` | Delete block — distortion filter removed |

After cleanup, `build_analyzer()` should be:

```rust
let mut analyzer = ImageAnalyzer::new()
    .with_detection_sigma(config.detection_sigma as f32)
    .with_min_star_area(config.min_star_area as usize)
    .with_max_star_area(config.max_star_area as usize)
    .with_saturation_fraction(config.saturation_fraction as f32)
    .with_max_stars(config.max_stars as usize)
    .with_trail_threshold(config.trail_threshold as f32)
    .with_mrs_layers(config.mrs_layers as usize);

if let Some(pool) = thread_pool {
    analyzer = analyzer.with_thread_pool(pool);
}
analyzer
```

### 9.3 Removed `snr_db` Field

`AnalysisResult.snr_db` was removed in v0.7.0. All references must be deleted or replaced.

**Rust files:**

| File | Line | Change |
|------|------|--------|
| `crates/athenaeum-core/src/analysis/analyzer.rs` | 79 | Remove `snr_db: result.snr_db as f64,` from `FrameAnalysis` construction |
| `crates/athenaeum-core/src/models.rs` | 719 | Remove `pub snr_db: f64,` from `FrameAnalysis` struct |
| `crates/athenaeum-core/src/rustafits_processor/mod.rs` | 161 | Remove `pub snr_db: f32,` from `AnnotationMetrics` |
| `crates/athenaeum-core/src/rustafits_processor/mod.rs` | 291 | Remove `snr_db: result.snr_db,` from `AnnotationMetrics` construction |
| `crates/athenaeum-core/src/db/schema.rs` | 662 | Remove `snr_db REAL NOT NULL,` from `frame_analysis` table DDL |
| `crates/athenaeum-core/src/db/analysis.rs` | all | Remove `snr_db` from all INSERT/SELECT/query-mapping code |

**TypeScript files:**

| File | Change |
|------|--------|
| `src/types/models.ts` | Remove `snr_db: number` from `FrameAnalysis` and `AnnotationMetrics` interfaces |
| `src/components/blink/FrameInfoPanel.tsx` | Remove the `<InfoRow label="SNR (dB)" ...>` line |
| `src/components/calibration/LightsAnalysisTable.tsx` | Remove `snr_db` from `SortField` type, averages calculation, sort switch case, column header, and table cell |
| `src/components/calibration/RejectionThresholdBar.tsx` | Remove `snr_db` from `RejectionThresholds` interface, `THRESHOLD_FIELDS` array, and `EMPTY_THRESHOLDS` |
| `src/components/LightsAnalysisView.tsx` | Remove `snr_db` rejection threshold logic and CSV export column |

**Database migration:** The `frame_analysis` table already has data with `snr_db`. Options:
- Add a migration to `ALTER TABLE frame_analysis DROP COLUMN snr_db` (SQLite 3.35.0+)
- Or keep the column but stop writing to it (simpler, backward-compatible)
- Or recreate the table without the column (safest for older SQLite)

### 9.4 New `measured_fwhm_kernel` Field (Optional)

`AnalysisResult.measured_fwhm_kernel: f32` is new in v0.7.0 — the FWHM used for the
final matched-filter kernel. Athenaeum can optionally expose this:

- Add `measured_fwhm_kernel: f64` to `FrameAnalysis` in `models.rs`
- Map from `result.measured_fwhm_kernel as f64` in `analyze_frame()`
- Add DB column and display column (or skip if not needed for UI)

### 9.5 New `FitMethod` Enum (Optional)

`StarMetrics.fit_method: FitMethod` is new — tracks which fitting method produced each
star's measurement (`FreeMoffat`, `FixedMoffat`, `Gaussian`, `Moments`).

- Add `FitMethod` to the `use astroimage::` import
- Athenaeum doesn't currently expose per-star metrics in the UI, so this is informational
- Could be useful for future per-star debugging or quality breakdown display

### 9.6 `AnalysisConfig` Simplification

**File:** `crates/athenaeum-core/src/analysis/config.rs`

Remove fields that no longer map to builder methods:

| Field | Action |
|-------|--------|
| `use_gaussian_fit: bool` | Remove — Gaussian is always available as fallback |
| `background_mesh_size: Option<u32>` | Remove — always-on with auto cell size |
| `use_moffat_fit: bool` | Remove — always-on (two-pass calibration) |
| `iterative_background: u32` | Remove — always-on (one re-estimation) |
| `moffat_beta: Option<f32>` | Remove — auto-derived from calibration |
| `max_distortion: Option<f32>` | Remove — distortion filter removed |
| `mrs_noise: u32` | Rename to `mrs_layers: u32`, default `0` (MAD noise) |

Remove the corresponding `default_*` helper functions (`default_true`, `default_one`,
`default_mesh_64`, `default_mrs_0`).

Update `validate()` to remove validation for deleted fields and rename `mrs_noise` checks
to `mrs_layers`.

**Serde migration note:** Existing JSON configs stored in the database will have the old
field names. Use `#[serde(default)]` on `mrs_layers` and add `#[serde(alias = "mrs_noise")]`
so old configs deserialize cleanly. Deleted fields are silently ignored during
deserialization (serde's default behavior when `deny_unknown_fields` is not active).

### 9.7 `AnalysisSettingsPanel` UI Simplification

**File:** `src/components/analysis/AnalysisSettingsPanel.tsx`

Remove UI controls for deleted config fields:
- "Use Gaussian Fit" checkbox (lines 197-208)
- "Use Moffat PSF Model" checkbox (lines 209-220)
- "Background Mesh Size" input (lines 221-235)
- "Iterative Background Passes" input (lines 236-246)
- "Fixed Moffat Beta" input (lines 258-272)
- "Max Distortion" input (lines 273-287)

Rename "MRS Noise Passes" to "MRS Wavelet Layers" and update `mrs_noise` → `mrs_layers`.
Update the description to note it defaults to 0 (MAD noise); set to 4 for MRS wavelet.

**File:** `src/types/analysis-config.ts`

Update the `AnalysisConfig` interface and `DEFAULT_ANALYSIS_CONFIG` to match:
- Remove `use_gaussian_fit`, `background_mesh_size`, `use_moffat_fit`,
  `iterative_background`, `moffat_beta`, `max_distortion`
- Rename `mrs_noise` → `mrs_layers`, default `0` (MAD noise)

### 9.8 `median_beta` is Now Always Present

In v0.6.x, `median_beta` was `Option<f32>` that was `None` when Moffat fitting was
disabled. In v0.7.0, Moffat fitting is always-on, so `median_beta` is always `Some(...)`.

The UI already handles `Option` correctly (conditional rendering), so no code change is
strictly required — but the conditional check can be simplified:
- `FrameInfoPanel.tsx`: The `{metrics.median_beta != null && ...}` guard can stay or
  be removed (beta is always present now)
- `LightsAnalysisTable.tsx`: `hasBeta` will always be `true` — the Beta column will
  always show

---

## 11. Migrating from v0.7.1 to v0.7.3

This section covers changes from rustafits v0.7.1 to v0.7.3.

### 11.0 Release Summary — What Changed

**Non-breaking additions:**
- `StarMetrics.fit_residual: f32` — new field, normalized fit quality metric
- `with_measure_cap(usize)` — new builder method (default 2000)
- `with_fit_max_iter(usize)` — new builder method (default 25)
- `with_fit_tolerance(f64)` — new builder method (default 1e-4)
- `with_fit_max_rejects(usize)` — new builder method (default 5)

**Behavioral changes (no API breakage):**
- Trail detection uses Rayleigh angle coherence test
- Minimum star count for Rayleigh test raised from 5 to 20
- Statistics use fit-residual-weighted sigma-clipped medians (was unweighted)
- FWHM/ecc/HFR statistics bypass ecc ≤ 0.8 filter on trailed frames
- `stars_detected` now reports raw detection count (before measure cap)
- FWHM R² improved from 0.972 to 0.995
- Eccentricity R² improved from 0.916 to 0.943

### 11.1 Dependency Version Bump

| File | Change |
|------|--------|
| `crates/athenaeum-core/Cargo.toml` | `rustafits = "0.7"` (no change needed if using `"0.7"`) |

### 11.2 New `StarMetrics.fit_residual` Field

`fit_residual: f32` is new on `StarMetrics`. Lower values = better fit quality.

- Moffat/Gaussian fits: `sqrt(Σ(data-model)² / n_pixels) / amplitude` (typically 0.01-0.50)
- Moments fallback: fixed `1.0`
- **Optional:** expose in per-star UI if Athenaeum shows individual star metrics
- **Optional:** store in DB if star-level quality analysis is desired

### 11.3 New Builder Methods (All Optional)

All new methods have sensible defaults — no changes required unless tuning is desired.

| Method | Default | When to expose in Athenaeum settings |
|--------|---------|--------------------------------------|
| `with_measure_cap(n)` | 2000 | If users analyze very dense fields (>10K stars) and want full coverage |
| `with_fit_max_iter(n)` | 25 | Advanced tuning only — increase for marginal accuracy gain |
| `with_fit_tolerance(tol)` | 1e-4 | Advanced tuning only |
| `with_fit_max_rejects(n)` | 5 | Advanced tuning only |

**Recommended:** Add `measure_cap` to `AnalysisConfig` with default 2000. The other
three are advanced and can be left at defaults unless users report fitting issues.

### 11.4 Trail Detection Changes

`possibly_trailed` now fires in two stages:
1. Rayleigh test on detection-stage angles (unchanged logic, but min stars raised to 20)
2. **New:** post-measurement check — if median PSF-fit eccentricity > 0.55

**UI impact:**
- More frames may be flagged as trailed (wind shake / vibration now caught)
- `trail_r_squared` may be low (< 0.1) on frames flagged only by Stage 2 — this is correct
  (non-coherent issues have random angles, hence low R²)
- Consider showing a tooltip: "High eccentricity detected (guiding issue)" when
  `possibly_trailed = true` but `trail_r_squared < 0.3`

### 11.5 Statistics Accuracy Changes

Frame-level statistics (FWHM, eccentricity, HFR) are now computed using
**fit-residual-weighted sigma-clipped medians**. This is a transparent improvement —
the fields and types are unchanged, only the values are more accurate.

**Impact on stored data:** If Athenaeum stores historical frame analysis results,
re-analyzing existing frames will produce slightly different FWHM/ecc/HFR values.
Consider noting the analysis version in the database to track this.

### 11.6 `AnalysisConfig` Additions

**File:** `crates/athenaeum-core/src/analysis/config.rs`

Add new optional fields:

```rust
pub measure_cap: u32,      // default 2000
```

**File:** `src/types/analysis-config.ts`

```typescript
measure_cap: number;  // default 2000
```

**File:** `src/components/analysis/AnalysisSettingsPanel.tsx`

Optional: Add "Measure Cap" input (range 0-10000, step 100). Description:
"Maximum stars to PSF-fit. 0 = measure all. Default 2000."

---

## 12. Migration from v0.7.5 to v0.8.0

### 12.1 Summary

v0.8.0 optimizes the analysis pipeline for speed and robustness. The API is fully
backward-compatible — no breaking changes. The behavioral changes are:

1. **Default `measure_cap` changed from 2000 to 500.** This is the biggest impact.
2. **New detection-stage filters** reject nebula knots, edge stars, and extended sources.
3. **Moments-based pre-screening** rejects non-stellar candidates before expensive LM fitting.
4. **Spatial grid selection** balances measured stars across the field.

### 12.2 measure_cap Default Change (Action Required)

The library default for `measure_cap` changed from **2000 to 500**.

**If athenaeum uses `ImageAnalyzer::new()` without calling `.with_measure_cap()`:**
The analysis will now measure ~500 stars instead of ~2000. This produces the same
quality aggregate statistics (median FWHM, eccentricity, HFR) but with fewer per-star
entries in `result.stars`. Frame ranking/sorting is unaffected.

**If athenaeum stores per-star catalogs or uses `result.stars` for star matching:**
Update the default in `AnalysisConfig` from 2000 to 500, or explicitly call
`.with_measure_cap(500)` to match the new library default.

**If athenaeum already calls `.with_measure_cap(n)` with an explicit value:**
No change needed — the explicit value overrides the default.

**Recommended action:** Update `AnalysisConfig.measure_cap` default from 2000 to 500.

### 12.3 Detection Filter Changes (No Action Required)

Three new morphological filters run automatically during Pass 2 detection:

| Filter | What it rejects | Threshold |
|--------|----------------|-----------|
| Sharpness | Nebula knots (too diffuse), cosmic rays (too sharp) | [0.1, 0.9] |
| Concentration index | Extended/diffuse sources | [0.15, 0.75] |
| Edge margin | Stars within 2*FWHM of image edge | max(2*FWHM, 8px) |

These filters are always active in Pass 2 and cannot be disabled. They improve
the quality of `result.stars` by removing non-stellar detections that previously
appeared as measured stars with poor fit residuals.

**Impact on `result.stars_detected`:** May be slightly lower than before on frames
with significant nebulosity (nebula knots no longer counted as stars).

**Impact on aggregate statistics:** Generally improved — fewer junk detections
means cleaner medians.

### 12.4 Moments Pre-Screening (No Action Required)

Before expensive LM fitting, each candidate is now screened using cheap image
moments (~500x faster than LM). Candidates are rejected if:

- Moment-FWHM outside [0.7x, 2.5x] the calibrated field FWHM
- Moment eccentricity > 0.8 (bypassed on trailed frames)
- Peak/mean sharpness outside [0.3, 5.0]

Safety fallback: if screening rejects >50% of candidates (and the frame is not
trailed and has >=3 stars), screening is disabled and all candidates are measured.
This prevents over-filtering on unusual frames (extreme defocus, wind shake).

**No API changes.** The screening is internal to `measure_stars`.

### 12.5 Spatial Grid Selection (No Action Required)

When `measure_cap` limits the number of stars measured, the selection now uses a
4x4 spatial grid with round-robin selection instead of pure brightness ranking.
This ensures spatial coverage across the field (prevents all measured stars from
clustering near the field center where optical quality is best).

**Impact:** `result.stars` will contain stars from all parts of the image, not
just the brightest ones near the center. This improves the representativeness of
aggregate statistics (median FWHM now reflects the full field, including edge
aberrations).

### 12.6 Performance Expectations

With the default `measure_cap=500` and moments pre-screening:

| Scenario | v0.7.5 | v0.8.0 | Speedup |
|----------|--------|--------|---------|
| Typical frame (6252x4176, ~7000 detected stars) | ~1.3s analysis | ~0.2-0.3s analysis | ~5x |
| Dense field (>10000 detected stars) | ~2.0s+ analysis | ~0.3-0.4s analysis | ~5-6x |
| Sparse field (<200 detected stars) | ~0.3s analysis | ~0.3s analysis | ~1x (no change) |

**Note:** These are analysis-only times. Total conversion time (FITS read + stretch +
JPEG encode) adds ~1.5s regardless of analysis settings.

### 12.7 Accuracy Impact

Accuracy comparison (77-frame M42 dataset):

| Metric | v0.7.5 | v0.8.0 |
|--------|--------|--------|
| FWHM R² | 0.995 | 0.977 |
| Ecc R² | 0.943 | 0.946 |

Eccentricity agreement improved slightly. FWHM agreement decreased due to the
spatial grid selection including more edge-of-field stars (which have larger FWHM
from optical aberrations). The per-frame FWHM tracking remains very close — the R²
reduction is driven by a few outlier frames.

### 12.8 Athenaeum Checklist

- [ ] Update rustafits dependency to v0.8.0
- [ ] Update `AnalysisConfig.measure_cap` default from 2000 to 500
- [ ] Update `measure_cap` UI control description: "Default 500 (was 2000)"
- [ ] Re-test batch analysis performance (expect ~5x speedup on analysis-heavy workloads)
- [ ] Verify frame ranking still works correctly with new star selection
- [ ] Optional: re-analyze stored frames to update statistics with new filtering
