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
use rustafits::{ImageConverter, ProcessedImage};

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
use rustafits::ImageAnalyzer;

let result = ImageAnalyzer::new()
    .with_max_stars(500)
    .with_thread_pool(pool.clone())
    .analyze("light.fits")?;
```

### Key `AnalysisResult` Fields

| Field | Type | Description |
|-------|------|-------------|
| `width`, `height` | `usize` | Image dimensions (native resolution). |
| `stars` | `Vec<StarMetrics>` | Per-star metrics, sorted by flux descending. |
| `stars_detected` | `usize` | Total stars found before max_stars cap. |
| `median_fwhm` | `f32` | Median FWHM (pixels). Good seeing: 2-3 px. |
| `median_eccentricity` | `f32` | Median eccentricity. 0 = round, >0.5 = elongated. |
| `median_snr` | `f32` | Median per-star SNR. |
| `median_hfr` | `f32` | Median half-flux radius (pixels). |
| `snr_db` | `f32` | Image-wide SNR in dB (comparable to PixInsight SNRViews). |
| `snr_weight` | `f32` | PixInsight-style SNR weight for stacking. |
| `psf_signal` | `f32` | PSF signal strength: median(star_peaks) / noise. |
| `trail_r_squared` | `f32` | Rayleigh R² for directional coherence. 0 = no trail, 1 = strong trail. |
| `possibly_trailed` | `bool` | Advisory trail flag (dual-path Rayleigh test). |

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

### Analyzer Builder Methods

| Method | Default | Description |
|--------|---------|-------------|
| `with_detection_sigma(f32)` | 5.0 | Detection threshold in sigma. |
| `with_min_star_area(usize)` | 5 | Min connected-component area (filters hot pixels). |
| `with_max_star_area(usize)` | 2000 | Max area (filters galaxies/nebulae). |
| `with_saturation_fraction(f32)` | 0.95 | Reject stars above this fraction of 65535. |
| `with_max_stars(usize)` | 200 | Keep only the brightest N stars. |
| `without_gaussian_fit()` | Gaussian on | Use fast moments instead of Gaussian fit. |
| `with_background_mesh(usize)` | global | Mesh-grid background with given cell size. |
| `without_debayer()` | debayer on | Skip green-channel interpolation for OSC. |
| `with_max_eccentricity(f32)` | 1.0 | Reject stars above this eccentricity. |
| `with_trail_threshold(f32)` | 0.5 | R² threshold for trail detection. |
| `with_thread_pool(pool)` | global | Route parallel work to a custom rayon pool. |

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
use rustafits::{compute_annotations, AnnotationConfig};

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
use rustafits::create_annotation_layer;

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
use rustafits::annotate_image;

annotate_image(&mut image, &result, &AnnotationConfig::default());
ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```

### `AnnotationConfig`

All fields are public — override only what you need:

```rust
use rustafits::{AnnotationConfig, ColorScheme};

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

## 4. Trail Detection

Trail detection is **advisory** — the pipeline always computes full metrics and
reports trail status via `trail_r_squared` and `possibly_trailed`.

### `trail_r_squared`

| Value | Meaning |
|-------|---------|
| 0.0 | Uniform angles (no trail) or < 5 stars |
| 0.01-0.15 | Typical undersampled round stars |
| 0.15-0.40 | Oversampled with field-angle effects (coma, curvature) |
| 0.50-0.70 | Suspicious — possible mild trailing |
| 0.70-1.00 | Strong trail |

### `possibly_trailed`

Set to `true` when either path fires:
- **Path A**: R² > threshold AND p < 0.01
- **Path B**: median detection-stage ecc > 0.6 AND p < 0.05

The threshold defaults to 0.5, configurable via `with_trail_threshold()`.

---

## 5. Thread Pool

Share a single pool across all pipelines to prevent oversubscription:

```rust
use std::sync::Arc;
use rustafits::{ImageConverter, ImageAnalyzer, ThreadPoolBuilder};

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

## 6. Full Integration Example

Complete workflow: convert, analyze, and annotate with toggleable overlay.

```rust
use std::sync::Arc;
use rustafits::{
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

let result = ImageAnalyzer::new()
    .with_max_stars(1000)
    .with_thread_pool(pool.clone())
    .analyze("light_001.fits")?;

// Display metrics in subframe inspector
println!("Stars: {}  FWHM: {:.2}  Ecc: {:.3}  SNR: {:.1} dB  HFR: {:.2}  R²: {:.3}",
    result.stars.len(),
    result.median_fwhm,
    result.median_eccentricity,
    result.snr_db,
    result.median_hfr,
    result.trail_r_squared,
);

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

---

## 7. Integration Recommendations

### Subframe inspector columns

Add these to the subframe table alongside existing FWHM/eccentricity/SNR:

| Column | Field | Format | Notes |
|--------|-------|--------|-------|
| Trail R² | `trail_r_squared` | `{:.3}` | Continuous quality indicator |
| Trailed | `possibly_trailed` | icon/badge | Warning, not auto-reject |
| Stars | `stars.len()` | integer | Detected star count |
| HFR | `median_hfr` | `{:.2}` | Half-flux radius |
| PSF Signal | `psf_signal` | `{:.0}` | Signal strength |
| SNR Weight | `snr_weight` | `{:.1}` | Stacking weight |

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

---

## 8. Public Exports

All types and functions are re-exported from the crate root:

```rust
use rustafits::{
    // Conversion
    ImageConverter,
    ProcessedImage,

    // Analysis
    ImageAnalyzer,
    AnalysisResult,
    StarMetrics,
    AnalysisConfig,

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
