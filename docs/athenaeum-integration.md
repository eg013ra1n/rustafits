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
// MRS wavelet noise (1 layer), detection sigma 5.0
```

### Key `AnalysisResult` Fields

| Field | Type | Description |
|-------|------|-------------|
| `width`, `height` | `usize` | Image dimensions (native resolution). |
| `stars` | `Vec<StarMetrics>` | Per-star metrics, sorted by flux descending. |
| `stars_detected` | `usize` | Total stars with valid PSF measurements. |
| `median_fwhm` | `f32` | Median FWHM (pixels). Good seeing: 2-3 px. |
| `median_eccentricity` | `f32` | Median eccentricity. 0 = round, >0.5 = elongated. |
| `median_snr` | `f32` | Median per-star SNR. |
| `median_hfr` | `f32` | Median half-flux radius (pixels). |
| `measured_fwhm_kernel` | `f32` | FWHM used for final matched filter kernel. |
| `snr_weight` | `f32` | PixInsight-style SNR weight for stacking. |
| `psf_signal` | `f32` | PSF signal strength: median(star_peaks) / noise. |
| `frame_snr` | `f32` | Per-frame SNR: background / noise. For stacking: `stacked_snr = sqrt(sum(frame_snr_i²))`. |
| `trail_r_squared` | `f32` | Rayleigh R² for directional coherence. 0 = no trail, 1 = strong trail. |
| `median_beta` | `Option<f32>` | Median Moffat β. Always `Some(...)` with default settings (Moffat is always-on). |
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
| `beta` | `Option<f32>` | Moffat β shape parameter. `None` if Gaussian fit. |
| `fit_method` | `FitMethod` | Which fitting method produced this measurement. |

### Analyzer Builder Methods

| Method | Default | Description |
|--------|---------|-------------|
| `with_detection_sigma(f32)` | 5.0 | Detection threshold in sigma. |
| `with_min_star_area(usize)` | 5 | Min connected-component area (filters hot pixels). |
| `with_max_star_area(usize)` | 2000 | Max area (filters galaxies/nebulae). |
| `with_saturation_fraction(f32)` | 0.95 | Reject stars above this fraction of 65535. |
| `with_max_stars(usize)` | 200 | Keep only the brightest N stars in returned result. |
| `with_mrs_layers(usize)` | 1 | MRS wavelet noise layers. |
| `without_debayer()` | debayer on | Skip green-channel interpolation for OSC. |
| `with_trail_threshold(f32)` | 0.5 | R² threshold for trail detection. |
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
size), MRS wavelet noise (1 layer), detection sigma 5.0. These defaults are well-tested
across all image types (mono, OSC, dense fields, nebulae).

### PixInsight-Compatible Settings

The default pipeline already produces metrics close to PixInsight SubframeSelector
v1.9.2. Based on file-by-file comparison of 77 M42 frames (300s + 30s exposures),
no special configuration is required:

```rust
let result = ImageAnalyzer::new()
    .with_max_stars(500)
    .with_thread_pool(pool.clone())
    .analyze("light.fits")?;
```

#### Expected accuracy vs PixInsight (file-by-file, 77 frames)

| Metric | Accuracy vs PI | Notes |
|--------|----------------|-------|
| **FWHM** | -0.3% median bias | Near-perfect match |
| **Eccentricity** | -0.067 offset, R²=0.984 | Systematic methodology difference (see below) |
| **Noise** | 0.92x PI | MRS wavelet noise (default 1 layer) |
| **Star count** | ~3x PI | We detect more faint stars |

#### Eccentricity methodology note

Our eccentricity comes from the **Moffat/Gaussian fit** axis ratio: `e = sqrt(1 - (b/a)^2)`.
PixInsight likely uses **image moments** for eccentricity, which gives systematically
higher values (mean offset ~+0.067). The frame-to-frame correlation is excellent
(R²=0.984), meaning relative rankings are correct — only the absolute scale differs.

If Athenaeum displays eccentricity alongside PI values, consider noting this
calibration difference, or applying a linear correction for display purposes.

#### When to use each setting

| Setting | Use when |
|---------|----------|
| `with_mrs_layers(n)` | Increase wavelet layers for nebula-rich fields (e.g., `with_mrs_layers(2)`). Default 1 layer is sufficient for most images. |

For general-purpose Athenaeum use, the defaults give the best overall results:
accurate FWHM, stable eccentricity, and high star counts for field coverage.

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

## 6. Full Integration Example

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
// mesh-grid background (auto cell size), MRS wavelet noise (1 layer)
let result = ImageAnalyzer::new()
    .with_max_stars(500)             // Return top 500 stars for annotation
    .with_thread_pool(pool.clone())
    .analyze("light_001.fits")?;

// Display metrics in subframe inspector
println!("Stars: {}  FWHM: {:.2}  Ecc: {:.3}  HFR: {:.2}  R²: {:.3}",
    result.stars.len(),
    result.median_fwhm,
    result.median_eccentricity,
    result.median_hfr,
    result.trail_r_squared,
);
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

### PixInsight-compatible variant

The defaults already produce metrics close to PixInsight. No special configuration
is needed. To increase wavelet layers for nebula-rich fields:

```rust
let result = ImageAnalyzer::new()
    .with_mrs_layers(2)              // Extra wavelet layer for heavy nebulosity
    .with_max_stars(500)
    .with_thread_pool(pool.clone())
    .analyze("light_001.fits")?;

// FWHM will be within ~0.3% of PI SubframeSelector
// Eccentricity will be ~0.067 lower than PI (methodology difference)
// Noise will be ~0.92x PI (MRS wavelet noise)
// Star count will be ~3x PI (we detect more faint stars)
```

---

## 7. Integration Recommendations

### Subframe inspector columns

Core metrics — always available:

| Column | Field | Format | Notes |
|--------|-------|--------|-------|
| FWHM | `median_fwhm` | `{:.2}` | Median FWHM in pixels. Good seeing: 2-3 px |
| Eccentricity | `median_eccentricity` | `{:.3}` | Median eccentricity. 0=round, >0.5=elongated |
| HFR | `median_hfr` | `{:.2}` | Half-flux radius in pixels |
| Stars | `stars.len()` | integer | Stars returned (capped at `max_stars`) |
| Stars Detected | `stars_detected` | integer | Total stars with valid PSF fits |
| Moffat β | `median_beta` | `{:.2}` | Moffat shape param. Typical: 2-5. Always present with defaults |
| Background | `background` | `{:.0}` | Global background level (ADU) |
| Noise | `noise` | `{:.1}` | Background noise sigma (ADU) |

Additional quality indicators:

| Column | Field | Format | Notes |
|--------|-------|--------|-------|
| Trail R² | `trail_r_squared` | `{:.3}` | 0 = no trail, >0.5 = suspicious, >0.7 = strong |
| Trailed | `possibly_trailed` | icon/badge | Warning, not auto-reject |
| PSF Signal | `psf_signal` | `{:.0}` | median(star_peaks) / noise |
| SNR Weight | `snr_weight` | `{:.1}` | PixInsight-style weight for stacking |
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

---

## 8. Public Exports

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

## 9. Migrating from v0.6.x to v0.7.0

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
- Source-mask background re-estimation (always one pass after calibration)
- MRS wavelet noise (default 1 layer)
- Moffat → Gaussian → Moments fallback chain

**New fields:**
- `AnalysisResult.measured_fwhm_kernel: f32` — FWHM used for final matched-filter kernel
- `StarMetrics.fit_method: FitMethod` — which PSF model produced the measurement
  (`FreeMoffat`, `FixedMoffat`, `Gaussian`, `Moments`)

**Behavioral changes:**
- `median_beta` is now always populated (`Some(...)`) — Moffat fitting can't be disabled
- Default noise estimation changed from MAD to MRS wavelet (1 layer)
- Default background changed from global to mesh-grid (auto cell size)
- Accuracy vs PixInsight improved: FWHM -0.3% median (was +1.1%), ecc R²=0.984

**Remaining builder methods (unchanged):**
- `with_detection_sigma(f32)` — default 5.0
- `with_min_star_area(usize)` — default 5
- `with_max_star_area(usize)` — default 2000
- `with_saturation_fraction(f32)` — default 0.95
- `with_max_stars(usize)` — default 200
- `with_mrs_layers(usize)` — default 1
- `with_trail_threshold(f32)` — default 0.5
- `without_debayer()` — skip OSC green interpolation
- `with_thread_pool(Arc<ThreadPool>)` — custom rayon pool

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
| `mrs_noise: u32` | Rename to `mrs_layers: u32`, change default from `0` to `1` |

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
Update the description to note it defaults to 1 (always-on).

**File:** `src/types/analysis-config.ts`

Update the `AnalysisConfig` interface and `DEFAULT_ANALYSIS_CONFIG` to match:
- Remove `use_gaussian_fit`, `background_mesh_size`, `use_moffat_fit`,
  `iterative_background`, `moffat_beta`, `max_distortion`
- Rename `mrs_noise` → `mrs_layers`, default `1`

### 9.8 `median_beta` is Now Always Present

In v0.6.x, `median_beta` was `Option<f32>` that was `None` when Moffat fitting was
disabled. In v0.7.0, Moffat fitting is always-on, so `median_beta` is always `Some(...)`.

The UI already handles `Option` correctly (conditional rendering), so no code change is
strictly required — but the conditional check can be simplified:
- `FrameInfoPanel.tsx`: The `{metrics.median_beta != null && ...}` guard can stay or
  be removed (beta is always present now)
- `LightsAnalysisTable.tsx`: `hasBeta` will always be `true` — the Beta column will
  always show
