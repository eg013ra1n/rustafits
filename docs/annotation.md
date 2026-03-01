# Star Annotation Overlay

The annotation module draws detected-star ellipses onto converted images, showing
PSF shape, elongation direction, and quality at a glance. Zero external dependencies —
all rasterization uses integer Bresenham algorithms on the existing u8 pixel buffer.

## Quick Start

### CLI

```bash
# Annotate with defaults (eccentricity color coding, 200 stars)
rustafits image.fits output.jpg --annotate

# More stars, with logging
rustafits image.fits output.jpg --annotate --max-stars 5000 --log

# Works with all other options
rustafits image.fits output.jpg --annotate --downscale 2 --preview --quality 90
```

### Library (Simplest Path)

```rust
use astroimage::{ImageConverter, ImageAnalyzer, annotate_image, AnnotationConfig};

// 1. Convert image
let mut image = ImageConverter::new().process("light.fits")?;

// 2. Analyze same file
let result = ImageAnalyzer::new()
    .with_max_stars(500)
    .analyze("light.fits")?;

// 3. Burn annotations onto image
annotate_image(&mut image, &result, &AnnotationConfig::default());

// 4. Save
ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```

## Three-Tier API

The annotation system provides three levels of control. Pick the one that fits
your application:

### Tier 1: `compute_annotations()` — Raw Geometry

Returns a `Vec<StarAnnotation>` with transformed coordinates and ellipse
parameters. Use this when your app has its own rendering system (HTML Canvas,
SwiftUI, OpenGL, etc.).

```rust
use astroimage::{compute_annotations, AnnotationConfig, ImageAnalyzer, ImageConverter};

let image = ImageConverter::new().process("light.fits")?;
let result = ImageAnalyzer::new().analyze("light.fits")?;

let annotations = compute_annotations(
    &result,
    image.width,
    image.height,
    image.flip_vertical,
    &AnnotationConfig::default(),
);

for ann in &annotations {
    // ann.x, ann.y         — centroid in output image coordinates
    // ann.semi_major       — ellipse semi-major axis (output pixels)
    // ann.semi_minor       — ellipse semi-minor axis (output pixels)
    // ann.theta            — rotation angle (radians, CCW from +X)
    // ann.eccentricity     — original eccentricity value
    // ann.fwhm             — original FWHM (analysis pixels)
    // ann.color            — [R, G, B] based on color scheme

    my_canvas.draw_ellipse(
        ann.x, ann.y,
        ann.semi_major, ann.semi_minor,
        ann.theta,
        ann.color,
    );
}
```

### Tier 2: `create_annotation_layer()` — RGBA Overlay

Returns an RGBA `Vec<u8>` buffer (same dimensions as the output image) with
transparent background and colored ellipses. Use this when you want to toggle
annotations on/off without re-rendering the base image.

```rust
use astroimage::{create_annotation_layer, AnnotationConfig, ImageAnalyzer, ImageConverter};

let image = ImageConverter::new().process("light.fits")?;
let result = ImageAnalyzer::new().analyze("light.fits")?;

let layer = create_annotation_layer(
    &result,
    image.width,
    image.height,
    image.flip_vertical,
    &AnnotationConfig::default(),
);

// layer is Vec<u8>, RGBA format, length = width * height * 4
// Transparent where no annotations, opaque colored pixels on ellipses
// Composite over base image in your UI when user toggles "show stars"
```

### Tier 3: `annotate_image()` — Burn-In

Draws directly onto a `ProcessedImage` buffer. Simplest path — one call,
modifies the image in place. Use for CLI tools or one-shot rendering.

```rust
use astroimage::{annotate_image, AnnotationConfig, ImageAnalyzer, ImageConverter};

let mut image = ImageConverter::new().process("light.fits")?;
let result = ImageAnalyzer::new().analyze("light.fits")?;

annotate_image(&mut image, &result, &AnnotationConfig::default());

// image.data now contains the annotations burned in
ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```

## Configuration

All fields on `AnnotationConfig` are public. Override only what you need:

```rust
use astroimage::{AnnotationConfig, ColorScheme};

let config = AnnotationConfig {
    color_scheme: ColorScheme::Eccentricity,
    ecc_good: 0.5,      // ≤ 0.5 → green (good)
    ecc_warn: 0.6,       // 0.51–0.6 → yellow (warning), > 0.6 → red (problem)
    ..AnnotationConfig::default()
};
```

### `AnnotationConfig` Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `color_scheme` | `ColorScheme` | `Eccentricity` | How stars are colored (see below). |
| `show_direction_tick` | `bool` | `true` | Draw ticks along the elongation axis. Only shown when eccentricity > 0.15. Length proportional to eccentricity. |
| `min_radius` | `f32` | `6.0` | Minimum ellipse semi-axis in output pixels. Prevents tiny annotations on undersampled stars. |
| `max_radius` | `f32` | `60.0` | Maximum ellipse semi-axis in output pixels. Prevents oversized annotations on saturated stars. |
| `line_width` | `u8` | `2` | Line thickness. `1` = single pixel, `2` = 3px cross kernel, `3` = 5px diamond. |
| `ecc_good` | `f32` | `0.5` | Eccentricity threshold for green (good). Stars with ecc ≤ this value are green. |
| `ecc_warn` | `f32` | `0.6` | Eccentricity threshold for yellow (warning). Stars with ecc between `ecc_good` and `ecc_warn` are yellow. Above `ecc_warn` is red. |
| `fwhm_good` | `f32` | `1.3` | FWHM ratio threshold for green. Ratio = star FWHM / median FWHM. Below this is green. |
| `fwhm_warn` | `f32` | `2.0` | FWHM ratio threshold for yellow. Between `fwhm_good` and this is yellow. At or above is red. |

### `ColorScheme` Enum

| Variant | Description |
|---------|-------------|
| `Eccentricity` | Color by roundness: green (round, good) → yellow (slightly elongated) → red (problem). Uses `ecc_good` and `ecc_warn` thresholds. Best for diagnosing tracking/tilt/optics. |
| `Fwhm` | Color by focus quality relative to median: green (tight) → yellow (somewhat bloated) → red (very bloated). Uses `fwhm_good` and `fwhm_warn` thresholds. Best for diagnosing focus or field curvature. |
| `Uniform` | All annotations are green. Use when you only care about star positions/shapes, not quality grading. |

### `StarAnnotation` Fields

Returned by `compute_annotations()`. All coordinates are in output image space.

| Field | Type | Description |
|-------|------|-------------|
| `x` | `f32` | Centroid X in output image coordinates. |
| `y` | `f32` | Centroid Y in output image coordinates. |
| `semi_major` | `f32` | Ellipse semi-major axis (output pixels). Derived from `fwhm_x * scale * 2.5`, clamped to `[min_radius, max_radius]`. |
| `semi_minor` | `f32` | Ellipse semi-minor axis (output pixels). Derived from `fwhm_y * scale * 2.5`, clamped to `[min_radius, max_radius]`. |
| `theta` | `f32` | Rotation angle in radians, counter-clockwise from +X axis. |
| `eccentricity` | `f32` | Original eccentricity value from analysis. |
| `fwhm` | `f32` | Original geometric mean FWHM in analysis pixels. |
| `color` | `[u8; 3]` | RGB color assigned by the color scheme. |

## Coordinate Transform

Star positions from `AnalysisResult` are in **native resolution** (the analysis
pipeline uses green-channel interpolation at full size for OSC). The output image
may be smaller due to debayer (2x), downscale, and/or preview binning.

The annotation system handles this automatically:

```
scale_x = output_width  / result.width
scale_y = output_height / result.height

x_out = star.x * scale_x
y_out = star.y * scale_y          (or flipped if flip_vertical)
```

`ProcessedImage` carries a `flip_vertical` field set during processing. The
`annotate_image()` function reads it automatically. For `compute_annotations()`
and `create_annotation_layer()`, pass `image.flip_vertical` explicitly.

### Why This Works for All Cases

| Scenario | result dimensions | output dimensions | scale |
|----------|-------------------|-------------------|-------|
| Mono, no downscale | 4096 × 3072 | 4096 × 3072 | 1.0 |
| OSC, debayer only | 4096 × 3072 (native green) | 2048 × 1536 (debayered) | 0.5 |
| OSC, debayer + 2x downscale | 4096 × 3072 | 1024 × 768 | 0.25 |
| Mono, preview binning | 4096 × 3072 | 2048 × 1536 | 0.5 |

## What Gets Drawn

For each detected star:

- **Rotated ellipse** centered on the star's centroid, using `fwhm_x`, `fwhm_y`,
  and `theta` from the analysis. Semi-axes are scaled by 2.5x for visibility and
  clamped to `[min_radius, max_radius]`.

- **Color** based on the chosen scheme and configurable thresholds. Default
  (eccentricity): green ≤ 0.5, yellow 0.51–0.6, red > 0.6.

- **Direction ticks** (optional, default on) extending from the ellipse edge
  along the major axis (`theta`). Length is proportional to eccentricity. Only
  drawn when eccentricity > 0.15. Shows tracking drift direction, coma pattern,
  or field curvature orientation.

No text labels — the color coding and ellipse shape convey quality visually.
Numeric values are available via `AnalysisResult` fields or `--log` in the CLI.

## Integration Examples

### Athenaeum / Desktop App — Toggleable Layer

```rust
use astroimage::{
    create_annotation_layer, compute_annotations,
    AnnotationConfig, ColorScheme,
    ImageAnalyzer, ImageConverter,
};

// Process and analyze
let image = ImageConverter::new()
    .with_rgba_output()
    .process("light.fits")?;
let result = ImageAnalyzer::new()
    .with_max_stars(1000)
    .analyze("light.fits")?;

// Create overlay (transparent RGBA, same size as image)
let config = AnnotationConfig {
    ecc_good: 0.4,
    ecc_warn: 0.55,
    ..AnnotationConfig::default()
};
let overlay = create_annotation_layer(
    &result, image.width, image.height, image.flip_vertical, &config,
);

// In your UI:
// - Base layer: image.data (RGBA)
// - Annotation layer: overlay (RGBA, transparent bg)
// - Toggle annotation visibility without re-rendering base
// - Switch color scheme by regenerating overlay with different config
```

### Web App — Canvas Overlay with JavaScript

```rust
// Rust/Wasm side: compute geometry, send to JS
let annotations = compute_annotations(
    &result, canvas_width, canvas_height, image.flip_vertical, &config,
);

// Serialize to JSON or pass via wasm-bindgen
// Each StarAnnotation has: x, y, semi_major, semi_minor, theta, color
```

```javascript
// JavaScript side: draw on a separate canvas overlaid on the image
const ctx = overlayCanvas.getContext('2d');
ctx.clearRect(0, 0, canvas.width, canvas.height);

for (const ann of annotations) {
    ctx.beginPath();
    ctx.ellipse(ann.x, ann.y, ann.semi_major, ann.semi_minor, ann.theta, 0, Math.PI * 2);
    ctx.strokeStyle = `rgb(${ann.color[0]}, ${ann.color[1]}, ${ann.color[2]})`;
    ctx.lineWidth = 2;
    ctx.stroke();
}
```

### Custom Thresholds for Different Optical Systems

```rust
// Fast Newtonian — higher baseline eccentricity from coma
let config = AnnotationConfig {
    ecc_good: 0.6,
    ecc_warn: 0.75,
    ..AnnotationConfig::default()
};

// Refractor — tight tolerances
let config = AnnotationConfig {
    ecc_good: 0.3,
    ecc_warn: 0.5,
    ..AnnotationConfig::default()
};

// Focus quality check — color by FWHM instead of eccentricity
let config = AnnotationConfig {
    color_scheme: ColorScheme::Fwhm,
    fwhm_good: 1.2,
    fwhm_warn: 1.8,
    ..AnnotationConfig::default()
};
```

### Batch Processing with Shared Thread Pool

```rust
use std::sync::Arc;
use astroimage::{
    annotate_image, AnnotationConfig,
    ImageAnalyzer, ImageConverter, ThreadPoolBuilder,
};

let pool = Arc::new(ThreadPoolBuilder::new().num_threads(8).build()?);
let config = AnnotationConfig::default();

for path in &fits_files {
    let mut image = ImageConverter::new()
        .with_thread_pool(pool.clone())
        .process(path)?;

    let result = ImageAnalyzer::new()
        .with_thread_pool(pool.clone())
        .with_max_stars(500)
        .analyze(path)?;

    annotate_image(&mut image, &result, &config);

    let out = path.with_extension("annotated.jpg");
    ImageConverter::save_processed(&image, &out, 95)?;
}
```

## CLI Reference

```
--annotate           Overlay star detection ellipses on the output image
--max-stars <N>      Max stars for annotation analysis (default: 200)
--log                Show analysis stats (star count, median FWHM, eccentricity)
```

`--annotate` runs the full analysis pipeline (`ImageAnalyzer`) on the input file,
then draws annotations onto the converted output. Combined with `--log`, it prints
summary statistics:

```
$ rustafits image.fits output.jpg --annotate --max-stars 1000 --log
Converting to JPEG...
  Input:  image.fits
  Output: output.jpg
  Downscale: 1x
  Quality: 95
  Debayer: enabled
  Preview mode: disabled
Analysis: 784 stars detected, median FWHM=2.60, median ecc=0.430
Conversion successful!
```

## Public Exports

All annotation types and functions are re-exported from the crate root:

```rust
use astroimage::{
    // Functions
    annotate_image,           // Tier 3: burn into ProcessedImage
    create_annotation_layer,  // Tier 2: RGBA overlay buffer
    compute_annotations,      // Tier 1: raw geometry

    // Types
    AnnotationConfig,         // Configuration struct
    ColorScheme,              // Eccentricity | Fwhm | Uniform
    StarAnnotation,           // Per-star output geometry
};
```
