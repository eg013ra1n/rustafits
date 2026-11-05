# Star Annotation Overlay

The annotation module draws measured-star ellipses onto converted images, showing
PSF shape, elongation direction, and quality at a glance. Zero external dependencies —
all rasterization uses integer Bresenham algorithms on the existing u8 pixel buffer.

The stars come from a `StarOverlay` (see [The star list](#the-star-list-staroverlay)).
`rustafits convert --annotate` builds one from the same measurement `rustafits analyze`
prints (`star_measure::analysis_mode`), so the ellipses sit on the stars the numbers
were computed on. A library caller builds one the same way, or from its own list.

## Quick Start

### CLI

```bash
# Annotate with defaults (eccentricity color coding, a budget of 500 stars)
rustafits image.fits output.jpg --annotate

# More stars, with logging
rustafits image.fits output.jpg --annotate --max-stars 1000 -v

# Works with all other options
rustafits image.fits output.jpg --annotate --downscale 2 --preview --quality 90
```

### Library (Simplest Path)

```rust
use std::path::Path;
use astroimage::{annotate_image_with, AnnotationConfig, ImageConverter, StarOverlay};
use astroimage::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
use astroimage::star_measure::raw::read_raw_frame;

// 1. Measure the file (see usage.md)
let mut frame = read_raw_frame(Path::new("light.fits"))?;
frame.drop_unmeasured_planes();
let opts = FrameOptions::default();
let pm = measure_raw_frame(&frame, &opts, None);
let summary = FrameSummary::from_measurement(&pm, &opts);

// 2. Convert the same file
let mut image = ImageConverter::new().process("light.fits")?;

// 3. Build the star list and burn it onto the image
let overlay = StarOverlay::from_measurement(frame.w, frame.h, &pm, &summary);
annotate_image_with(&mut image, &overlay, &AnnotationConfig::default());

// 4. Save
ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```

The functions that take an `AnalysisResult` (`annotate_image`, `compute_annotations`,
`create_annotation_layer`) still work. They belong to the legacy analyzer, which is
deprecated since 1.8.0 and removed in 2.0: see [the legacy functions](#the-legacy-functions).

## The star list: `StarOverlay`

A `StarOverlay` is the list of stars to draw and the grid they were measured on. It is
a plain struct with public fields. Its field set is fixed for 1.x, so a caller can build
one from its own measurements.

| Field | Type | Description |
|-------|------|-------------|
| `stars` | `Vec<OverlayStar>` | The stars to draw. |
| `width`, `height` | `usize` | Size of the pixel grid the stars were measured on. A mosaic is measured at its native size, and rendered at half size. |
| `median_fwhm` | `f32` | The frame's median FWHM in source pixels, for `ColorScheme::Fwhm`. 0.0 when the frame has none: every star then reads green. |

`OverlayStar`:

| Field | Type | Description |
|-------|------|-------------|
| `x`, `y` | `f32` | Centroid in source pixels, in the file's row order (not flipped). |
| `fwhm_x` | `f32` | FWHM along the **major** axis, source pixels. |
| `fwhm_y` | `f32` | FWHM along the minor axis, source pixels. |
| `theta` | `f32` | Angle of the major axis, radians, counter-clockwise from +X. |
| `eccentricity` | `f32` | 0 = round, approaching 1 = elongated. |
| `fwhm` | `f32` | Geometric mean FWHM, source pixels. |

Two constructors:

- `StarOverlay::from_measurement(frame_w, frame_h, &pm, &summary)`: one star per fit in
  `pm.stars`, and `summary.median_fwhm`. Pass the measured frame's `RawFrame::w` and
  `RawFrame::h`, also for a mosaic.
- `StarOverlay::from(&result)`: the stars and the median FWHM of a legacy
  `AnalysisResult`, as they are.

## Three-Tier API

The annotation system provides three levels of control. Pick the one that fits
your application:

### Tier 1: `compute_overlay()` — Raw Geometry

Returns a `Vec<StarAnnotation>` with transformed coordinates and ellipse
parameters. It is empty when the overlay has no stars or a zero size. Use this when
your app has its own rendering system (HTML Canvas, SwiftUI, OpenGL, etc.).

```rust
use astroimage::{compute_overlay, AnnotationConfig, ImageConverter};

let image = ImageConverter::new().process("light.fits")?;
// `overlay` is a StarOverlay, built as shown above

let annotations = compute_overlay(
    &overlay,
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
    // ann.fwhm             — original FWHM (source pixels)
    // ann.color            — [R, G, B] based on color scheme

    my_canvas.draw_ellipse(
        ann.x, ann.y,
        ann.semi_major, ann.semi_minor,
        ann.theta,
        ann.color,
    );
}
```

### Tier 2: `create_overlay_layer()` — RGBA Overlay

Returns an RGBA `Vec<u8>` buffer (same dimensions as the output image) with
transparent background and colored ellipses. Use this when you want to toggle
annotations on/off without re-rendering the base image.

```rust
use astroimage::{create_overlay_layer, AnnotationConfig, ImageConverter};

let image = ImageConverter::new().process("light.fits")?;

let layer = create_overlay_layer(
    &overlay,
    image.width,
    image.height,
    image.flip_vertical,
    &AnnotationConfig::default(),
);

// layer is Vec<u8>, RGBA format, length = width * height * 4
// Transparent where no annotations, opaque colored pixels on ellipses
// Composite over base image in your UI when user toggles "show stars"
```

### Tier 3: `annotate_image_with()` — Burn-In

Draws directly onto a `ProcessedImage` buffer. Simplest path — one call,
modifies the image in place. Use for CLI tools or one-shot rendering.

```rust
use astroimage::{annotate_image_with, AnnotationConfig, ImageConverter};

let mut image = ImageConverter::new().process("light.fits")?;

annotate_image_with(&mut image, &overlay, &AnnotationConfig::default());

// image.data now contains the annotations burned in
ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```

### The legacy functions

`compute_annotations()`, `create_annotation_layer()` and `annotate_image()` take an
`AnalysisResult` from the deprecated `ImageAnalyzer::analyze*`. They are wrappers: each
builds `StarOverlay::from(&result)` and calls its twin above, so both draw the same
pixels. Their signatures have not changed.

| Legacy | Current |
|--------|---------|
| `compute_annotations(&result, w, h, flip, &config)` | `compute_overlay(&overlay, w, h, flip, &config)` |
| `create_annotation_layer(&result, w, h, flip, &config)` | `create_overlay_layer(&overlay, w, h, flip, &config)` |
| `annotate_image(&mut image, &result, &config)` | `annotate_image_with(&mut image, &overlay, &config)` |

The wrappers stay as long as `AnalysisResult` does, which is until 2.0.

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
| `ellipse_scale` | `f32` | `1.2` | Semi-axis as a multiple of the measured FWHM, before clamping. Was a hardcoded `2.5`: on oversampled frames that drew semi-axes wide enough for a clean single star to cover its neighbours and read as a blend. |
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

Returned by `compute_overlay()` (and `compute_annotations()`). All coordinates are in
output image space. The struct derives `Debug`, `Clone`, `Copy` and `PartialEq`.

| Field | Type | Description |
|-------|------|-------------|
| `x` | `f32` | Centroid X in output image coordinates. |
| `y` | `f32` | Centroid Y in output image coordinates. |
| `semi_major` | `f32` | Ellipse semi-major axis (output pixels). `fwhm_x * scale_x * ellipse_scale`, clamped to `[min_radius, max_radius]`. |
| `semi_minor` | `f32` | Ellipse semi-minor axis (output pixels). `fwhm_y * scale_y * ellipse_scale`, clamped to `[min_radius, max_radius]`. |
| `theta` | `f32` | Rotation angle in radians, counter-clockwise from +X axis. |
| `eccentricity` | `f32` | Original eccentricity value from the measurement. |
| `fwhm` | `f32` | Original geometric mean FWHM in source pixels. |
| `color` | `[u8; 3]` | RGB color assigned by the color scheme. |

## Coordinate Transform

Star positions in a `StarOverlay` are in the pixels the frame was **measured** on,
which is the file's native size (a Bayer mosaic is measured on its green sites at full
size). The output image may be smaller due to debayer (2x), downscale, and/or preview
binning.

The annotation system handles this automatically:

```
scale_x = output_width  / overlay.width
scale_y = output_height / overlay.height

x_out = star.x * scale_x
y_out = star.y * scale_y          (or flipped if flip_vertical)
theta = -star.theta               (negated if flip_vertical, else unchanged)
```

When `flip_vertical` is true (FITS images are Y-up, output is Y-down), the
Y-axis reflection reverses the handedness of the coordinate system. Angles
measured counter-clockwise from +X become clockwise, so theta is negated to
keep the annotation ellipse and direction ticks aligned with the actual star
elongation in the output image.

`ProcessedImage` carries a `flip_vertical` field set during processing. The
`annotate_image_with()` function reads it automatically. For `compute_overlay()`
and `create_overlay_layer()`, pass `image.flip_vertical` explicitly. The stars of an
overlay are always in file row order: the flip is applied here, never before.

### Why This Works for All Cases

| Scenario | overlay dimensions | output dimensions | scale |
|----------|-------------------|-------------------|-------|
| Mono, no downscale | 4096 × 3072 | 4096 × 3072 | 1.0 |
| OSC, debayer only | 4096 × 3072 (native mosaic) | 2048 × 1536 (debayered) | 0.5 |
| OSC, debayer + 2x downscale | 4096 × 3072 | 1024 × 768 | 0.25 |
| Mono, preview binning | 4096 × 3072 | 2048 × 1536 | 0.5 |

## What Gets Drawn

For each star in the overlay:

- **Rotated ellipse** centered on the star's centroid, using `fwhm_x`, `fwhm_y`,
  and `theta` from the fit. Semi-axes are scaled by `ellipse_scale` (1.2x by default) and
  clamped to `[min_radius, max_radius]`.

- **Color** based on the chosen scheme and configurable thresholds. Default
  (eccentricity): green ≤ 0.5, yellow 0.51–0.6, red > 0.6.

- **Direction ticks** (optional, default on) extending from the ellipse edge
  along the major axis (`theta`). Length is proportional to eccentricity. Only
  drawn when eccentricity > 0.15. Shows tracking drift direction, coma pattern,
  or field curvature orientation.

No text labels — the color coding and ellipse shape convey quality visually.
Numeric values are available from `FrameSummary` and the star list, or from `-v` in the CLI.

## Integration Examples

### Desktop App — Toggleable Layer

```rust
use std::path::Path;
use astroimage::{
    create_overlay_layer, AnnotationConfig, ImageConverter, StarOverlay,
};
use astroimage::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
use astroimage::star_measure::raw::read_raw_frame;

// Measure with a budget of 1000 stars, then render
let mut frame = read_raw_frame(Path::new("light.fits"))?;
frame.drop_unmeasured_planes();
let mut opts = FrameOptions::default();
opts.max_stars = 1000;
let pm = measure_raw_frame(&frame, &opts, None);
let summary = FrameSummary::from_measurement(&pm, &opts);
let stars = StarOverlay::from_measurement(frame.w, frame.h, &pm, &summary);

let image = ImageConverter::new()
    .with_rgba_output()
    .process("light.fits")?;

// Create the layer (transparent RGBA, same size as image)
let config = AnnotationConfig {
    ecc_good: 0.4,
    ecc_warn: 0.55,
    ..AnnotationConfig::default()
};
let layer = create_overlay_layer(
    &stars, image.width, image.height, image.flip_vertical, &config,
);

// In your UI:
// - Base layer: image.data (RGBA)
// - Annotation layer: layer (RGBA, transparent bg)
// - Toggle annotation visibility without re-rendering base
// - Switch color scheme by regenerating the layer with a different config
```

### Web App — Canvas Overlay with JavaScript

```rust
// Rust/Wasm side: compute geometry, send to JS
let annotations = compute_overlay(
    &stars, canvas_width, canvas_height, image.flip_vertical, &config,
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
    annotate_image_with, AnnotationConfig,
    ImageConverter, StarOverlay, ThreadPoolBuilder,
};
use astroimage::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
use astroimage::star_measure::raw::read_raw_frame;

let pool = Arc::new(ThreadPoolBuilder::new().num_threads(8).build()?);
let config = AnnotationConfig::default();
let opts = FrameOptions::default();

for path in &fits_files {
    let mut image = ImageConverter::new()
        .with_thread_pool(pool.clone())
        .process(path)?;

    let mut frame = read_raw_frame(path)?;
    frame.drop_unmeasured_planes();
    let pm = measure_raw_frame(&frame, &opts, Some(&pool));
    let summary = FrameSummary::from_measurement(&pm, &opts);
    let overlay = StarOverlay::from_measurement(frame.w, frame.h, &pm, &summary);

    annotate_image_with(&mut image, &overlay, &config);

    let out = path.with_extension("annotated.jpg");
    ImageConverter::save_processed(&image, &out, 95)?;
}
```

## CLI Reference

```
-a, --annotate       Overlay star ellipses on the output image
--max-stars <N>      Star budget for --annotate: the brightest good stars, 10 to 2000 (default: 500)
-v, --verbose        Per-file details and timings on stderr (--log is an alias)
```

`--annotate` measures the input with `star_measure::analysis_mode`, the measurement
`rustafits analyze` prints and Athenaeum's Analysis tab shows, with its defaults (6.0 σ,
Moffat 4). It then draws the brightest good stars of that measurement, up to
`--max-stars`, onto the converted output. A file that measurement cannot read fails with
an error that names it, even if the plain render could have read it. With `-v` it
prints the measurement's summary (the star count is the number the detector found; at
most `--max-stars` are drawn):

```
$ rustafits tests/mono.fits annotated.jpg -a -v
tests/mono.fits -> annotated.jpg  6252x4176  read 5 ms | analysis 374 ms (4775 stars, FWHM 2.01) | process 10 ms | encode 82 ms | write 2 ms
1 converted in 0.5 s
```

`--max-stars` below 10 or above 2000 is refused with exit code 2. `--measure-cap` is
still accepted so old scripts run; it prints a warning and does nothing (it is hidden from
`--help` and goes in 2.0). `-j` defaults to 2 with `--annotate`: on 26 MP frames more
jobs were no faster and held more memory (about 1.6 GB at two jobs).

## Public Exports

All annotation types and functions are re-exported from the crate root:

```rust
use astroimage::{
    // Functions over a StarOverlay
    annotate_image_with,      // Tier 3: burn into ProcessedImage
    create_overlay_layer,     // Tier 2: RGBA overlay buffer
    compute_overlay,          // Tier 1: raw geometry

    // Legacy twins over an AnalysisResult (deprecated analyzer)
    annotate_image,
    create_annotation_layer,
    compute_annotations,

    // Types
    StarOverlay,              // The star list and the grid it was measured on
    OverlayStar,              // One star of it
    AnnotationConfig,         // Configuration struct
    ColorScheme,              // Eccentricity | Fwhm | Uniform
    StarAnnotation,           // Per-star output geometry
};
```
