# rustafits

High-performance FITS/XISF to JPEG/PNG converter for astronomical images with auto-stretch, Bayer debayering, and SIMD acceleration. Pure Rust — no system dependencies.

## Features

- **FITS & XISF Support**: Native readers for both formats (no external libraries)
- **Auto-Stretch**: Median-based statistical stretching (PixInsight STF compatible)
- **Bayer Debayering**: Super-pixel 2x2 block averaging (RGGB, BGGR, GBRG, GRBG)
- **Preview Mode**: 2x2 binning for fast previews
- **SIMD Optimized**: SSE2/AVX2 (x86_64) and NEON (aarch64) with automatic detection
- **RGBA Output**: Optional RGBA pixel data for canvas/web display
- **In-Memory API**: Get raw pixel data without file I/O — ideal for GUI apps
- **Image Analysis**: Star detection, FWHM/HFR/eccentricity measurement, and SNR computation (PixInsight-comparable)
- **Star Annotation**: Color-coded ellipse overlay showing PSF shape, elongation direction, and quality grading

## Supported Formats

| Format | Extensions | Data Types |
|--------|-----------|------------|
| FITS | `.fits`, `.fit` | 8/16/32-bit int, 32/64-bit float |
| XISF | `.xisf` | All sample formats, zlib/LZ4/Zstd compression |

## Installation

### Cargo (Recommended)

No system dependencies needed — everything is pure Rust:

```bash
cargo install rustafits
```

### From Source

```bash
git clone https://github.com/eg013ra1n/rustafits
cd rustafits
cargo build --release
sudo cp target/release/rustafits /usr/local/bin/
```

### Homebrew (macOS/Linux)

```bash
brew tap eg013ra1n/rustafits
brew install rustafits
```

## CLI Usage

```bash
# Basic conversion
rustafits image.fits output.jpg
rustafits image.xisf output.png

# Fast preview (2x2 binning)
rustafits large.fits preview.jpg --preview

# Downscaled output
rustafits large.fits preview.jpg --downscale 4

# Star annotation overlay
rustafits image.fits annotated.jpg --annotate --max-stars 500 --log

# Options
rustafits <input> <output> [OPTIONS]
  --downscale <N>   Downscale factor (default: 1)
                    For Bayer/OSC images, the super-pixel debayer
                    inherently halves resolution, so --downscale 2
                    equals debayer only, --downscale 4 = debayer + 2x
                    extra downscale, etc.
  --quality <Q>     JPEG quality 1-100 (default: 95)
  --no-debayer      Disable Bayer debayering
  --preview         2x2 binning for mono images
  --annotate        Overlay star detection ellipses on the output
  --max-stars <N>   Max stars for annotation analysis (default: 200)
  --log             Show detailed information
```

## Library Usage

Add to your `Cargo.toml`:

```toml
[dependencies]
rustafits = "0.4"
```

### File output

```rust
use astroimage::ImageConverter;

ImageConverter::new()
    .with_preview_mode()
    .with_quality(90)
    .convert("input.fits", "output.jpg")?;
```

### In-memory processing

Get raw RGB pixel data without writing to disk — useful for GUI viewers, web backends, and Tauri apps:

```rust
use astroimage::{ImageConverter, ProcessedImage};

let image: ProcessedImage = ImageConverter::new()
    .with_downscale(2)
    .process("input.fits")?;

// image.data     - Vec<u8>, interleaved RGB or RGBA bytes
// image.width    - pixel width
// image.height   - pixel height
// image.channels - 3 (RGB) or 4 (RGBA)
// image.is_color - true if debayered/RGB, false if mono (gray replicated to RGB)
```

### Star annotation overlay

Analyze an image for stars and draw color-coded ellipses showing PSF shape and quality:

```rust
use astroimage::{
    ImageConverter, ImageAnalyzer,
    annotate_image, AnnotationConfig, ColorScheme,
};

let mut image = ImageConverter::new().process("light.fits")?;
let result = ImageAnalyzer::new()
    .with_max_stars(500)
    .analyze("light.fits")?;

// Burn annotations with default settings (eccentricity color coding)
annotate_image(&mut image, &result, &AnnotationConfig::default());

// Or customize thresholds and color scheme
let config = AnnotationConfig {
    color_scheme: ColorScheme::Eccentricity,
    ecc_good: 0.5,   // ≤ 0.5 → green
    ecc_warn: 0.6,   // 0.51–0.6 → yellow, > 0.6 → red
    ..AnnotationConfig::default()
};
annotate_image(&mut image, &result, &config);

ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```

Three API tiers for different integration needs:

| Function | Returns | Use Case |
|----------|---------|----------|
| `compute_annotations()` | `Vec<StarAnnotation>` | Raw geometry for custom rendering (Canvas2D, SwiftUI, SVG) |
| `create_annotation_layer()` | `Vec<u8>` (RGBA) | Transparent overlay for toggleable layer compositing |
| `annotate_image()` | modifies `ProcessedImage` | Burn-in for CLI or one-shot use |

**`compute_annotations(result, width, height, flip_vertical, config)`** — Transforms star positions from analysis coordinates to output image coordinates (handling debayer scaling, downscale, and vertical flip), computes ellipse semi-axes from `fwhm_x`/`fwhm_y`, and assigns colors. Returns `Vec<StarAnnotation>` where each entry contains `x`, `y`, `semi_major`, `semi_minor`, `theta`, `eccentricity`, `fwhm`, and `color` — everything needed to draw the ellipse in any rendering system.

**`create_annotation_layer(result, width, height, flip_vertical, config)`** — Calls `compute_annotations()` internally, then rasterizes all ellipses and direction ticks onto a transparent RGBA buffer (same dimensions as the output image). Use as a compositable layer that can be toggled on/off without re-rendering the base image.

**`annotate_image(image, result, config)`** — Calls `compute_annotations()` internally, then draws directly onto the `ProcessedImage.data` buffer (RGB or RGBA). Reads `image.flip_vertical` automatically. Simplest path — one call, image modified in place.

**`ImageConverter::save_processed(image, path, quality)`** — Saves a `ProcessedImage` to disk as JPEG or PNG. Use after `annotate_image()` or any other post-processing on the pixel buffer.

#### AnnotationConfig fields

| Field | Default | Description |
|-------|---------|-------------|
| `color_scheme` | `Eccentricity` | `Eccentricity` (tracking/optics), `Fwhm` (focus), or `Uniform` (all green) |
| `show_direction_tick` | `true` | Draw ticks along elongation axis (visible when ecc > 0.15) |
| `min_radius` | `6.0` | Minimum ellipse semi-axis in output pixels |
| `max_radius` | `60.0` | Maximum ellipse semi-axis in output pixels |
| `line_width` | `2` | Line thickness: `1` = 1px, `2` = 3px cross, `3` = 5px diamond |
| `ecc_good` | `0.5` | Eccentricity at or below this is green (good) |
| `ecc_warn` | `0.6` | Eccentricity between good and warn is yellow; above is red |
| `fwhm_good` | `1.3` | FWHM ratio (star/median) below this is green |
| `fwhm_warn` | `2.0` | FWHM ratio between good and warn is yellow; above is red |

See [Annotation Documentation](docs/annotation.md) for full API reference, integration examples, and coordinate transform details.

### Builder methods

| Method | Description |
|--------|-------------|
| `with_downscale(n)` | Downscale by factor n (Bayer images: debayer counts as 2x, extra downscale applied for n > 2) |
| `with_quality(q)` | JPEG quality 1-100 |
| `without_debayer()` | Skip Bayer debayering |
| `with_preview_mode()` | 2x2 binning for fast previews |
| `with_rgba_output()` | Output RGBA instead of RGB (adds alpha=255 channel) |
| `with_thread_pool(pool)` | Use a custom rayon thread pool (see below) |

### Multi-image concurrent processing

By default, all parallel work (debayering, stretch, binning, byte conversion) runs on rayon's global thread pool. This works well for single-image processing, but when processing multiple images concurrently from separate threads, they all compete for the same pool — causing thread oversubscription and degraded throughput.

Use `with_thread_pool()` to route all parallel work to a dedicated or shared pool:

```rust
use std::sync::Arc;
use astroimage::{ImageConverter, ThreadPoolBuilder};

// Create a shared pool once at startup
let pool = Arc::new(
    ThreadPoolBuilder::new()
        .num_threads(num_cpus::get())
        .build()
        .unwrap()
);

// Process multiple images concurrently
let handles: Vec<_> = paths.iter().map(|path| {
    let pool = Arc::clone(&pool);
    let path = path.clone();
    std::thread::spawn(move || {
        ImageConverter::new()
            .with_thread_pool(pool)
            .process(&path)
    })
}).collect();

let results: Vec<_> = handles.into_iter()
    .map(|h| h.join().unwrap())
    .collect();
```

**Recommendations by concurrency level:**

| Concurrent images | Strategy |
|-------------------|----------|
| 1-3 | Default global pool is fine |
| 4-8 | Shared pool via `with_thread_pool()` with `num_cpus` threads |
| 8+ | Shared pool + limit concurrency with a semaphore or channel |

**Memory budget:** Each full-resolution image (e.g. 4096x3072 16-bit) uses ~150 MB peak. For 10 concurrent images, budget ~1.5 GB. Use `with_preview_mode()` or `with_downscale()` to reduce memory usage.

## Performance

Benchmarks on Apple M4 (6252x4176 16-bit images):

| Mode | Time |
|------|------|
| FITS | ~460ms |
| FITS (preview) | ~130ms |
| XISF (LZ4 compressed) | ~290ms |

### SIMD Acceleration

SIMD is used across the processing pipeline with automatic runtime dispatch:

| Operation | SSE2 | AVX2 | NEON |
|-----------|------|------|------|
| Stretch | 4 px/iter | 8 px/iter | 4 px/iter |
| Binning | yes | yes | yes |
| u16 to f32 | yes | yes | yes |
| Gray to RGB | SSSE3 pshufb | AVX2 pshufb | yes |
| Debayer (f32) | yes | — | yes |

## Architecture

```
rustafits/
├── src/
│   ├── lib.rs              # Library entry + public API
│   ├── types.rs            # Core types (PixelData, ProcessedImage, etc.)
│   ├── annotate.rs         # Star annotation overlay (3-tier API)
│   ├── converter.rs        # ImageConverter builder
│   ├── pipeline.rs         # Processing pipeline
│   ├── output.rs           # JPEG/PNG file output
│   ├── bin/rustafits.rs    # CLI tool
│   ├── formats/
│   │   ├── mod.rs          # Format dispatch
│   │   ├── fits.rs         # FITS reader
│   │   └── xisf.rs         # XISF reader (zlib/LZ4/Zstd)
│   ├── analysis/
│   │   ├── mod.rs            # Analyzer builder + pipeline orchestration
│   │   ├── background.rs     # Background estimation (global + mesh-grid)
│   │   ├── detection.rs      # Star detection (DAOFIND + CCL)
│   │   ├── fitting.rs        # Levenberg-Marquardt Gaussian fitting
│   │   ├── metrics.rs        # FWHM, eccentricity, HFR measurement
│   │   └── snr.rs            # Per-star and image-wide SNR
│   └── processing/
│       ├── mod.rs           # Processing module
│       ├── stretch.rs       # Auto-stretch (SIMD)
│       ├── debayer.rs       # Bayer debayering (SIMD)
│       ├── binning.rs       # 2x2 binning (SIMD)
│       ├── downscale.rs     # Integer downscaling
│       └── color.rs         # Color conversions (SIMD)
```

**Dependencies** (all pure Rust): anyhow, flate2 (rust_backend), lz4_flex, ruzstd, image, quick-xml, base64

## Troubleshooting

**Slow conversion**: Use `--preview` for mono images or `--downscale 2`

**Black/white output**: Run with `--log` to check stretch parameters

**Downscale + Bayer/OSC**: The super-pixel debayer already halves resolution (2x). A `--downscale 2` on a Bayer image produces debayer-only output with no extra downscale. Use `--downscale 4` or higher for additional reduction beyond debayering.

## References

- PixInsight — Screen Transfer Function documentation
- [FITS Standard](https://fits.gsfc.nasa.gov/)
- [XISF Specification](https://pixinsight.com/xisf/)

## License

Apache-2.0
