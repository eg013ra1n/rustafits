# rustafits

High-performance FITS/XISF to JPEG/PNG converter for astronomical images with auto-stretch, Bayer debayering, and SIMD acceleration.

## Features

- **FITS & XISF Support**: Native readers for both formats (no external libraries)
- **Auto-Stretch**: Median-based statistical stretching (STF-compatible midtones transfer)
- **Bayer Debayering**: Super-pixel 2x2 block averaging (RGGB, BGGR, GBRG, GRBG)
- **Preview Mode**: 2x2 binning for fast previews
- **SIMD Optimized**: SSE2/AVX2 (x86_64) and NEON (aarch64) with automatic detection
- **RGBA Output**: Optional RGBA pixel data for canvas/web display
- **In-Memory API**: Get raw pixel data without file I/O — ideal for GUI apps
- **Image Analysis**: Two-pass Moffat-primary PSF calibration with adaptive moments screening, star detection, FWHM/HFR/eccentricity measurement, SNR computation, auto-tuned mesh-grid background, and MAD noise estimation (optional MRS wavelet)
- **JPEG via libjpeg-turbo**: SIMD-accelerated JPEG encoding (NEON on aarch64, AVX2 on x86_64) via turbojpeg
- **Star Annotation**: Color-coded ellipse overlay showing PSF shape, elongation direction, and quality grading

## Supported Formats

| Format | Extensions | Data Types |
|--------|-----------|------------|
| FITS | `.fits`, `.fit` | 8/16/32-bit int, 32/64-bit float |
| XISF | `.xisf` | All sample formats, zlib/LZ4/Zstd compression |

## Installation

### Cargo (Recommended)

```bash
cargo install rustafits
```

**Build requirements:** cmake (and nasm on x86_64) for the turbojpeg SIMD JPEG encoder.

```bash
# macOS
brew install cmake nasm

# Debian/Ubuntu
sudo apt install cmake nasm

# Arch Linux
sudo pacman -S cmake nasm
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
rustafits = "0.8"
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

### Image analysis

Detect stars, measure PSF shape, and compute image quality metrics:

```rust
use astroimage::ImageAnalyzer;

let result = ImageAnalyzer::new()
    .with_max_stars(500)
    .with_optics(620.0, 3.76)  // focal length mm, pixel size µm → arcsec output
    .analyze("light.fits")?;

println!("Stars: {}  FWHM: {:.2} px ({:.1}\")  Ecc: {:.3}  Seeing: {:.1}\"",
    result.stars_detected, result.median_fwhm,
    result.median_fwhm_arcsec.unwrap_or(0.0),
    result.median_eccentricity,
    result.median_fwhm_arcsec.unwrap_or(0.0));

// Per-stage timing breakdown
let t = &result.stage_timing;
println!("Timing: bg={:.0}ms det={:.0}ms cal={:.0}ms meas={:.0}ms total={:.0}ms",
    t.background_ms, t.detection_pass1_ms, t.calibration_ms,
    t.measurement_ms, t.total_ms);
```

Default configuration uses a two-pass calibration pipeline: pass 1 fits free-beta Moffat
on bright calibration stars to derive the field PSF model (beta, FWHM). Pass 2 applies
fixed-beta Moffat to all detected stars with Gaussian and moments fallbacks. Background
uses parallelized mesh-grid with auto-tuned cell size and MAD noise estimation (MRS wavelet
available via `with_mrs_layers(4)`). OSC/Bayer images are green-interpolated before
detection and PSF fitting.

### Batch analysis

Analyze multiple images in parallel with progress reporting:

```rust
use astroimage::ImageAnalyzer;

let analyzer = ImageAnalyzer::new()
    .with_optics(620.0, 3.76);

let paths: Vec<&str> = vec!["frame001.fits", "frame002.fits", /* ... */];

let results = analyzer.analyze_batch(&paths, 4, |done, total, path| {
    println!("[{}/{}] {}", done, total, path.display());
});

for (path, result) in &results {
    match result {
        Ok(r) => println!("{}: FWHM={:.2}\" ecc={:.3}",
            path.display(),
            r.median_fwhm_arcsec.unwrap_or(0.0),
            r.median_eccentricity),
        Err(e) => eprintln!("{}: {}", path.display(), e),
    }
}
```

The `concurrency` parameter controls how many frames are analyzed simultaneously.
Results are returned in approximate completion order with their paths.

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

### ImageConverter builder methods

| Method | Description |
|--------|-------------|
| `with_downscale(n)` | Downscale by factor n (Bayer images: debayer counts as 2x, extra downscale applied for n > 2) |
| `with_quality(q)` | JPEG quality 1-100 |
| `without_debayer()` | Skip Bayer debayering |
| `with_preview_mode()` | 2x2 binning for fast previews |
| `with_rgba_output()` | Output RGBA instead of RGB (adds alpha=255 channel) |
| `with_thread_pool(pool)` | Use a custom rayon thread pool (see below) |

### ImageAnalyzer builder methods

| Method | Description |
|--------|-------------|
| `with_detection_sigma(f32)` | Detection threshold in sigma above background (default 5.0) |
| `with_min_star_area(usize)` | Minimum star area in stamp (default 5 px) |
| `with_max_star_area(usize)` | Maximum star area in stamp (default 2000 px) |
| `with_saturation_fraction(f32)` | Reject stars above this fraction of 65535 (default 0.95) |
| `with_max_stars(usize)` | Keep only the brightest N stars (default 200) |
| `with_measure_cap(usize)` | Max stars to PSF-fit for statistics (default 500, 0 = all) |
| `with_mrs_layers(usize)` | Noise layers: 0 = fast MAD (default), 1-6 = MRS wavelet |
| `with_trail_threshold(f32)` | R² threshold for Rayleigh trail detection (default 0.5) |
| `with_optics(f64, f64)` | Focal length (mm) + pixel size (µm) → enables arcsec output |
| `without_debayer()` | Skip green-channel interpolation for OSC images |
| `with_thread_pool(pool)` | Use a custom rayon thread pool |

### AnalysisResult fields

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `width`, `height` | usize | pixels | Image dimensions |
| `background` | f32 | ADU | Global background level |
| `noise` | f32 | ADU | Background noise sigma |
| `stars_detected` | usize | — | Total detections (before measure cap) |
| `stars` | Vec\<StarMetrics\> | — | Per-star metrics (brightest N) |
| `median_fwhm` | f32 | pixels | Median FWHM across measured stars |
| `median_fwhm_arcsec` | Option\<f32\> | arcsec | Median FWHM (requires `with_optics`) |
| `median_eccentricity` | f32 | — | 0 = round, →1 = elongated |
| `median_hfr` | f32 | pixels | Median half-flux radius |
| `median_hfr_arcsec` | Option\<f32\> | arcsec | Median HFR (requires `with_optics`) |
| `median_snr` | f32 | — | Median per-star SNR |
| `plate_scale` | Option\<f32\> | arcsec/px | Plate scale (requires `with_optics`) |
| `trail_r_squared` | f32 | — | Rayleigh R̄² for directional coherence |
| `possibly_trailed` | bool | — | True if coherent trailing detected |
| `median_beta` | Option\<f32\> | — | Moffat β (None if Gaussian/moments) |
| `pass1_detections` | usize | — | Pass 1 detection count (before calibration) |
| `calibrated_fwhm` | f32 | pixels | Calibrated field FWHM from Moffat pass |
| `stars_measured` | usize | — | Stars that survived PSF fitting |
| `moffat_count` | usize | — | Moffat fits among measured stars |
| `gaussian_count` | usize | — | Gaussian fits among measured stars |
| `stage_timing` | StageTiming | ms | Per-stage timing breakdown |

### StarMetrics fields

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `x`, `y` | f32 | pixels | Subpixel centroid position |
| `fwhm_x`, `fwhm_y` | f32 | pixels | FWHM along major/minor axis |
| `fwhm` | f32 | pixels | Geometric mean FWHM |
| `fwhm_arcsec` | Option\<f32\> | arcsec | FWHM (requires `with_optics`) |
| `eccentricity` | f32 | — | 0 = round, →1 = elongated |
| `theta` | f32 | radians | Position angle of major axis |
| `hfr` | f32 | pixels | Half-flux radius |
| `hfr_arcsec` | Option\<f32\> | arcsec | HFR (requires `with_optics`) |
| `snr` | f32 | — | Per-star aperture photometry SNR |
| `peak`, `flux` | f32 | ADU | Peak and total flux |
| `beta` | Option\<f32\> | — | Moffat β parameter |
| `fit_method` | FitMethod | — | FreeMoffat, FixedMoffat, Gaussian, or Moments |
| `fit_residual` | f32 | — | Normalized fit quality (lower = better) |

### StageTiming fields

| Field | Type | Description |
|-------|------|-------------|
| `background_ms` | f64 | Background mesh + noise estimation |
| `detection_pass1_ms` | f64 | Pass 1 star detection |
| `calibration_ms` | f64 | Free-beta Moffat calibration |
| `detection_pass2_ms` | f64 | Pass 2 detection with refined kernel |
| `measurement_ms` | f64 | PSF measurement on measured stars |
| `snr_ms` | f64 | Per-star SNR computation |
| `statistics_ms` | f64 | Statistics aggregation |
| `total_ms` | f64 | Total pipeline wall time |

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
| Mono FITS → JPEG | ~107ms |
| OSC FITS → JPEG | ~67ms |
| XISF → JPEG | ~106ms |
| Mono FITS + annotate | ~970ms |
| Analysis only | ~250-750ms |

### SIMD Acceleration

SIMD is used across the processing pipeline with automatic runtime dispatch:

| Operation | SSE2 | AVX2 | NEON |
|-----------|------|------|------|
| Stretch | 4 px/iter | 8 px/iter | 4 px/iter |
| Binning | yes | yes | yes |
| u16 to f32 | yes | yes | yes |
| Gray to RGB | SSSE3 pshufb | AVX2 pshufb | yes |
| Debayer (f32) | yes | — | yes |
| JPEG encode | — | AVX2 DCT | NEON DCT |

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
│   │   ├── background.rs     # Background estimation (global, mesh-grid, MRS wavelet)
│   │   ├── convolution.rs    # Separable convolution + B3-spline smoothing
│   │   ├── detection.rs      # Star detection (DAOFIND + proximity blend rejection)
│   │   ├── fitting.rs        # LM Gaussian & Moffat PSF fitting (free/fixed beta)
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

**Dependencies**: anyhow, flate2 (rust_backend), lz4_flex, ruzstd, image (PNG only), turbojpeg (libjpeg-turbo), quick-xml, base64, rayon

## Troubleshooting

**Slow conversion**: Use `--preview` for mono images or `--downscale 2`

**Black/white output**: Run with `--log` to check stretch parameters

**Downscale + Bayer/OSC**: The super-pixel debayer already halves resolution (2x). A `--downscale 2` on a Bayer image produces debayer-only output with no extra downscale. Use `--downscale 4` or higher for additional reduction beyond debayering.

## References

- [FITS Standard](https://fits.gsfc.nasa.gov/)
- [XISF Specification](https://pixinsight.com/xisf/)
- Stetson, P.B. (1987) — DAOFIND star detection algorithm
- SExtractor — Background estimation methodology

## License

Apache-2.0
