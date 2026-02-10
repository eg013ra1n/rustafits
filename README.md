# rustafits

High-performance FITS/XISF to JPEG/PNG converter for astronomical images with auto-stretch, Bayer debayering, and SIMD acceleration. Pure Rust — no system dependencies.

## Features

- **FITS & XISF Support**: Native readers for both formats (no external libraries)
- **Auto-Stretch**: Median-based statistical stretching (QuickFits/PixInsight compatible)
- **Bayer Debayering**: Super-pixel 2x2 block averaging (RGGB, BGGR, GBRG, GRBG)
- **Preview Mode**: 2x2 binning for fast previews
- **SIMD Optimized**: SSE2/AVX2 (x86_64) and NEON (aarch64) with automatic detection
- **In-Memory API**: Get raw pixel data without file I/O — ideal for GUI apps

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

# Options
rustafits <input> <output> [OPTIONS]
  --downscale <N>   Downscale factor (default: 1)
  --quality <Q>     JPEG quality 1-100 (default: 95)
  --no-debayer      Disable Bayer debayering
  --preview         2x2 binning for mono images
  --log             Show detailed information
```

## Library Usage

Add to your `Cargo.toml`:

```toml
[dependencies]
rustafits = "0.3"
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

// image.data     - Vec<u8>, interleaved RGB bytes
// image.width    - pixel width
// image.height   - pixel height
// image.is_color - true if debayered/RGB, false if mono (gray replicated to RGB)
```

### Builder methods

| Method | Description |
|--------|-------------|
| `with_downscale(n)` | Downscale by factor n |
| `with_quality(q)` | JPEG quality 1-100 |
| `without_debayer()` | Skip Bayer debayering |
| `with_preview_mode()` | 2x2 binning for fast previews |

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
| Gray to RGB | scalar | pshufb | yes |
| Debayer (f32) | yes | — | yes |

## Architecture

```
rustafits/
├── src/
│   ├── lib.rs              # Library entry + public API
│   ├── types.rs            # Core types (PixelData, ProcessedImage, etc.)
│   ├── converter.rs        # ImageConverter builder
│   ├── pipeline.rs         # Processing pipeline
│   ├── output.rs           # JPEG/PNG file output
│   ├── bin/rustafits.rs    # CLI tool
│   ├── formats/
│   │   ├── mod.rs          # Format dispatch
│   │   ├── fits.rs         # FITS reader
│   │   └── xisf.rs         # XISF reader (zlib/LZ4/Zstd)
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

## References

- QuickLook.Plugin.FitsViewer (Siyu Zhang) — Stretch algorithm reference
- PixInsight — Screen Transfer Function documentation
- [FITS Standard](https://fits.gsfc.nasa.gov/)
- [XISF Specification](https://pixinsight.com/xisf/)

## License

Apache-2.0
