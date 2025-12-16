# rustafits

High-performance FITS/XISF to JPEG converter for astronomical images with auto-stretch and Bayer debayering.

## Features

- **FITS & XISF Support**: Native readers for both formats (no external FITS library)
- **Auto-Stretch**: Median-based statistical stretching (QuickFits/PixInsight compatible)
- **Bayer Debayering**: Super-pixel 2x2 block averaging (RGGB, BGGR, GBRG, GRBG)
- **Preview Mode**: 2x2 binning for fast previews (~90ms)
- **SIMD Optimized**: ARM NEON + x86_64 SSE2 with automatic detection

## Supported Formats

| Format | Extensions | Data Types |
|--------|-----------|------------|
| FITS | `.fits`, `.fit` | 8/16/32-bit int, 32/64-bit float |
| XISF | `.xisf` | All sample formats, zlib/LZ4/Zstd compression |

## Installation

### From Source (Recommended)

```bash
# Install dependencies
# macOS
brew install lz4 zstd

# Debian/Ubuntu
sudo apt-get install liblz4-dev libzstd-dev zlib1g-dev pkg-config

# Fedora/RHEL
sudo dnf install lz4-devel libzstd-devel zlib-devel

# Build
git clone https://github.com/eg013ra1n/rustafits
cd rustafits
cargo build --release
sudo cp target/release/rustafits /usr/local/bin/
```

### Cargo

```bash
cargo install rustafits
```

### Homebrew (macOS/Linux)

```bash
brew tap eg013ra1n/rustafits
brew install rustafits
```

## Usage

```bash
# Basic conversion
rustafits image.fits output.jpg
rustafits image.xisf output.jpg

# Fast preview (2x2 binning)
rustafits large.fits preview.jpg --preview

# Downscaled preview
rustafits large.fits preview.jpg --downscale 4

# Options
rustafits <input> <output.jpg> [OPTIONS]
  --downscale <N>   Downscale factor (default: 1)
  --quality <Q>     JPEG quality 1-100 (default: 95)
  --no-debayer      Disable Bayer debayering
  --preview         2x2 binning for mono images
  --log             Show detailed information
```

### Rust Library

```rust
use fits_converter::FitsConverter;

FitsConverter::new()
    .with_preview_mode()
    .with_quality(90)
    .convert("input.fits", "output.jpg")?;
```

## Performance

Benchmarks on Apple M4 (6252x4176 16-bit images):

| Mode | Time |
|------|------|
| FITS | ~460ms |
| FITS (preview) | ~130ms |
| XISF (LZ4 compressed) | ~290ms |

## Architecture

```
rustafits/
├── src/
│   ├── lib.rs            # Library entry
│   ├── lib_c.rs          # Rust FFI wrapper
│   └── bin/rustafits.rs  # CLI tool
├── c_src/
│   ├── fits_processor.c  # Processing pipeline (SIMD stretch, debayer)
│   ├── fits_reader.c     # Custom FITS parser
│   ├── xisf_reader.c     # XISF parser with decompression
│   ├── jpeg_writer.c     # stb_image_write JPEG encoder
│   └── base64.c          # Base64 decoder for XISF
└── build.rs              # C compilation
```

**Dependencies**: zlib, LZ4, Zstandard (for XISF compression)

## Troubleshooting

**Build errors about lz4/zstd**: Install compression libraries (see Installation)

**Slow conversion**: Use `--preview` for mono images or `--downscale 2`

**Black/white output**: Run with `--log` to check stretch parameters

## References

- QuickLook.Plugin.FitsViewer (Siyu Zhang, GPL-3.0) - Auto-stretch algorithm
- PixInsight - Screen Transfer Function documentation
- [FITS Standard](https://fits.gsfc.nasa.gov/)
- [XISF Specification](https://pixinsight.com/xisf/)

## License

GPL-3.0 (implements algorithms from QuickLook.Plugin.FitsViewer)
