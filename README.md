# rustafits

High-performance FITS to JPEG converter for astronomical images with QuickFits/PixInsight-compatible auto-stretch and Bayer debayering.

## Features

- **Auto-Stretch**: Median-based statistical stretching (robust to outliers)
- **Bayer Debayering**: Super-pixel 2×2 block averaging (RGGB, BGGR, GBRG, GRBG)
- **Preview Mode**: 2×2 binning for mono images (4× fewer pixels, ~90ms processing)
- **Fast**: 60-100ms for 4096×4096 images (Quickselect + SIMD + compiler optimizations)
- **Multi-Platform SIMD**: ARM NEON + x86_64 SSE2 with automatic selection

## Installation

### From crates.io (Recommended - All Platforms)

```bash
cargo install rustafits
```

### Homebrew (macOS/Linux)

```bash
brew tap eg013ra1n/rustafits
brew install rustafits
```

### Debian/Ubuntu (x86_64)

```bash
# Download .deb from releases
wget https://github.com/eg013ra1n/rustafits/releases/latest/download/rustafits-VERSION-amd64.deb
sudo dpkg -i rustafits-VERSION-amd64.deb
```

### Raspberry Pi / ARM64 Linux

```bash
# Install dependencies
sudo apt-get update
sudo apt-get install -y libcfitsio-dev pkg-config build-essential

# Build from source
cargo install rustafits
```

### Arch Linux (AUR)

```bash
yay -S rustafits
# or
paru -S rustafits
```

### Fedora/RHEL/CentOS

```bash
# Download .rpm or add COPR repository (coming soon)
sudo dnf install rustafits-VERSION.rpm
```

### From Source

**Prerequisites**:
```bash
# macOS
brew install cfitsio

# Debian/Ubuntu
sudo apt-get install libcfitsio-dev pkg-config build-essential

# Fedora/RHEL
sudo dnf install cfitsio-devel pkg-config

# Rust (1.70+)
# Install from https://rustup.rs/
```

**Build**:
```bash
git clone https://github.com/eg013ra1n/rustafits
cd rustafits
cargo build --release
sudo cp target/release/rustafits /usr/local/bin/
```

## Usage

### CLI

```bash
# Basic conversion
rustafits input.fits output.jpg

# Fast preview (mono with binning)
rustafits large.fits preview.jpg --preview --quality 30

# Quick preview (downscaled)
rustafits large.fits preview.jpg --downscale 4 --quality 85

# Full quality
rustafits image.fits output.jpg --quality 100

# Options
rustafits <input.fits> <output.jpg> [OPTIONS]
  --downscale <N>      Downscale factor (default: 1)
  --quality <Q>        JPEG quality 1-100 (default: 95)
  --no-debayer         Disable Bayer debayering
  --preview            Enable preview mode (2x2 binning for mono)
  --log                Show detailed information
```

### Rust Library

Add to `Cargo.toml`:
```toml
[dependencies]
rustafits = { git = "<repo-url>" }
```

Basic usage:
```rust
use fits_converter::FitsConverter;

fn main() -> anyhow::Result<()> {
    FitsConverter::new()
        .with_preview_mode()  // Fast 2x2 binning for mono
        .with_quality(90)
        .convert("input.fits", "output.jpg")?;
    Ok(())
}
```

## Processing Pipeline

```
FITS Reading → Downscaling (optional) → Debayering (if Bayer) →
Auto-Stretch → Orientation Correction → JPEG Encoding
```

**Key Steps**:
1. **FITS Reading**: CFITSIO loads 16-bit/32-bit data + metadata (BAYERPAT, ROWORDER)
2. **Debayering**: Super-pixel 2×2 → RGB (no interpolation, 3-5× faster than bilinear)
3. **Auto-Stretch**: Median + MADN-based parameter calculation, per-channel transformation
4. **JPEG Encoding**: stb_image_write with configurable quality

## Performance

### Benchmarks

Test: 6252×4176 16-bit images on Apple M4

| Mode | Resolution | Time | Speedup |
|------|-----------|------|---------|
| Mono (full) | 6252×4176 | ~220ms | Baseline |
| **Mono (preview)** | **3126×2088** | **~90ms** | **2.4×** |
| RGB (debayered) | 3124×2088 | ~100ms | 2.2× |

### Optimizations

1. **Algorithmic**: Quickselect O(n) median vs qsort O(n log n) — 20-30% faster
2. **SIMD**: ARM NEON / x86_64 SSE2 vectorization (4 floats/cycle) — 15-25% faster
3. **Compiler**: `-O3 -march=native -ffast-math`
4. **Preview Mode**: 2×2 binning for mono images — 2.4× faster

**Libraries**:
- **CFITSIO**: NASA's FITS I/O library (handles all formats, compression)
- **stb_image_write**: Single-header JPEG encoder (MIT/public domain)

## Architecture

```
Rust CLI (rustafits) → Rust Wrapper (src/lib_c.rs) → C Core (c_src/)
```

**Why C + Rust?**
- C: Performance-critical processing, direct CFITSIO access
- Rust: Safe API, error handling, easy distribution

**Project Structure**:
```
rustafits/
├── src/
│   ├── lib.rs           # Library entry point
│   ├── lib_c.rs         # Rust FFI wrapper
│   └── bin/
│       └── rustafits.rs # CLI tool
├── c_src/
│   ├── fits_processor.c # Main pipeline
│   ├── fits_processor.h # C API
│   └── jpeg_writer.c    # JPEG encoding
├── ffi/
│   ├── debayer.c        # Bayer algorithms
│   └── stretch.c        # Auto-stretch
├── Cargo.toml
└── build.rs             # Compiles C code
```

## Troubleshooting

**"Failed to find cfitsio library"**
```bash
brew install cfitsio  # macOS
export PKG_CONFIG_PATH=/opt/homebrew/lib/pkgconfig:$PKG_CONFIG_PATH
```

**Slow conversion (>500ms)**
- Use `--preview` for mono images
- Use `--downscale 2` or `4` for quick previews

**Black/white output**
- Run with `--log` to see stretch parameters
- Input may be all zeros or very uniform

## References

- **QuickLook.Plugin.FitsViewer**: Original auto-stretch implementation (Siyu Zhang, GPL-3.0)
- **PixInsight**: Screen Transfer Function documentation
- **CFITSIO**: https://heasarc.gsfc.nasa.gov/fitsio/
- **FITS Standard**: https://fits.gsfc.nasa.gov/

## License

Implements algorithms from GPL-3.0 licensed QuickLook.Plugin.FitsViewer.
