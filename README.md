# FITS to JPEG Converter

A fast Rust library and CLI tool for converting FITS (Flexible Image Transport System) astronomical images to JPEG format. Based on the QuickLook.Plugin.FitsViewer implementation, this converter uses CFITSIO for fast FITS file reading and implements advanced image processing techniques.

## Features

- **Fast FITS Reading**: Uses the `fitsio` crate (CFITSIO Rust bindings)
- **Auto-Stretch Algorithm**: Implements the KStars/PixInsight median-based statistical stretching for optimal contrast
- **Bayer Pattern Debayering**: Super-pixel algorithm for RGGB, BGGR, GBRG, and GRBG patterns
- **Configurable Downscaling**: Nearest-neighbor downscaling for faster processing
- **High Performance**: Parallel processing with Rayon
- **Multiple Interfaces**: Rust library, CLI tool, and C FFI

## Installation

### Prerequisites

Make sure you have CFITSIO installed on your system:

**macOS:**
```bash
brew install cfitsio
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install libcfitsio-dev
```

### Building from Source

```bash
cargo build --release
```

The compiled binary will be at `target/release/fits-to-jpeg`

## Usage

### Command Line Interface

Basic conversion:
```bash
fits-to-jpeg input.fits output.jpg
```

With options:
```bash
fits-to-jpeg input.fits output.jpg --downscale 2 --quality 90
```

Available options:
- `--downscale <N>`: Downscale factor (1 = no downscaling, default: 1)
- `--quality <Q>`: JPEG quality 1-100 (default: 95)
- `--no-debayer`: Disable Bayer pattern debayering

### Rust Library

Add to your `Cargo.toml`:
```toml
[dependencies]
fits-converter = { path = "path/to/fits-converter" }
```

Example usage:
```rust
use fits_converter::FitsConverter;

fn main() -> anyhow::Result<()> {
    // Simple conversion
    FitsConverter::new()
        .convert("input.fits", "output.jpg")?;

    // With custom settings
    FitsConverter::new()
        .with_downscale(2)
        .with_quality(90)
        .without_debayer()
        .convert("input.fits", "output.jpg")?;

    Ok(())
}
```

### C FFI

The library provides a C-compatible interface for integration with other languages.

Example C code:
```c
#include "fits_converter.h"

int main() {
    FitsConverterHandle* converter = fits_converter_create();

    fits_converter_set_downscale(converter, 2);
    fits_converter_set_quality(converter, 90);

    int result = fits_converter_convert(
        converter,
        "input.fits",
        "output.jpg"
    );

    fits_converter_destroy(converter);

    return result == FITS_OK ? 0 : 1;
}
```

Build the C library:
```bash
cargo build --release -p fits-converter-ffi
```

## Image Processing Pipeline

The converter applies the following processing steps:

1. **FITS Reading**: Loads image data using CFITSIO
   - Supports u16 and f32 data types
   - Reads BAYERPAT and ROWORDER header keywords

2. **Downscaling** (optional): Nearest-neighbor sampling for speed
   - Reduces dimensions by specified factor
   - Applied before other processing for efficiency

3. **Debayering** (if Bayer pattern detected):
   - Super-pixel algorithm: 2×2 block → 1 RGB pixel
   - 50% resolution reduction but no interpolation artifacts
   - Faster than bilinear/VNG debayering

4. **Auto-Stretch**: Statistical contrast enhancement
   - Samples up to 500k pixels for parameter calculation
   - Computes median and MADN (Median Absolute Deviation Normalized)
   - Non-linear midtones transformation
   - Per-channel independent processing

5. **Orientation Correction**: Handles ROWORDER keyword
   - Flips image vertically if needed

6. **JPEG Encoding**: High-quality JPEG output
   - Configurable quality (default: 95)
   - RGB or grayscale output

## Performance

This implementation matches the speed of QuickLook.Plugin.FitsViewer by:

- **Statistical sampling**: Only processes 500k pixels for stretch parameters
- **Nearest-neighbor downscaling**: Fast but trades quality for speed
- **Super-pixel debayering**: Much faster than interpolation
- **Parallel processing**: Uses Rayon for multi-threaded stretching
- **Efficient memory layout**: Planar image format for cache efficiency

## Supported FITS Formats

- **Bit depths**: 16-bit unsigned (BITPIX=16), 32-bit float (BITPIX=-32)
- **Color spaces**: Mono, RGB (3-channel)
- **Bayer patterns**: RGGB, BGGR, GBRG, GRBG
- **Dimensions**: 2D and 3D images

## Development

Run tests:
```bash
cargo test
```

Build documentation:
```bash
cargo doc --open
```

## Project Structure

```
fits-converter/
├── src/
│   ├── lib.rs           # Library API
│   ├── fits.rs          # FITS file reading
│   ├── downscale.rs     # Downscaling algorithms
│   ├── debayer.rs       # Bayer pattern debayering
│   ├── stretch.rs       # Auto-stretch algorithm
│   ├── processor.rs     # Image processing pipeline
│   ├── output.rs        # JPEG encoding
│   └── bin/
│       └── fits-to-jpeg.rs  # CLI tool
├── ffi/                 # C FFI wrapper
│   ├── src/lib.rs
│   └── fits_converter.h
└── Cargo.toml
```

## References

This implementation is based on:
- [QuickLook.Plugin.FitsViewer](https://github.com/siyu6974/QuickLook.Plugin.FitsViewer) - Original C++/C# implementation
- CFITSIO - NASA's FITS I/O library
- KStars/PixInsight auto-stretch algorithm

## License

This project is provided as-is for astronomical image processing.
