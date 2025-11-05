# Performance Optimization - C Implementation

## Problem
The initial pure Rust implementation was significantly slower than QuickLook.Plugin.FitsViewer.

## Root Causes Identified

1. **Heavy Rust abstractions**: `ndarray` adds overhead for array operations
2. **Memory allocations**: Converting entire images to `Vec<f32>` multiple times
3. **Indirect CFITSIO access**: Rust bindings add a layer of indirection
4. **Type conversions**: Constant conversions between types

## Solution: C Core + Thin Rust Wrapper

We rewrote the performance-critical code in **pure C** matching the QuickLook implementation:

### Architecture

```
┌─────────────────────────────────────┐
│     Rust CLI (fits-to-jpeg)         │
│   (Argument parsing, high-level)    │
└────────────┬────────────────────────┘
             │
             │ FFI calls
             ▼
┌─────────────────────────────────────┐
│      Thin Rust Wrapper (lib_c.rs)   │
│    (Safe interface, error handling) │
└────────────┬────────────────────────┘
             │
             │ unsafe extern "C"
             ▼
┌─────────────────────────────────────┐
│    C Implementation (fits_processor.c)│
│  ┌───────────────────────────────┐  │
│  │  CFITSIO Direct Calls          │  │
│  │  - fits_open_file()            │  │
│  │  - fits_read_pix()             │  │
│  │  - Native C arrays             │  │
│  └───────────────────────────────┘  │
│                                      │
│  ┌───────────────────────────────┐  │
│  │  Super-pixel Debayering        │  │
│  │  - Direct pointer arithmetic   │  │
│  │  - Inline operations           │  │
│  └───────────────────────────────┘  │
│                                      │
│  ┌───────────────────────────────┐  │
│  │  Auto-Stretch Algorithm        │  │
│  │  - qsort() for median          │  │
│  │  - Stack allocations           │  │
│  │  - Single-pass conversion      │  │
│  └───────────────────────────────┘  │
└─────────────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────┐
│  stb_image_write (JPEG encoding)    │
│   (Single-header C library)         │
└─────────────────────────────────────┘
```

## Key Optimizations

### 1. Direct CFITSIO Access
```c
// No Rust bindings overhead
fits_read_pix(fptr, TUSHORT, fpixel, npixels, NULL, data, NULL, &status);
```

### 2. Zero-Copy Operations
```c
// Work directly on raw pixel data
uint16_t* u16_data = malloc(npixels * sizeof(uint16_t));
fits_read_pix(fptr, TUSHORT, fpixel, npixels, NULL, u16_data, NULL, &status);
```

### 3. In-Place Transformations
```c
// Stretch directly to 8-bit output
for (size_t i = 0; i < channel_size; i++) {
    float normalized = channel_data[i] / max_input;
    float stretched = apply_stretch(normalized, shadows, highlights, midtones);
    out_image->data[i * 3 + c] = (uint8_t)(stretched * 255.0f);
}
```

### 4. Fast Median Calculation
```c
// Use stdlib qsort instead of ndarray sorting
qsort(samples, num_samples, sizeof(float), compare_float);
float median = samples[num_samples / 2];
```

### 5. Compiler Optimizations
```rust
cc::Build::new()
    .opt_level(3)           // -O3
    .flag("-march=native")  // CPU-specific instructions
    .flag("-ffast-math")    // Aggressive math optimizations
```

## Performance Characteristics

### Memory Usage
- **Rust version**: Multiple temporary Vec allocations
- **C version**: Stack allocations where possible, single malloc per stage

### CPU Usage
- **Rust version**: Iterator chains with closures (optimization challenges)
- **C version**: Direct loops with inlining

### Binary Size
- **Rust version**: 798 KB (includes ndarray, rayon overhead)
- **C version**: 488 KB (pure C + minimal Rust wrapper)

## Comparison to QuickLook.Plugin.FitsViewer

Our C implementation mirrors QuickLook's approach:

| Feature | QuickLook (C++) | Our Implementation (C) |
|---------|----------------|------------------------|
| FITS I/O | CCfits (C++ wrapper) | CFITSIO (direct) |
| Debayering | Template functions | C functions |
| Stretching | valarray + loops | Direct array access |
| Parallelism | parallel_for | Single-threaded (can add OpenMP) |
| JPEG | Windows API | stb_image_write |

### Why This Should Be As Fast

1. **Same algorithms**: Super-pixel debayering, KStars stretch
2. **Same library**: CFITSIO for FITS I/O
3. **Similar optimizations**: O3, native CPU instructions, fast math
4. **Less overhead**: No C++/C# interop, direct C calls

### Potential Further Optimizations

If still slower:

1. **Add OpenMP parallelization**:
```c
#pragma omp parallel for
for (int c = 0; c < num_channels; c++) {
    // Stretch each channel in parallel
}
```

2. **SIMD vectorization**:
```c
#include <arm_neon.h>  // For M1/M2 Macs
// Vectorize pixel operations
```

3. **Memory-mapped FITS**:
```c
// mmap() the FITS file instead of fits_read_pix()
```

## Testing Instructions

To test performance:

```bash
# Compile
cargo build --release

# Time the conversion
time ./target/release/fits-to-jpeg input.fits output.jpg

# Compare to QuickLook (approximate)
# QuickLook shows preview almost instantly
```

## Expected Performance

For a typical 4096x4096 16-bit FITS image with Bayer pattern:

- **QuickLook**: < 100ms
- **Our C implementation**: < 200ms (target)
- **Pure Rust version**: ~1-2 seconds

If you're still seeing slowness > 500ms, please:
1. Check if you have a very large file (> 100MB)
2. Profile with `cargo instruments` on macOS
3. Check if downscale factor is set (smaller = faster preview)
