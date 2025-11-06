# Performance Optimization

## Current Performance

**Test System**: Apple M4, 16 GB RAM, macOS

| Mode | Resolution | Time | Speedup |
|------|-----------|------|---------|
| Mono (full) | 6252×4176 | ~220ms | Baseline |
| **Mono (preview)** | **3126×2088** | **~90ms** | **2.4×** |
| RGB (debayered) | 3124×2088 | ~100ms | 2.2× |

## Architecture

```
Rust CLI (rustafits) → Rust Wrapper (src/lib_c.rs) → C Core (c_src/)
```

**Why C + Rust?**
- **C**: Performance-critical processing, direct CFITSIO access, zero-copy operations
- **Rust**: Safe API, error handling, builder pattern, easy distribution

## Optimization Layers

### 1. Algorithmic Optimizations

**Quickselect for Median (O(n) vs O(n log n))**:
```c
// Replaced qsort() with quickselect algorithm
float find_median(float* data, size_t len) {
    return quickselect(data, 0, len - 1, len / 2);
}
```
- **Impact**: 20-30% faster statistical parameter computation
- **Implementation**: Median-of-three pivot selection for robustness

**Statistical Sampling**:
- Max 500k pixels for stretch parameters
- Uniform distribution across image
- No quality loss

**Super-Pixel Debayering**:
- Direct 2×2 → RGB mapping (no interpolation)
- 3-5× faster than bilinear debayering
- Artifact-free output

**Preview Mode (Mono Images)**:
```c
// 2×2 binning reduces pixels by 4×
static void bin_2x2_float(const float* in, float* out, size_t w, size_t h) {
    size_t out_w = w / 2;
    size_t out_h = h / 2;
    for (size_t y = 0; y < out_h; y++) {
        for (size_t x = 0; x < out_w; x++) {
            // Average 2×2 block
            out[y * out_w + x] = (p00 + p01 + p10 + p11) * 0.25f;
        }
    }
}
```
- **Impact**: 2.4× faster for mono images
- **Usage**: `--preview` flag or `.with_preview_mode()`

### 2. SIMD Vectorization (Multi-Platform)

**ARM NEON (Apple Silicon, Raspberry Pi)**:
```c
#ifdef HAS_SIMD_NEON
#include <arm_neon.h>

// Process 4 floats per cycle
float32x4_t input = vld1q_f32(&data[i]);
float32x4_t stretched = apply_stretch_vector(input);
uint32x4_t output = vcvtq_u32_f32(vminq_f32(vmaxq_f32(stretched, zero), max255));
#endif
```

**x86_64 SSE2 (Intel/AMD)**:
```c
#ifdef HAS_SIMD_SSE2
#include <emmintrin.h>

// Process 4 floats per cycle
__m128 input = _mm_loadu_ps(&data[i]);
__m128 stretched = apply_stretch_vector_sse2(input);
__m128i output = _mm_cvtps_epi32(_mm_max_ps(_mm_min_ps(stretched, max255), zero));
#endif
```

**Automatic Platform Detection**:
```c
#if defined(__ARM_NEON) || defined(__aarch64__)
    #define HAS_SIMD_NEON 1
#elif defined(__x86_64__) || defined(_M_X64)
    #define HAS_SIMD_SSE2 1
#endif
```
- **Impact**: 15-25% faster stretch application
- **Zero runtime overhead**: Compile-time selection
- **Fallback**: Optimized scalar code for other platforms

### 3. Compiler Optimizations

```rust
// build.rs
cc::Build::new()
    .opt_level(3)           // -O3: Maximum optimization
    .flag("-march=native")  // Enable NEON/SSE2/AVX
    .flag("-ffast-math")    // Relaxed floating-point semantics
```

**Flags**:
- **-O3**: Auto-vectorization, loop unrolling, inlining
- **-march=native**: All CPU-specific instructions (NEON on M4, AVX on Intel)
- **-ffast-math**: Assumes no NaN/Inf, allows reassociation

### 4. Memory Optimizations

**Cache-Friendly Layout**:
```c
// Planar RGB: [RRR...][GGG...][BBB...]
// Better cache locality for per-channel operations
```

**Efficient Writes**:
```c
// Grayscale→RGB: Single memset for 3 bytes
memset(&out_image->data[i * 3], val, 3);
// ~8% faster than 3 separate writes
```

**Zero-Copy Operations**:
```c
// Direct CFITSIO to float buffer
fits_read_pix(fptr, TFLOAT, fpixel, npixels, NULL, float_data, NULL, &status);
```

### 5. I/O Optimizations

**Direct CFITSIO Access**:
- No binding layer overhead
- Single-pass FITS reading
- Automatic BZERO/BSCALE handling

**Efficient JPEG Encoding**:
- stb_image_write: Compact, optimized library
- Adjustable quality for speed/size tradeoff
- ~10-20ms for 2048×2048

## Profiling Results

**Typical Hotspots** (after all optimizations):
- `fits_read_img()`: 25-35% (I/O bound - CFITSIO)
- `super_pixel_debayer()`: 15-20% (memory bound)
- `find_median()`: 15-20% (CPU bound - quickselect)
- `apply_stretch_simd()`: 10-15% (SIMD optimized)
- `save_jpeg()`: 10-15% (JPEG encoding)

**Previous Hotspots** (now optimized):
- ~~`qsort()`: Was 30-40%~~ → Quickselect: 15-20%
- ~~`apply_stretch()`: Was 15-20%~~ → SIMD: 10-15%

## Libraries Used

**CFITSIO** (NASA HEASARC):
- Industry-standard FITS I/O library
- Handles all FITS formats (BITPIX 16, -32, etc.)
- Automatic data type conversions
- Compressed FITS support
- License: NASA Open Source Agreement

**stb_image_write** (Sean Barrett):
- Single-header JPEG encoder
- Public domain / MIT license
- Optimized for code size and simplicity
- ~10-20ms JPEG encoding

## Further Optimization Opportunities

### OpenMP Parallelization

Requires Homebrew LLVM on macOS:
```bash
brew install llvm
```

```c
#pragma omp parallel for
for (int c = 0; c < num_channels; c++) {
    compute_stretch_params(...);
    apply_stretch(...);
}
```
- **Potential**: 30-40% speedup for RGB images
- **Trade-off**: Adds dependency on OpenMP runtime

### libjpeg-turbo

Replace stb_image_write with libjpeg-turbo:
- Uses SIMD instructions (AVX2/NEON)
- 2-5× faster JPEG encoding
- **Potential**: 5-15ms saved per image
- **Trade-off**: Larger binary, more complex build

### Memory-Mapped I/O

```c
// mmap() FITS files instead of fits_read_img()
void* mapped = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
```
- **Potential**: 20-30% I/O reduction
- **Trade-off**: Increases code complexity, platform-specific

## Profiling Commands

**macOS Instruments**:
```bash
cargo instruments --release --bin rustafits -- input.fits output.jpg
```

**Linux perf**:
```bash
perf record ./target/release/rustafits input.fits output.jpg
perf report
```

**Benchmark**:
```bash
hyperfine --warmup 3 './target/release/rustafits input.fits output.jpg'
```

## Comparison to QuickLook.Plugin.FitsViewer

Our implementation achieves comparable performance:

| Metric | QuickLook (C++) | rustafits (C) |
|--------|----------------|---------------|
| 4K RGB | < 100ms | ~100ms |
| 6K Mono | ~200ms | ~220ms |
| 6K Mono (preview) | N/A | ~90ms |

**Why We're Competitive**:
1. Same algorithms (super-pixel debayering, KStars stretch)
2. Same library (CFITSIO)
3. Similar optimizations (O3, native CPU, fast math)
4. **Additional**: Quickselect, SIMD, preview mode

## Performance Tips

**For Fastest Previews**:
```bash
rustafits input.fits output.jpg --preview --quality 30
```

**For Quick Overviews**:
```bash
rustafits input.fits output.jpg --downscale 4
```

**For Production Quality**:
```bash
rustafits input.fits output.jpg --quality 100
```
