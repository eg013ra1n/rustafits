# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

rustafits is a high-performance FITS/XISF to JPEG/PNG converter for astronomical images, written in pure Rust. The crate name is `rustafits` but the library is exported as `astroimage` (`use astroimage::...`).

## Build & Test Commands

```bash
# Build
cargo build
cargo build --release

# Run tests
cargo test
cargo test <test_name>          # run a single test

# Build debug binary (exposes analysis internals)
cargo build --features debug-pipeline

# Run CLI
cargo run -- input.fits output.jpg [--downscale 2] [--quality 90] [--preview] [--annotate] [--no-debayer] [--log]

# Check without building
cargo check

# Clippy lint
cargo clippy
```

## Architecture

### Crate Structure

The package produces two binaries (`rustafits`, `rustafits-debug`) and a library (`astroimage`). The `debug-pipeline` feature flag makes `analysis` and `formats` modules public.

### Processing Pipeline

`converter.rs` (builder API) → `pipeline.rs` (orchestration) → processing modules → `output.rs`

Pipeline flow for u16 data:
1. Format reader (FITS big-endian / XISF little-endian) → `(ImageMetadata, PixelData)`
2. Debayer (if Bayer pattern detected) — 2x2 super-pixel, produces planar f32 RGB
3. Downscale (integer factor, operates on raw u16 before float conversion when possible)
4. u16→f32 conversion (SIMD-accelerated)
5. Preview binning (optional 2x2)
6. Auto-stretch (median-based STF: `((m-1)*x) / ((2m-1)*x - m)`, max_input always 65536.0)
7. Color conversion to interleaved u8 RGB/RGBA

### Internal Pixel Format

- **PixelData**: `Uint16(Vec<u16>) | Float32(Vec<f32>)` — owns allocations, no raw pointers
- **Planar f32**: channels stored contiguously as RRRGGGBBB (not interleaved)
- XISF float32 [0,1] is scaled by 65535 to match FITS u16 range

### Module Map

- `types.rs` — `PixelData`, `ImageMetadata`, `ProcessedImage`, `BayerPattern`
- `converter.rs` — `ImageConverter` builder with optional `Arc<rayon::ThreadPool>`
- `pipeline.rs` — orchestrates read → process → encode
- `output.rs` — JPEG/PNG file writing
- `formats/` — `fits.rs` (FITS reader), `xisf.rs` (XISF reader with zlib/LZ4/Zstd decompression)
- `processing/` — `stretch.rs`, `debayer.rs`, `binning.rs`, `downscale.rs`, `color.rs`
- `analysis/` — `background.rs` (mesh-grid + MRS wavelet), `detection.rs` (DAOFIND), `fitting.rs` (two-pass Moffat-primary PSF calibration, LmResult + fit_residual), `metrics.rs` (fit_residual per star), `snr.rs`, `convolution.rs`, `render.rs`, `mod.rs` (two-stage trail detection, residual-weighted statistics)
- `annotate.rs` — 3-tier annotation API (raw geometry / RGBA layer / burn-in)

### SIMD Strategy

SSE2 baseline on x86_64, AVX2 via runtime detection, NEON on aarch64. All operations have scalar fallbacks. SIMD-accelerated: stretch, binning, u16→f32, gray→RGB, debayer.

Key gotcha: `_mm_shuffle_epi8` (pshufb) requires SSSE3, not SSE2 — gray→RGB on SSE2 falls back to scalar.

### Parallelism (rayon)

Analysis pipeline stages parallelized via rayon: background mesh cells, separable convolution (row-parallel), peak detection (row-parallel), NMS sort, post-NMS stamp/filter processing, PSF measurement, per-star SNR. NMS grid scan stays sequential (greedy dependency). Blend rejection stays sequential (mutual marking).

### Thread Pool

`ImageConverter` accepts an optional `Arc<rayon::ThreadPool>`. Uses `pool.install()` to redirect all nested `par_*` calls. The library re-exports `rayon::ThreadPool` and `rayon::ThreadPoolBuilder`.

## Test Data

Test files in `tests/`:
- `cocoon.fits` — dense star field
- `mono.fits` — monochrome 16-bit (comet, bright extended object)
- `osc.fits` — OSC Bayer pattern
- `test.xisf` — XISF format

Always test changes against all four files when modifying processing code.

## Release

CI builds Linux x86_64 on tag push (`v*.*.*`). macOS is built manually (see `BUILD_MACOS.md`). Packages: Homebrew tap, AUR, Debian, RPM spec.
