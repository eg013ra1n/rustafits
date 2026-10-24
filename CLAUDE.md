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

# Run CLI (classic form, or batch: several inputs/dirs, -o DIR, -j N; see --help)
cargo run --release -- input.fits output.jpg [-d 2] [-q 90] [-p] [-a] [--no-debayer] [--vng] [-v]
cargo run --release -- dir/ -r -o out/ -n

# Check without building
cargo check

# Clippy lint
cargo clippy
```

## Architecture

### Crate Structure

The package produces two binaries (`rustafits`, `rustafits-debug`) and a library (`astroimage`). The `debug-pipeline` feature flag makes `analysis` and `formats` modules public. The default `cli` feature carries the `rustafits` binary's only extra deps (clap, indicatif); the binary requires it, and library users (Athenaeum) can build with `default-features = false`. `tests/cli.rs` drives the binary end to end on synthetic FITS.

### Processing Pipeline

`converter.rs` (builder API) → `pipeline.rs` (orchestration) → processing modules → `output.rs`

Pipeline flow for u16 data:
1. Format reader (FITS big-endian / XISF little-endian) → `(ImageMetadata, PixelData)`
2. Debayer (if Bayer pattern detected) — by default 2x2 super-pixel, producing
   planar f32 RGB at HALF the input width and height (so it consumes a factor 2
   of `--downscale`). `ImageConverter::with_vng_debayer()` switches the render
   pipeline to `vng_debayer_f32` instead: full native resolution, and the whole
   downscale factor still applies. `vng_debayer_f32` is also a library entry
   point callers invoke directly (Athenaeum's calibrated-lights export does).
   The CLI does not expose VNG.
3. Downscale (integer factor, operates on raw u16 before float conversion when possible)
4. u16→f32 conversion (SIMD-accelerated)
5. Preview binning (optional 2x2)
6. Auto-stretch (median-based STF: `((m-1)*x) / ((2m-1)*x - m)`, max_input always 65536.0)
7. Color conversion to interleaved u8 RGB/RGBA

### Internal Pixel Format

- **PixelData**: `Uint16(Vec<u16>) | Float32(Vec<f32>)` — owns allocations, no raw pointers
- **Planar f32**: channels stored contiguously as RRRGGGBBB (not interleaved)
- XISF float32 [0,1] is scaled by 65535 to match FITS u16 range, after
  `bounds="lo:hi"` normalization (identity for the default `0:1`); multi-
  `<Image>` headers resolve to the largest image, ties to the first. This
  is a PRODUCTION CONTRACT, not just a rendering convenience —
  `athenaeum-core`'s `integration::banded::spill_via_read_raw` (master
  builds, light calibration) and `analysis::analyzer::analyze_frame` both
  take `PixelData::Float32` straight into ADU-domain code with no rescale
  of their own, so removing this scale silently breaks them (a 65535×-too-
  small master; a detector that finds nothing) even though
  `processing::stretch`'s auto-stretch is separately provably invariant to
  a uniform rescale of its float input (`unit_range_floats_stretch_like_u16`).
  A caller that wants the file's own native units (M4a's
  `stacking::measure`, which assumes calibrated frames are float32 in
  `[0, 1]`) divides by 65535 itself — see `read_xisf_image`'s doc comment
  (M4a Task 1, controller ruling R-M4a-11) for the full chain of evidence,
  including the one-time attempt to remove this scale and why it was
  reverted.

### Module Map

- `types.rs` — `PixelData`, `ImageMetadata`, `ProcessedImage`, `BayerPattern`
- `converter.rs` — `ImageConverter` builder with optional `Arc<rayon::ThreadPool>`
- `pipeline.rs` — orchestrates read → process → encode
- `output.rs` — JPEG/PNG file writing, plus the public `encode_jpeg` in-memory
  entry point. JPEG goes through `libjpeg-turbo-rs` (pure Rust, NEON/AVX2, no
  cmake/nasm); PNG through `image`. RGBA is passed to the JPEG encoder
  unchanged — it drops alpha itself, so do NOT de-interleave to RGB first.
- `formats/` — `fits.rs` (FITS reader), `xisf.rs` (XISF reader with zlib/LZ4/Zstd decompression)
- `processing/` — `stretch.rs`, `debayer.rs` (super-pixel 2x2 + native-resolution green interpolation), `vng.rs` (full-resolution 8-gradient VNG demosaic, rayon over 64-row bands), `binning.rs`, `downscale.rs`, `color.rs`
- `platesolving/` — quad building and matching (`pattern_matcher.rs`, incl. `build_quads_multi`), `ransac.rs`, `wcs.rs`, `projection.rs` (gnomonic), `proper_motion.rs`, `transform.rs`, `types.rs`. Reference: `docs/platesolving.md`.
- `analysis/` — `background.rs` (mesh-grid + MRS wavelet), `detection.rs` (DAOFIND), `fitting.rs` (two-pass Moffat-primary PSF calibration, LmResult + fit_residual), `metrics.rs` (fit_residual per star), `snr.rs`, `convolution.rs`, `render.rs`, `mod.rs` (two-stage trail detection, residual-weighted statistics)
- `annotate.rs` — 3-tier annotation API (raw geometry / RGBA layer / burn-in)

### SIMD Strategy

SSE2 baseline on x86_64, AVX2 via runtime detection, NEON on aarch64. All operations have scalar fallbacks. The JPEG encoder dependency follows the same rule (its own NEON/AVX2 paths, scalar fallback) — its `simd` feature is ON by default, so never add `default-features = false` to it. SIMD-accelerated: stretch, binning, u16→f32, gray→RGB, debayer.

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

`processing::vng` carries its own guards: a serial reference implementation the
row-banded parallel version is diffed against, and a pin against an external
reference demosaic. Both live in the module's test section — a change to the
band split must keep the serial diff bit-identical.

`tests/jpeg_encoder.rs` needs no fixtures: it guards the JPEG encoder with
synthetic images and decodes them with `jpeg-decoder` — deliberately an
*independent* decoder, not the decode side of `libjpeg-turbo-rs`, so a symmetric
bug there cannot pass a self-round-trip. When bumping the pinned
`libjpeg-turbo-rs` version, also run it under an x86_64 target
(`cargo test --target x86_64-apple-darwin`): the encoder bug documented in
`docs/libjpeg-turbo-rs-issue.md` was byte-clean on aarch64 and only visible on
x86_64.

## Release

No system build dependencies: every dep is pure Rust, so packaging recipes and CI
jobs must not reintroduce cmake/nasm. CI builds Linux x86_64 on tag push
(`v*.*.*`). macOS is built manually (see `BUILD_MACOS.md`). Packages: Homebrew
tap, AUR, Debian, RPM spec.
