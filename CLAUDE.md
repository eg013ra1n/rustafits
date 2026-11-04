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
cargo test --features serde     # + star_measure's serde names and pass-1 encoding pins

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

The package produces two binaries (`rustafits`, `rustafits-debug`) and a library (`astroimage`). The `debug-pipeline` feature flag makes `analysis` and `formats` modules public. The default `cli` feature carries the `rustafits` binary's only extra deps (clap, clap_complete, indicatif); the binary requires it, and library users (Athenaeum) can build with `default-features = false`. `serde` (off by default) derives `Serialize` / `Deserialize` on the `star_measure` and `stats` types; `test-support` exposes `star_measure::StarFit::test_at` for downstream tests; `lm-counters` adds the debug-only LM solve counters, re-exported as `astroimage::lm_counters`. The binary lives in `src/bin/rustafits/`: `main.rs` (clap; conversion is the default command, plus `info`, `analyze`, `completions`), `batch.rs` (input collection, shared pool, worker loop, progress bar), `convert.rs`, `info.rs`, `analyze.rs`, `report.rs` (table/JSON/CSV). `analyze` runs `ImageAnalyzer` with its library defaults, and `tests/cli.rs` pins the CLI's numbers to `ImageAnalyzer::new().analyze()`. These are no longer the numbers Athenaeum shows: its Analysis tab moved to its own `star_measure` module (see Downstream Users). `tests/cli.rs` drives the binary end to end on synthetic FITS.

### Processing Pipeline

`converter.rs` (builder API) → `pipeline.rs` (orchestration) → processing modules → `output.rs`

Pipeline flow for u16 data:
1. Format reader (FITS big-endian / XISF little-endian) → `(ImageMetadata, PixelData)`.
   `bayer_pattern`, FITS and XISF alike, has odd `XBAYROFF` / `YBAYROFF`
   folded in (column / row swap of the 2×2 cell: `formats::shift_bayer_pattern`),
   so it addresses the stored pixels as read — in file row order, whatever
   `flip_vertical` says. FITS folds them onto `BAYERPAT`; XISF folds them onto
   whatever pattern its scan found (`BAYERPAT` or `<ColorFilterArray>`, both
   `<FITSKeyword>` offsets), and reads an offset as a number rounded (`0.6` is
   1) where FITS takes the integer prefix (`0.6` is 0): the XISF rule is the
   one Athenaeum's Analysis tab applied before phase 4, kept so its numbers do
   not move. `star_measure::raw` folds nothing itself.
2. Debayer (if Bayer pattern detected) — by default 2x2 super-pixel, producing
   planar f32 RGB at HALF the input width and height (so it consumes a factor 2
   of `--downscale`). `ImageConverter::with_vng_debayer()` switches the render
   pipeline to `vng_debayer_f32` instead: full native resolution, and the whole
   downscale factor still applies. `vng_debayer_f32` is also a library entry
   point callers invoke directly (Athenaeum's calibrated-lights export does).
   The CLI exposes it as `--vng`.
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
  builds, light calibration) and this crate's own
  `star_measure::raw::read_raw_frame` (the Analysis tab's reader, in this
  crate since phase 4) both take `PixelData::Float32` straight into
  ADU-domain code with no rescale of their own, so removing this scale
  silently breaks them (a 65535×-too-small master; a detector that finds
  nothing) even though
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
- `processing/` — `stretch.rs`, `debayer.rs` (super-pixel 2x2 + native-resolution green interpolation), `vng.rs` (full-resolution 8-gradient VNG demosaic, rayon over 64-row bands), `cfa.rs` (green-site mask / `G1 + G2` sum / single green phase, for measuring a raw mosaic without interpolation — Athenaeum's Analysis tab), `binning.rs`, `downscale.rs`, `color.rs`
- `platesolving/` — quad building and matching (`pattern_matcher.rs`, incl. `build_quads_multi`), `ransac.rs`, `wcs.rs`, `projection.rs` (gnomonic), `proper_motion.rs`, `transform.rs`, `types.rs`. Reference: `docs/platesolving.md`.
- `analysis/` — `background.rs` (mesh-grid + MRS wavelet), `detection.rs` (DAOFIND), `fitting.rs` (two-pass Moffat-primary PSF calibration, LmResult + fit_residual), `metrics.rs` (fit_residual per star), `snr.rs`, `convolution.rs`, `render.rs`, `mod.rs` (two-stage trail detection, residual-weighted statistics)
- `annotate.rs` — 3-tier annotation API (raw geometry / RGBA layer / burn-in)
- `star_measure/` — Athenaeum's one star-measurement path (per-star Moffat fit, per-frame β, frame aggregates, sky, noise, the Analysis budget, trail / defocus verdicts), copied from `athenaeum-core` at 0b6fb183 (port plan phase 2) with paths, serde gates and visibility changed, nothing else. `adapter` lists the `#[doc(hidden)]` helpers Athenaeum's stacking calls (outside semver); `testing` holds `StarFit::test_at`. Port plan phase 4 (copied at e6344ddd, the plan names what changed beyond paths) added the Analysis-mode entry, the semver-covered part of the module: `raw` (`read_raw_frame`: a FITS / XISF file as a `RawFrame` in the file's ADU, header facts from `ImageConverter::read_header`), `analysis_mode` (`FrameOptions`, `measure_raw_frame`, `FrameSummary`: the measurement and the frame figures the Analysis tab shows), and under them `green` (a mosaic measured on its green sites) and `trail` (trail R²). Its docs cite Athenaeum's specs, rulings and modules.
- `stats.rs` — robust statistics (median / MAD / BWMV, the 1/16 stratified sample, noise-scale factors) from Athenaeum's `integration::stats`, which re-exports them.
- `geometry/` (crate-private) — the kd-tree `star_measure` marks blends with, and the test RNG; `test_support.rs` (tests only) — Athenaeum's fixture fields and log captures, byte-identical.

### Downstream Users

**Athenaeum** (`athenaeum-core`, through its `rustafits` submodule): `default-features = false, features = ["serde"]`, plus `lm-counters` when its own optional `lm-counters` feature is on (for its probe), and `test-support` as a dev-dependency. It uses:
- the readers (`ImageConverter::read_raw` throughout);
- `star_measure` (`pub use astroimage::star_measure`): its stacking and LN measure through the module, and, since port plan phase 4, its Analysis tab, Blink and plate-solve gate measure through `analysis_mode` (`raw::read_raw_frame`, `measure_raw_frame`, `FrameSummary`), keeping only the row mapping, the ids and `config_hash`;
- the render pipeline and `encode_jpeg`;
- `processing::vng` (calibrated-lights export) and `processing::cfa` (the Analysis tab's green-site path);
- the `platesolving` types (`WcsSolution`, `SipCoefficients`, `SolveHints`, `ProperMotionCorrector`);
- `ImageAnalyzer::detect_fast_data`, for its registration detector and its `Peak` seed source.

Its own `star_measure` module measures stars for stacking, LN and, since its Plan C (2026-10-24), the Analysis tab, Blink and the plate-solve gate. It is built on `analysis::fitting` (`fit_moffat_2d_fixed_beta_conditioned`, `fit_moffat_2d_bounded`, `Moffat2DResult`, `PixelSample`) and `analysis::background` (`estimate_background_mesh`, `estimate_noise_mrs`). Those functions are therefore a de facto public API even though `analysis` is private without `debug-pipeline`: a change that moves one bit of their output moves Athenaeum's stacking weights and Analysis numbers (its `MEASURE_VERSION` cache key). That is what the bit pins in `analysis::fitting`'s tests guard.

This crate now holds a copy of that module (`star_measure`, `stats`). Until Athenaeum switches to it (port plan phase 3), Athenaeum's own copy is frozen and the two must stay identical: a change to one goes into both in one pair of commits. After the switch this copy is the only one, and `MEASURE_VERSION`'s policy applies here: one bump per merged plan that changes what a fit accepts or measures, never one per commit; a change that moves no number keeps it. Athenaeum's truth lab (athenaeum-lab `star_probe` scorecard) stays the referee for any change to these numbers.

The full `ImageAnalyzer::analyze*` is no longer called by the app, only by Athenaeum's lab (`star_probe`'s `legacy-analysis` arm) and a probe-only stacking arm. The plan for moving `star_measure` here and switching Athenaeum and the CLI to it is `docs/superpowers/plans/2026-11-03-star-measure-port.md`.

**solvemyastro** (the plate solver Athenaeum uses; path dependency, default features) uses rustafits only to decode frames and to detect stars: `ImageAnalyzer::detect_fast` and `FastStar`.

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

The bit-for-bit pins in `analysis::fitting`'s tests are recorded on aarch64
and compared exactly there; elsewhere `assert_bits_pinned` allows
1e-12 · max(1, |pin|) per field (x86_64 libm rounds differently; round-star θ
and near-zero fields drift most). Re-record pins on aarch64.

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
