# Changelog

All notable changes to rustafits will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `encode_jpeg` — in-memory RGB/RGBA → baseline JPEG (4:2:0) encoding on the
  library surface, for embedders that need JPEG bytes without a file write
  (Athenaeum's Perseus preview). Backed by the same encoder `save_image` uses;
  RGBA alpha is read and discarded, so callers pass frames through without a
  de-interleaving copy.

### Changed

- **JPEG encoding moved off the `turbojpeg` C binding to a pure-Rust
  libjpeg-turbo reimplementation** (`libjpeg-turbo-rs`, pinned to `=0.8.0`).
  Output is **byte-identical** — verified on the full-frame `cocoon`, `mono` and
  `osc` test images, whose encoded JPEGs hash the same under both backends — and
  encode time is unchanged (69.1 ms vs 69.0 ms for a 6252×4176 frame at q90 on
  Apple Silicon).

  The point is the build, not the speed: `turbojpeg-sys` compiled libjpeg-turbo
  from source and therefore required **cmake and nasm**. Those are now gone from
  the README install instructions, the Homebrew formula, the PKGBUILD, the RPM
  spec and the CI jobs — `cargo install rustafits`, cross-compilation and distro
  packaging need nothing beyond the Rust toolchain. The crate is now what
  `rustafits.spec` already claimed: no C dependencies.

  0.8.0 is the first release we build on because it carries the fix for a 4:2:0
  trailing-MCU divergence from C libjpeg-turbo that we reported upstream
  (issue #362, fixed in 0.7.0). See `docs/libjpeg-turbo-rs-issue.md` for the
  verification table and the re-check procedure for future version bumps.

- `save_image` no longer builds an intermediate RGB copy for RGBA frames: the
  encoder reads 4 bytes per pixel and discards alpha itself, so a `W*H*3`
  allocation per frame is gone.

## [1.0.1] — 2026-05-07

### Fixed

- **FITS `ROWORDER` interpretation now matches astronomical-software
  convention.** The `flip_vertical` flag computed by the FITS reader had
  its truth condition inverted: it was set when `ROWORDER='TOP-DOWN'`
  (the case that *doesn't* need a flip) and cleared when `ROWORDER` was
  absent or `'BOTTOM-UP'` (the cases that *do*). The corrected rule:

  | `ROWORDER`     | `flip_vertical` |
  | -------------- | --------------- |
  | `TOP-DOWN`     | `false`         |
  | `BOTTOM-UP`    | `true`          |
  | (missing)      | `true`          |
  | (other)        | `true`          |

  This matches PixInsight, ds9, Siril, and AstroImageJ — all of which
  treat the absence of `ROWORDER` as bottom-up by default and flip for
  screen display. The comparison is now case-insensitive so non-canonical
  spellings (`'top-down'`) are still recognized.

### Behavior change for downstream consumers

This is a **semantics change** for the existing `ImageMetadata.flip_vertical`
field, not an API change. Files that previously rendered "right-side up"
under the inverted rule will now render flipped 180° vertically (and vice
versa). For libraries written by N.I.N.A. and similar capture software
(which set `ROWORDER='TOP-DOWN'` explicitly), the displayed JPEG will
*no longer* be flipped — the on-screen orientation will match what
PixInsight shows for the same file.

Consumers that pass `flip_vertical` through to a UI overlay (e.g., star
annotations) automatically stay aligned with the rendered image because
the flag still flows from a single source.

## [1.0.0] — 2026-04-14

First stable release. Establishes a public API commitment: additive changes may
land in 1.x; breaking changes will require 2.0.

### Added

- **Plate solving module** (`astroimage::platesolving`) — quad-based pattern
  matching, gnomonic (TAN) projection, proper-motion propagation, RANSAC
  outlier filter, WCS solution (`WcsSolution`), and similarity/affine/SIP
  transform fitting (`TransformFitter`). Reusable for both catalog-matched
  plate solving and frame-to-frame star registration. rustafits does not
  touch the filesystem — callers pass pre-loaded star lists.
- **Fast star detection** — `ImageAnalyzer::detect_fast`, `detect_fast_data`,
  and `detect_fast_raw` methods that produce lean `FastStar { x, y, peak, flux }`
  centroids in ~300–500 ms on a full-frame image (release build). The pipeline
  skips PSF calibration, second-pass detection, Levenberg–Marquardt Moffat
  fitting, SNR photometry, and trail detection. Intended for pipelines that
  only need positions and brightness ordering — blind plate solving, quad
  hash matching, quick previews.
- New public types: `FastStar`, `FastAnalysisResult`, `FastDetectTiming`.
- New dependency: `nalgebra = "0.33"` (used by the plate-solving transform fit).
- `criterion` as a dev-dependency for future benchmarks.

### Changed

- README gains a "Fast star detection" section with a code example, method
  table, pipeline description, and field tables for the three new result types.

### Removed

- Orphan benchmark `benches/plate_solve.rs` (referenced an API from an earlier
  design iteration that was superseded during development).

### Notes

- No migration required. All changes are additive. Existing consumers of
  `ImageAnalyzer::analyze` and the precise pipeline are unaffected.
- Fast detection produces pass-1 centroids (~0.3–0.5 px accuracy). Use
  `analyze` when you need FWHM, eccentricity, HFR, SNR, or any PSF metric.

[1.0.0]: https://github.com/eg013ra1n/rustafits/releases/tag/v1.0.0
