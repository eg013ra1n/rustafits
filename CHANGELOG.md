# Changelog

All notable changes to rustafits will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
