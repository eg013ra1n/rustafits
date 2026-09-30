# Changelog

All notable changes to rustafits will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- FITS headers and data layout are now parsed by
  [`fitsio-pure`](https://github.com/OrbitalCommons/fitsio-pure) (pure Rust);
  the byte-order and scaling conversions are unchanged, and output is
  bit-identical on every file the old reader accepted.

### Added

- Tile-compressed FITS (`.fz`: RICE_1, GZIP_1, GZIP_2), gzipped FITS
  (`.fits.gz`), the `.fts` extension, and images stored in the first
  extension HDU behind an empty primary.

## [1.3.0] — 2026-09-20

Fifteen commits since 1.2.0, from three cycles of the stacking work in the
host application: the M4a quality round (detection levels, the XISF reader),
the full-resolution debayer option, and two performance tiers on the analysis
kernels. Minor, not patch: the crate gains public API. One feature-gated struct
is marked `#[non_exhaustive]` — see **Changed** for why that is not a major.

### Added

- **`DetectionLevels` and `ImageAnalyzer::with_detection_levels`.** The fast
  detector used to pick its two ladder levels from a bright-pixel rank budget
  derived from `max_stars` — blind to sky brightness — and then descended
  through a fixed 30-σ arm and a per-tile adaptive pass until the cap was
  reached, so a caller that wanted a threshold-defined population had no way to
  ask for one. `DetectionLevels::NoiseRelative { k1, k2 }` puts the two levels
  at `background + k·noise` and runs only those two arms; `DetectionLevels::
  Absolute` takes the two levels verbatim in ADU above background, for a caller
  that pre-filters the image and must keep the threshold still while the filter
  reshapes the population. `RankBudget` stays the default and reproduces the
  previous behaviour exactly.
- **`ProcessConfig.vng_debayer` / `ImageConverter::with_vng_debayer`** route a
  CFA frame through the gradient debayer (`processing::vng`) instead of the
  super-pixel one, keeping the native pixel grid. Off by default; the two
  debayers reduce differently (super-pixel *is* a 2× downscale), and the new
  path applies the whole `downscale_factor` rather than inheriting the `/2`
  compensation.
- **`LmResult::iters`** (under `debug-pipeline`): the Levenberg–Marquardt
  solvers report their iteration count, so a caller can pin that an arithmetic
  change left the control flow untouched.

### Changed

- **`LmResult` is `#[non_exhaustive]`** (under `debug-pipeline`). It is a
  solver *result*: nothing outside the crate has a reason to construct one, and
  every in-tree destructure already uses `..`. The attribute lets the solvers
  grow report fields without a major bump each time. If you destructure it
  exhaustively, add `..`.
- **XISF float samples keep the u16-domain convention on purpose.** The reader
  briefly stopped scaling Float32/Float64 samples by 65 535 (the render pipeline
  is provably invariant to a uniform rescale), then restored it: other consumers
  of this reader spill the samples into ADU-domain scratch and hand them to
  ADU-domain detectors. The convention is now a documented cross-crate contract
  rather than an accident.

### Fixed

- **XISF: the largest `<Image>` wins, `byteOrder="big"` and `bounds="lo:hi"` are
  honoured, and an un-parsable `<Image>` candidate no longer fails the whole
  read.** A master written by an external tool carries a same-size weight map
  after the data; the reader used to take the first image it saw. The
  multi-image scan introduced a regression — a parse error in any candidate
  failed the read — which is closed with fixture tests; Float64 samples keep
  f64 arithmetic until the final cast.

### Performance

All bit-identical to the previous code unless stated; measured on real
26-megapixel frames on a 10-core machine.

- **The adaptive detector takes the caller's background/noise pair** instead of
  recomputing it, skips the dead noise-map interpolation, and borrows a
  single-channel luminance rather than copying it (detect −15…−18 ms per plane,
  background −40 ms).
- **`estimate_noise_mrs`** runs on ping-pong buffers (peak residency unchanged
  — the à-trous layers are the floor — but no per-layer allocation).
- **HFD measurement** (`hfd_at`) uses a single-copy scratch window with an
  interior fast path, a static r² table and selection-based medians; the
  histograms are built in parallel. Detection on a 6224×4168 frame −45 %.
- **The Moffat solver does one transcendental per sample-iteration**:
  `base^-β` is evaluated once and its derivative is `−β·power/base`; `powi`
  replaces `powf` when β is integral; the residual vector is reused from the
  last accepted Jacobian pass; the Cholesky solve runs on a thread-local
  scratch. Control flow untouched — over 180 solves (36 stamps × fixed
  β 2.5/4/6/10 + free β) the iteration count and convergence verdict are
  identical, positions move ≤ 7e-15 px and fluxes ≤ 2e-14 relative — while a
  fit takes 30–35 % less time. A `debug_assert` pins the one contract the reuse
  relies on (`params[7]` equals the fixed β).

## [1.2.0] — 2026-09-06

### Added

- **`FastStar::eccentricity` — the fast detector measures star shape.** It was
  computed only by the full analysis path, so anything reading `detect_fast`
  could not tell a star from a streak; a plate solver would happily build
  quads out of streak fragments and return a confident, wrong position. The
  moment maths is now shared (`shape_from_moments`) and runs for every fast
  detection: 0 is round, approaching 1 is a streak.

  The stamp follows the star's own size (2 × HFD, clamped) rather than a fixed
  width. This is the whole point of the change: a window narrower than the
  object sees only its bright core and reports it round, which is exactly what
  happens on wind-shaken frames whose stars are 13 px across.

  Read this rather than `sx`/`sy` when the question is "star or trail".
  Those come from the optional Moffat fit, which declines almost everything on
  a trailed frame — leaving them zero precisely when the shape matters most.
  Measured on real frames: healthy fields sit at median eccentricity 0.50,
  frames whose stars are streaks at 0.97.

  Minor, not patch, on the same reasoning as v1.1.0's `FastStar` additions:
  the struct gains a public field.

### Documentation

- `docs/platesolving.md` is new — what the crate does and does not provide for
  a solver (detection, quads, the shape gate above), written against the
  shipped API rather than an intended one. `docs/detection.md` and
  `docs/fitting.md` gained the sections that had drifted out of them.
- README, `CLAUDE.md` and `docs/usage.md` describe the API as it ships: the
  VNG entry point, `encode_jpeg`, and the fact that a one-shot-colour frame
  rendered through `ImageConverter` still goes through the half-resolution
  super-pixel path.

## [1.1.0] — 2026-09-05

### Added

- **`processing::vng::vng_debayer_f32` — full-resolution 8-gradient VNG
  demosaic.** Unlike `super_pixel_debayer_f32`, which averages each 2×2 CFA
  block and so halves both axes and broadens every star profile, VNG keeps the
  native pixel grid and returns planar RGB at the input geometry. This is what a
  one-shot-colour frame needs when the output is going into stacking rather than
  onto a screen. Validated against an external reference demosaic on real
  26 MP frames: median |difference| on interior pixels is 0 in green and ~3e-5
  in red and blue on a 0..65535 scale. Bit-identity with another implementation
  is not a goal — gradient thresholds are an implementation choice.

  The interior pass runs on rayon in 64-row bands. The split is by output row
  only and every band reads the same immutable mosaic, so the result does not
  depend on the thread count or on where a band boundary falls; a serial
  reference in the test module pins that. Measured on a 6248×4176 frame:
  3.9 s serial → 0.74 s.

- `encode_jpeg` — in-memory RGB/RGBA → baseline JPEG (4:2:0) encoding on the
  library surface, for embedders that need JPEG bytes without a file write
  (Athenaeum's Perseus preview). Backed by the same encoder `save_image` uses;
  RGBA alpha is read and discarded, so callers pass frames through without a
  de-interleaving copy.

- **Adaptive multi-level fast detection** replaces the single global threshold:
  a falling-threshold ladder with an occupancy mask, a saturation-aware level
  clip and a per-tile deep background, so long-focal-length and
  nebula-swamped fields still yield enough well-centroided stars to solve.
  `detect_fast`/`_data`/`_raw` use it; the analysis-UI detector is untouched.
- `build_quads_multi` — quad-pool densification for sparse fields via a
  configurable per-star group size. `build_quads` is preserved exactly as the
  `group_size = 4` case.
- Per-star SNR on `FastStar`. `detect_stars_adaptive` already computed an
  aperture-photometry SNR and discarded it; exposing it lets a solver re-rank
  its quad pool, because extended nebula and galaxy knots carry huge flux at
  low SNR and would otherwise crowd the top of a flux-ranked pool. No behaviour
  change on its own — the consumer opts in.
- Opt-in PSF centroid refinement: `ImageAnalyzer::with_centroid_refine`, plus
  `FastStar { sx, sy, fwhm }`. Default OFF, so `detect_fast` output is
  unchanged unless asked. `raw_x`/`raw_y` keep the pass-1 intensity-weighted
  centroid alongside the refined one, so a caller can compare the two per frame.
- `AnnotationConfig::ellipse_scale` (default 1.2× FWHM, was a hardcoded 2.5×).
  The old scale drew ~50 px lassos on oversampled frames, where clean single
  stars visually covered their neighbours and read as blends.
- `tracing` as a facade (no subscriber in library code): the analysis
  pipeline's already-measured stage durations are now also emitted as
  structured `debug!` events instead of being available only as returned
  `StageTiming`/`FastDetectTiming` data.
- `DataType` is re-exported, so consumers and tests can build an
  `ImageMetadata` without reaching into a private module.

### Changed

- **JPEG encoding moved off the `turbojpeg` C binding to a pure-Rust
  libjpeg-turbo reimplementation** (`libjpeg-turbo-rs`, pinned to `=0.8.0`).
  Output was **byte-identical** in every comparison run — the full-frame
  `cocoon`, `mono` and `osc` test images hash the same under both backends on
  aarch64, and a width sweep matches C libjpeg-turbo byte-for-byte on aarch64
  and x86_64 — and encode time is unchanged (69.1 ms vs 69.0 ms for a 6252×4176
  frame at q90 on Apple Silicon).

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

- **Out-of-range JPEG quality now clamps to 1..=100 instead of returning an
  error.** The `turbojpeg` binding rejected a quality outside 1..=100; the new
  encoder clamps, and `encode_jpeg` now makes that explicit. This also changes
  `ImageConverter::save_processed` (shipped in 1.0.1), which passes its
  `quality` argument straight through — a call that used to fail with
  `Image save failed` now succeeds at the clamped quality.
  `ImageConverter::with_quality` and the CLI already clamped/validated, so they
  are unaffected.

- Pass-2 measurement now rejects blended close pairs. Tight pairs used to
  collapse into one elongated pseudo-star, and slightly wider ones fragmented
  into two survivors each measuring the other's light: on a dense oversampled
  bench frame ~31 % of measured stars had a second significant peak inside
  their fitting window.

### Fixed

- **A frame containing NaN no longer renders black.** `compute_stretch_params`
  skips non-finite samples; quickselect is comparison-based, so a single NaN
  in the sample set poisoned the median and MADN and every stretch parameter
  came back NaN — with no error anywhere. Float FITS/XISF routinely carries NaN
  in registration borders and rejected pixels. A fully non-finite frame now
  falls back to neutral parameters; scattered NaN pixels still render black
  individually, which is the right display for missing data.
- **The fast detector's star cap is applied after a full scan, not in scan
  order.** `scan_region` returned as soon as `max_stars` detections had
  accumulated — row-major, so on a dense frame the fast path only ever saw a
  horizontal stripe (on `cocoon.fits`, all 500 stars sat above y = 1470 of
  4176, a 3 % overlap with the slow pipeline's top 100). Plate solving received
  spatially degenerate star sets, and registration fitted a transform on one
  stripe and extrapolated it across the frame.
- Per-star fit sigma is floored at the field sigma from bright-star
  calibration. On defocused frames the per-star half-max estimate collapsed and
  shrank the fit window to ~5 px, reporting 0.5–3 px FWHM for true 9–16 px
  stars. Defocused reference sets now land within ~5 % of an external
  measurement tool.
- `moments_ellipse` confines its initialising moments to the circular fitting
  disc instead of summing over the whole square stamp, where corners and
  neighbouring stars skewed the initial axis ratio and angle enough to strand
  the LM fit in a rounder local minimum. Matched-set eccentricity delta against
  an external per-star reference improves from a median of −0.019 to −0.007,
  and theta RMS from 19.7° to 13.8°.
- The FITS reader rejects `NAXIS3` outside {1, 3} with an error naming the
  offending value, instead of letting it flow through and corrupt downstream
  channel-count assumptions. `NAXIS < 3` keeps the mono default.
- XISF headers containing a dangling `&` in element text load again.
  quick-xml 0.41 rejects a lone ampersand where 0.36/0.37 read straight
  through; the reader opts back into permissive parsing. XISF headers are
  machine-written so this is rare, but a malformed one should still render
  rather than fail to load.
- XISF `Uint8` samples scale by 257, not 256, so 255 maps to 65535 — 8-bit data
  rendered about 0.4 % darker than the same image stored as `u16`.
- `ImageConverter::with_downscale` clamps its factor to ≥ 1; zero divided by
  zero inside the downscale kernels.

### Dependencies

- quick-xml 0.41, nalgebra 0.35, lz4_flex 0.14, ruzstd 0.9, base64 0.23.
- The unused `criterion` dev-dependency is gone.

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

[1.1.0]: https://github.com/eg013ra1n/rustafits/releases/tag/v1.1.0
[1.0.0]: https://github.com/eg013ra1n/rustafits/releases/tag/v1.0.0
