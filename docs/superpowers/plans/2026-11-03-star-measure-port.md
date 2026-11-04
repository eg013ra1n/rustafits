# Port Athenaeum's `star_measure` into rustafits — Plan

> **Status:** approved, implementation not started. The owner took decisions D-1…D-9 on 2026-11-03, every one as
> recommended (section 4). Facts below were read on 2026-11-03 at Athenaeum `main` 0b6fb183, rustafits `main`
> 66b97b1 (Athenaeum's submodule pin) and solvemyastro `main`.

## Background

On 2026-10-24 Athenaeum's Plan C switched the Analysis tab, Blink and the plate-solve input gate from rustafits'
`ImageAnalyzer` to `athenaeum-core`'s own `star_measure` module, the module its stacking Measure already used.
The owner picked **option A**: rustafits is not touched, its analyzer stays here unused by the app, and "the owner
switches rustafits himself later". This plan is that switch.

Sources in the Athenaeum repo:
- `docs/superpowers/plans/2026-10-24-one-star-path-plan-c.md`, l.22–23 and l.45;
- `docs/stacking/README.md`, R-OSP-17;
- `docs/superpowers/specs/2026-10-07-one-star-measurement-path-design.md`, D-1 (l.66–68) and §2.2 (l.270).

The texts of options A/B/C themselves are in Athenaeum's git-ignored
`.superpowers/sdd/2026-10-24-plan-c/rulings.md`. This plan does not depend on them.

## 1. Goal

There is one star-measurement implementation, and it lives in rustafits:

- Athenaeum's stacking Measure, the LN refit, the Analysis tab, Blink and the plate-solve gate call
  `astroimage::star_measure` instead of `athenaeum_core::star_measure`, and every number stays bit-identical.
- `rustafits analyze` and `--annotate` report what Athenaeum's Analysis tab shows for the same file. 1.6.0 built
  the CLI around that property and lost it the same day, when Athenaeum switched.

**Non-goals:**
- Changing any number.
- Changing the fast detector (`ImageAnalyzer::detect_fast*` / `FastStar`). solvemyastro uses it, and so do
  Athenaeum's registration detector and its `Peak` seed source.
- Changing the plate solver.
- Plan D's open items.
- Changing the truth lab. athenaeum-lab stays the referee.

## 2. Where things stand

### 2.1 The module

`crates/athenaeum-core/src/star_measure/` has 12 files: 12 132 lines before the test modules, plus 291 `#[test]`s.

| file | lines (non-test) | role |
| ---- | ---- | ---- |
| `mod.rs` | 2180 | `measure_plane`, `prepare_plane` / `finish_plane` / `plane_beta`, `cross_plane_flags`, sibling seeds, `fit_stars_at`, `fit_stars_with_model`, `correct_at_ring`, `MeasureOptions`, `Mode`, `SeedSpec`, `PlaneMeasurement`, `PixelNoise`, `SnrNoise`, `MEASURE_VERSION` (= 8) |
| `seeds.rs` | 2763 | the matched filter, field width, knot flags, noise correlation |
| `fit.rs` | 2709 | per-star Moffat fit, acceptance, canonical axes, blend marking |
| `photometry.rs` | 1667 | HFR, SNR, ring sky, aperture correction `c`, saturation verdict |
| `robust.rs` | 753 | robust location / deviation helpers |
| `domain.rs` | 404 | working units, sky, noise model (mesh + MRS) |
| `population.rs` | 381 | population / list selection (probe and lab) |
| `substeps.rs` | 342 | sub-step timers |
| `budget.rs` | 287 | the Analysis budget (≤ N stars over a 4 × 4 grid) |
| `beta.rs` | 246 | per-frame β |
| `frame.rs` | 203 | frame aggregates |
| `diagnostics.rs` | 197 | refusal reasons and counts |

### 2.2 What it depends on

**Production code:**

| dependency | used in | after the port |
| ---- | ---- | ---- |
| rustafits `analysis::fitting`: `fit_moffat_2d_fixed_beta_conditioned`, `fit_moffat_2d_bounded`, `Moffat2DResult`, `PixelSample` | fit, photometry, population | `crate::analysis::fitting`. It is already here; Athenaeum reaches it today through `debug-pipeline`. |
| rustafits `analysis::background`: `estimate_background_mesh`, `estimate_noise_mrs` | domain | `crate::analysis::background` |
| Athenaeum `integration::stats`: `MAD_TO_SIGMA`, `CLIP_LO/HI`, `NOISE_CLIP_LO/HI`, `median_in_place`, `median_of`, `mad_about`, `avg_dev_about`, `clip_sample`, `location_scale`, `LocationScale`, `ScaleEstimator`, `for_each_stratified`, `stratified_sample`, `noise_scale_factors_within` | domain, fit, mod, photometry, seeds | a new `astroimage::stats`. It takes the first ~230 lines of the 625-line file; the normalization pairs after them stay in Athenaeum. |
| Athenaeum `geometry::kdtree::KdTree2` (227 lines) | fit (`mark_blend_neighbours`) | a private copy here. Athenaeum keeps its own for registration. |
| `serde` derives on 24 types; `ts_rs::TS` on `PsfModel` (and on `ScaleEstimator`, in `integration::stats`) | 8 files | D-4, D-5 |
| `libc` (`clock_gettime` of the thread CPU clock, unix) | substeps | a regular rustafits dependency (it was the `cli` feature's) |
| `rayon`, `tracing` | everywhere | already rustafits dependencies |

**Thirteen `pub(crate)` items** are called from elsewhere in Athenaeum, so they must stay reachable from it:

| item | called from |
| ---- | ---- |
| `install_in_span` | `stacking::{measure, ln}` |
| `PlanePrep::seeds`, `verdict_time`, `verdict_iterations` | `stacking::measure` |
| `PlanePrep::is_storable`, `PlanePrep::for_pin` | `stacking::pass1_artifact` |
| `frame::finite_median` | `analysis::analyzer` |
| `fit::inside_half_max` | `ln::scale` |
| `fit::mark_blend_neighbours` | `stacking::common_set` |
| `seeds::envelope_floor_factor` | `stacking::integrate` |
| `seeds::read_noise_correlation`, `NoiseCorrelation`, `NoiseCorrelation::rho` | `stacking::measure` (tests) |

(Corrected in phase 2: the first count listed `aperture_correction`, `detect_with` and `ring_sky`, whose matches
outside the module are a struct field, a doc link and a step name, and missed `PlanePrep::seeds`.)

**Test code:**
- **Field generators** in `test_support`: `add_noise` (21 uses), `fixture_field` (7), `width_spread_field` (+ `_beta`),
  `trailed_field`, `nan_border_field`, `SpreadStar`, `trail_plane_fixture`, `spread_moffat_field`, `gaussian_field`,
  `moffat_field` / `MoffatStar`, plus `geometry::ransac::SplitMix64`. These get a `#[cfg(test)]` copy here.
- **Log capture**: `stacking::solve::tests::captured_logs` and `stacking::telemetry::capture_logs` (13 call sites).
  These get a local helper, and `tracing-subscriber` becomes a rustafits dev-dependency.
- **Three tests** need Athenaeum's own code, so they stay there (phase 3 moves them out of `star_measure/`):
  - two adapter-identity tests, into `stacking::measure`'s tests: `module_measure_equals_the_adapter_bit_for_bit`
    and `adapter_matched_filter_seeds_fit_exactly_as_the_module_matched_filter` (both feature `solver`);
  - `robust::integration_copies_agree_with_this_module`, which compares `robust` with `integration::combine` and
    `integration::student_t`, into the integration engine's tests.
- **Serialization pins**: `serde_json` names and the `postcard` encoding of `PlanePrep` (byte length + xxh3), the
  first-pass artifact Athenaeum stores. They run under `--features serde`, with `postcard` and `xxhash-rust` as
  dev-dependencies.
- **`impl StarFit { fn test_at(..) }`** in Athenaeum's `test_support.rs` is an inherent impl on what becomes a
  foreign type, so it no longer compiles there. It has 6 call sites, in `fit.rs` and `stacking/fits_artifact.rs`.
  See D-3.

### 2.3 Athenaeum's Analysis-mode layer (Plan C)

`crates/athenaeum-core/src/analysis/`:

| file | lines (non-test) | what | port? |
| ---- | ---- | ---- | ---- |
| `input.rs` | 285 | `read_raw_frame` → `RawFrame { planes: Mono / Cfa / Rgb, w, h, domain, flip_vertical, channels }`: full scale, clip and gain (PC-2, R-PC-2); XISF `XBAYROFF` / `YBAYROFF` fold | yes. Today its header facts (BITPIX, BZERO, BSCALE, DATAMAX, EGAIN, X/YBAYROFF) come from Athenaeum's `fits_parser`; after the port they come from `ImageConverter::read_header`. The XISF fold moves into the XISF reader, as #7 did for FITS. |
| `cfa.rs` | 125 | `measure_green_sites` on `processing::cfa` (PC-4, R-PC-3, R-PC-8) | yes |
| `trail.rs` | 26 | `trail_r_squared`: Rayleigh R² over the canonical θ | yes |
| `analyzer.rs` | 219 | `measure_raw_frame(_with)` → `PlaneMeasurement`; `map_measurement` → `FrameAnalysis` / `StarMetric` rows | the frame arithmetic yes, as a library `FrameSummary`: point-like `stars_detected`, medians with PC-7's zero rule, background, noise, frame SNR, SNR weight, PSF signal, trail R², verdict, refused count, β under PC-9. The row mapping (ids, `config_hash`, `analyzed_at`) stays in Athenaeum. |
| `config.rs` | 142 | `AnalysisConfig` with sqlite load / save and a sha2 hash | no. Rustafits takes a plain options struct. |

### 2.4 Who calls rustafits' analysis today

| caller | API | role |
| ---- | ---- | ---- |
| Athenaeum `star_measure` | `analysis::fitting`, `analysis::background` (through `debug-pipeline`) | the Moffat fits, background mesh and MRS noise behind every stacking and Analysis number |
| Athenaeum `stacking::measure` | `ImageAnalyzer::detect_fast_data`; `analyze_data` | the stored `Peak` seed source; `SeedOverride::Full` (probe-only arm) |
| Athenaeum `stacking::register::detect` | `ImageAnalyzer::detect_fast_data` | the registration detector |
| Athenaeum examples: `star_probe` (`legacy-analysis` arm), `detection_audit`, `integrate_probe` | `ImageAnalyzer`; `analysis::fitting::lm_counters` | the lab |
| solvemyastro | `ImageAnalyzer::detect_fast`, `FastStar` | plate-solve detection (decode + detect only) |
| rustafits CLI `analyze`, `--annotate`, `rustafits-debug` | `ImageAnalyzer::analyze*` | — |

Athenaeum's app no longer calls `ImageAnalyzer::analyze*`; its lab and its probe-only `Full` arm do.

## 3. Target layout

```
astroimage::stats              robust statistics (pub; Athenaeum's integration::stats re-exports it)
astroimage::star_measure       the module: Athenaeum's submodule tree and names, unchanged
    ::adapter                  #[doc(hidden)] re-exports of the 13 helpers Athenaeum's stacking calls (D-3)
    ::raw                      RawFrame, RawPlanes, read_raw_frame          (from analysis/input.rs)
    ::green                    measure_green_sites, native_seed             (from analysis/cfa.rs)
    ::trail                    trail_r_squared                              (from analysis/trail.rs)
    ::analysis_mode            FrameOptions, measure_raw_frame(_with), FrameSummary  (from analysis/analyzer.rs)
astroimage::analysis           unchanged: fitting, background, detection / detect_fast*, legacy ImageAnalyzer
```

The submodule names that Athenaeum does not have (`adapter`, `raw`, `green`, `trail`, `analysis_mode`) are
proposals. After the port, Athenaeum has:
- `pub use astroimage::star_measure;` in `lib.rs`, so every `crate::star_measure::…` path compiles unchanged;
- `integration::stats` re-exporting the moved primitives;
- in `analysis/`: `config.rs` and the row mapping from `FrameSummary`.

## 4. Owner decisions

Taken by the owner on 2026-11-03, every one as recommended.

- **D-1 Scope: (b)**, the module plus the Analysis-mode entry and the CLI switch (phases 0–5). The alternative was
  (a), the module only (phases 0–3); it would not give the CLI Athenaeum's numbers.
- **D-2 Home after the port.**
  - rustafits is the only home: changes to `star_measure` land as rustafits PRs, and Athenaeum bumps the submodule;
  - `MEASURE_VERSION` moves with the code, and its policy moves into rustafits' CLAUDE.md ("one bump per merged
    plan that changes what a fit accepts or measures"); Athenaeum's cache keys keep reading the constant;
  - the truth lab stays in athenaeum-lab and stays the referee.
- **D-3 API surface.** The module's surface is shaped by one consumer: every `MeasureOptions` field is public, plus
  the 13 helpers. Therefore:
  - `pub mod star_measure`, always compiled (the CLI needs it from phase 5);
  - the helpers re-exported from a `#[doc(hidden)] pub mod adapter` and documented as outside semver;
  - `StarFit::test_at` behind a `test-support` feature that Athenaeum enables for its tests only;
  - the stable, semver-covered entry is `analysis_mode`.
- **D-4 serde.** 24 types derive `Serialize` / `Deserialize`; some are on Athenaeum's wire (`PsfModel` camelCase in
  the stacking config, the pass-1 / `.athf` artifacts). They get a `serde` feature, off by default and enabled by
  Athenaeum, with every derive and every `#[serde(...)]` attribute under `cfg_attr`. Athenaeum's round-trip tests
  pin the wire.
- **D-5 `ts_rs::TS` on `PsfModel`** (and on `ScaleEstimator`, found in phase 1). Athenaeum cannot derive a foreign
  trait on what becomes a foreign type, so it keeps wire enums of its own, with `From` both ways; rustafits gets no
  `ts-rs` feature. Regenerating the TS contract must leave `src/types` unchanged.
- **D-6 The legacy `ImageAnalyzer::analyze*`.**
  - kept through 1.x, marked `#[deprecated]` after phase 5;
  - removed in 2.0 (phase 6), once the lab retires `legacy-analysis` and Athenaeum retires `SeedOverride::Full`;
  - `detect_fast*` and `FastStar` stay in every case.
- **D-7 CLI `analyze` output.** Today: STARS FWHM HFR ECC SNR WEIGHT BG NOISE TRAIL. Analysis mode adds
  `stars_refused`, `frame_verdict`, `measure_version`, `psf_signal`, `frame_snr` (a ratio, PC-11) and `median_beta`
  (NULL under Moffat 4, PC-9). Therefore:
  - VERDICT replaces TRAIL in the table (PC-15's rule);
  - JSON / CSV carry every field;
  - `--max-stars` becomes the budget (default 500);
  - new `--detection-threshold` (σ, default 6.0, range 2–20) and `--psf-model moffat4|auto`.
- **D-8 Version.**
  - 1.7.0: phases 1–2, additive;
  - 1.8.0: phase 5, with a **Changed** entry for the CLI's numbers. Phase 4's rustafits additions ride with it,
    unless Athenaeum needs them first; then they get a minor release of their own;
  - 2.0: phase 6, the removal under D-6.
- **D-9 History.** The files are copied at a fixed Athenaeum sha, named in the commit message, rather than imported
  with history. Athenaeum's history stays the record, and its plans and rulings cite Athenaeum shas anyway.

## 5. Gates

- **G-1 Bit identity.** After Athenaeum switches (phase 3), and again after phase 4:
  - `cargo test -p athenaeum-core --all-targets` passes with no pin edited, run as at Plan B B6;
  - `MEASURE_VERSION` is unchanged;
  - `star_probe`'s scorecard shows 0 rows moved on every arm (stacking all-stars, analysis-mode budget) over the
    star-test corpus;
  - Athenaeum's stage hashes are unchanged, so a stored frame set is not re-measured.
- **G-2 Freeze.** From the phase-0 sha until phase 3 merges, Athenaeum does not change:
  - `star_measure/`;
  - the ported part of `integration/stats.rs`;
  - `geometry/kdtree.rs`;
  - `analysis/{input,cfa,trail,analyzer}.rs` (until phase 4 merges).

  A change that cannot wait goes into both copies in one pair of commits.
- **G-3 rustafits.** The following all pass with no new warnings, and render hashes are unchanged:
  - `cargo test`;
  - `cargo test --features debug-pipeline`;
  - `cargo build --no-default-features`;
  - `cargo clippy`.

  No clippy-driven rewrite of moved float code: a lint the moved code trips is allowed at the module top, never
  "fixed" in a way that could move a bit.
- **G-4 Kernels.** `analysis::fitting` and `analysis::background` do not change during the port: their numbers are
  Athenaeum's numbers.
- **G-5 Pins.** Exact on aarch64, `assert_bits_pinned` tolerance elsewhere (CLAUDE.md).
- **G-6 Thread count.** No result reads a thread count (Athenaeum Plan A's rule). The rayon reductions move
  unchanged.

## 6. Phases

### Phase 0: freeze and baseline (Athenaeum, lab)

- [ ] 0.1 Pick the source sha: Athenaeum `main` once no open branch touches the files of G-2. Record it here.
- [ ] 0.2 Record the baseline at that sha:
  - the full Athenaeum test run;
  - the lab scorecard `before` (stacking and analysis-mode arms, star-test corpus);
  - the stage hashes of one stored frame set.
- [ ] 0.3 Start the freeze (G-2).

### Phase 1: `astroimage::stats` (rustafits, additive)

Source sha: Athenaeum `main` 0b6fb183 (owner, 2026-11-04; the freeze of G-2 runs from it).

- [x] 1.1 Copy the primitives listed in 2.2 verbatim from `integration/stats.rs`, with their unit tests (they need
  `add_noise` / `SplitMix64` test copies). Done: lines 1–236 (everything above the normalization pairs), the 8
  primitive tests, and a `ScaleEstimator` serde-name test cut from `serde_names_match_the_spec`.
- [x] 1.2 Make it a public module whose doc carries the same contract text as Athenaeum's.
- Verify: G-3.

### Phase 2: `astroimage::star_measure` (rustafits, additive)

- [x] 2.1 Copy the 12 files verbatim at the sha. The only edits are paths:
  - `crate::integration::stats` → `crate::stats`;
  - `astroimage::analysis::` → `crate::analysis::`;
  - `crate::geometry::kdtree` → a private `kdtree` copy.

  Review with `git diff --no-index --word-diff` against the originals: nothing but those paths may differ.
  Done by `scripts/star_measure_port/port_star_measure.py`, so the copy can be regenerated: the paths above (the kd-tree and the test RNG keep Athenaeum's
  `crate::geometry::…` paths in a crate-private `geometry` module), the two log captures →
  `crate::test_support::…`, the module doc's provenance paragraph, and the clippy allow-list (G-3).
  `analysis::{fitting, background}` become `pub(crate)` without `debug-pipeline`, so the module reaches them.
- [x] 2.2 serde under `cfg_attr` (D-4); `ts_rs::TS` dropped (D-5). `fit::not_measured`, which only serde's
  `default = …` calls, is gated with it.
- [x] 2.3 The 13 helpers become `pub` and `#[doc(hidden)]` at their own paths (Athenaeum's paths keep working) and
  are re-exported from `star_measure::adapter`; `StarFit::test_at` goes behind `test-support` in
  `star_measure::testing` (D-3).
- [x] 2.4 Tests: every test except the three that stay in Athenaeum (2.2), with paths edited and nothing else.
  They use the local log captures, a `tracing-subscriber` dev-dependency, and Athenaeum's generators copied into a
  `#[cfg(test)]` `test_support.rs`. The serialization tests carry `#[cfg(feature = "serde")]`; one test carries
  `#[allow(clippy::almost_swapped)]` (its second assignment is a no-op).
- [x] 2.5 Make `analysis::fitting::lm_counters` public under `lm-counters` without `debug-pipeline`
  (`astroimage::lm_counters`), so phase 3 can drop `debug-pipeline` from Athenaeum.
- [x] 2.6 CLAUDE.md:
  - the module map line;
  - the `MEASURE_VERSION` policy (D-2);
  - "this module's numbers are Athenaeum's stacking numbers".
- Verify: G-3, and note the release binary's size change.
  - `scripts/star_measure_port/port_digest.rs` runs every entry point (`measure_plane` in `All`, `Budget` and
    `Population`, `prepare_plane` / `plane_beta` / `finish_plane`, `fit_stars_at` / `correct_at_ring`) over three
    synthetic fields and crops of the three test frames, and writes the results in full. Phase 3 reruns it.
  - Result at 0b6fb183 (x86_64, toolchain 1.96.1): the 288 copied tests match Athenaeum's own run name for name and
    outcome for outcome; the digest (31 MB, ~15k fits on the real frames) is byte-identical between Athenaeum's copy
    (dev profile, opt-level 2) and rustafits' (release, LTO), sha256 `4f907cab…`. Renders, `info` and `analyze` are
    byte-identical to `main`; the release binary grows by 752 bytes (the CLI does not call the module yet).
- Release 1.7.0, or have Athenaeum pin the git sha.

### Phase 3: Athenaeum switches the module (Athenaeum branch)

- [ ] 3.1 Bump the submodule with `default-features = false, features = ["serde"]`, plus `test-support` for its
  tests.
- [ ] 3.2 Code moves:
  - `pub use astroimage::star_measure;` in `lib.rs`;
  - delete `src/star_measure/`;
  - move the two adapter-identity tests into `stacking::measure`.
- [ ] 3.3 `integration::stats` re-exports the moved primitives. The normalization code stays.
- [ ] 3.4 The `PsfModel` wire enum (D-5); TS contract unchanged.
- [ ] 3.5 Logging. `star_measure` events now carry the target `astroimage::star_measure::…` instead of
  `athenaeum_core::star_measure::…`. Check:
  - the logging field dictionary (`docs/logging/README.md`), log-mcp queries, and every test that asserts a target;
  - that the sub-step timers still switch on. `substeps.rs` keys them on
    `tracing::enabled!(target: "athenaeum_core::stacking::measure", …)`, the adapter's target, on purpose.
- [ ] 3.6 Drop `debug-pipeline` if nothing but `lm_counters` needed it (2.5).
- Verify: G-1.

### Phase 4: the Analysis-mode entry (rustafits, then Athenaeum)

- [ ] 4.1 The XISF reader folds `XBAYROFF` / `YBAYROFF` from its FITS keywords, as #7 did for FITS. Athenaeum's own
  XISF fold must go in the same submodule bump: two folds would cancel on an odd offset.
- [ ] 4.2 `star_measure::raw` from `input.rs`. Header facts come from `ImageConverter::read_header`; the domain
  rules (PC-2, R-PC-2) and the tests are unchanged.
- [ ] 4.3 `star_measure::green` and `star_measure::trail` from `cfa.rs` and `trail.rs`.
- [ ] 4.4 `star_measure::analysis_mode`:
  - `FrameOptions { detection_threshold: 6.0, psf_model: Moffat4, max_stars: 500 }`;
  - `measure_raw_frame(_with)`;
  - `FrameSummary` with `map_measurement`'s frame arithmetic (2.3).
- [ ] 4.5 Athenaeum: `analysis/` keeps `config.rs` and the mapping from `FrameSummary` to `FrameAnalysis` /
  `StarMetric`; `input.rs`, `cfa.rs` and `trail.rs` go.
- Verify: G-1 on the analysis-mode arm (every row of the corpus identical, verdicts included), and Athenaeum's
  `analysis` tests.

### Phase 5: rustafits' own consumers

- [ ] 5.1 `rustafits analyze` on `analysis_mode` (D-7):
  - `tests/cli.rs` pins the CLI to `measure_raw_frame` + `FrameSummary` instead of `ImageAnalyzer::new().analyze()`;
  - README and CHANGELOG (**Changed**: the numbers).
- [ ] 5.2 `--annotate`: `annotate.rs` takes `StarMetrics`. Either add an adapter from `StarFit` (x, y, fwhm_x / y,
  θ, eccentricity), or keep the legacy analyzer for overlays until D-6. Positions are in file row order, so the
  render's flip applies, as in Blink.
- [ ] 5.3 `rustafits-debug`, `docs/usage.md`, `docs/athenaeum-integration.md`, CLAUDE.md.
- [ ] 5.4 `#[deprecated]` on `ImageAnalyzer::analyze*`, but not on `detect_fast*` (D-6).

### Phase 6 (later, 2.0): remove the legacy full analysis

This waits until the lab retires `legacy-analysis` and Athenaeum retires `SeedOverride::Full`. It removes:
- `ImageAnalyzer::analyze`, `analyze_data`, `analyze_raw` and `analyze_batch`;
- everything only they reach (a reachability pass lists it first).

`detect_fast*`, `fitting`, `background`, `convolution` and `detection` stay.

## 7. Risks

- **R-1 Bit drift across the crate boundary.** Rust neither contracts nor reassociates float operations without an
  explicit `mul_add` or intrinsic. So moving code between crates built with the same profile cannot change a result
  by itself; Athenaeum builds both at opt-level 2 in dev and with thin LTO in release. G-1 is the proof, not this
  argument.
- **R-2 The freeze collides with Athenaeum's active plans.** Phase 0 waits for them, or every later change lands
  twice.
- **R-3 Semver.** A large public surface shaped by one consumer (D-3).
- **R-4 Wire formats.** Drift in serde attributes would change stored stacking configs and artifacts. Athenaeum's
  round-trip tests catch it (G-1).
- **R-5 Log targets change** (3.5).
- **R-6 Platform pins.** libm differs by platform: x86_64 vs aarch64 (CLAUDE.md), and macOS vs MSVC `pow`.
  Athenaeum's `docs/superpowers/open-items.md` (2026-10-18) records rustafits' `byte_identical_to_before` pins and
  one athenaeum-core pin failing on Windows. rustafits 1.6.0 relaxed its pins off aarch64. Re-check that entry
  before reading a Windows failure as a port regression.
- **R-7 Build cost.** The CLI compiles about 12k lines it does not call until phase 5. Measure it at 2.x.

## 8. Verification commands

rustafits:

```bash
cargo test
cargo test --features debug-pipeline
cargo build --no-default-features
cargo clippy --all-targets --message-format=short 2>&1 | grep -E "src/(stats|star_measure)/"
```

Athenaeum (from its CLAUDE.md and Plan C's constraints):

```bash
CARGO_TARGET_DIR=.superpowers/target-acc CARGO_BUILD_JOBS=4 cargo test -p athenaeum-core --all-targets -- --test-threads=4
CARGO_TARGET_DIR=.superpowers/target-acc cargo check -p athenaeum-core --no-default-features
```

The lab: athenaeum-lab `star-lab` `star_probe` scorecard, `before` vs `after`. "0 moved" is the bar.
