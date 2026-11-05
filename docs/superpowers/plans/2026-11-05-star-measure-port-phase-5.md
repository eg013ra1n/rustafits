# star_measure port, phase 5: rustafits' own consumers — task plan

> **For agentic workers:** REQUIRED SUB-SKILL: Athenaeum's `lean-execution` skill (it overrides
> superpowers:subagent-driven-development where they conflict). One fresh implementer per task; the controller re-gates
> every task. Steps use checkbox (`- [ ]`) syntax.
>
> **Status:** approved 2026-11-05 (r2; plan review r1 APPROVE WITH FIXES, 3 blocking / 7 important / 8 minor, folded).
> The owner: «Да, оба пункта» — Q-1 yes (T4 is in), Q-2 yes (PR 13 merges before the 1.8.0 tag).

**Goal:** the `rustafits` CLI measures stars the way Athenaeum's Analysis tab does (`star_measure::analysis_mode`), its
`--annotate` overlay draws those same stars, the legacy `ImageAnalyzer::analyze*` is marked `#[deprecated]`, and
1.8.0 ships phases 4 and 5 together.

**Architecture:** `analyze` and `convert --annotate` call `raw::read_raw_frame` → `measure_raw_frame` →
`FrameSummary::from_measurement`, the stable entry phase 4 built. The annotation API gains a neutral star list
(`StarOverlay`) that both the new measurement and the legacy `AnalysisResult` convert into, so the drawing code is shared
and the legacy functions keep their signatures. No float code in `star_measure`, `analysis::fitting` or
`analysis::background` changes.

**Tech stack:** Rust 2021, clap 4 (the `cli` feature), the crate's own `star_measure`, `annotate`, `formats`.

**Spec:** the port plan `docs/superpowers/plans/2026-11-03-star-measure-port.md`: §4 decisions D-3 (API surface), D-6
(legacy analyzer), D-7 (CLI output), D-8 (version), D-10 (branching); §5 gates; §6 phase 5 items 5.1–5.4. Facts were
read at rustafits `main` 8b25b85 (= Athenaeum `main` b3cbeb85's pin) on 2026-11-05.

## Global constraints

- Numbers: `MEASURE_VERSION` stays 8; nothing under `src/star_measure/`, `src/stats.rs`, `src/analysis/{fitting,background}.rs`
  changes (G-4). The port digest stays byte-identical (`a07688d0…` on aarch64).
- Semver: 1.8.0 is additive for library callers. The legacy `ImageAnalyzer::analyze`, `analyze_data`, `analyze_raw`,
  `analyze_batch` stay and get `#[deprecated(since = "1.8.0", note = …)]`; `detect_fast*` and `FastStar` are not
  deprecated (D-6). The CLI's numbers change by design (CHANGELOG **Changed**).
- G-3 on every task: `cargo test -p rustafits`, `--features serde`, `--features debug-pipeline`,
  `cargo build -p rustafits --no-default-features` pass with **no new warnings** (deprecation warnings included);
  `cargo clippy -p rustafits` prints the same 11 lines as `main`.
- Commits: author `eg013ra1n <vilen.sharifov@gmail.com>` with the harness's attribution trailer; check
  `git log -1 --format=%an` after every commit. Never name another program in code, docs or messages (the comparison
  program is "the ERP").
- Builds: from the Athenaeum worktree root (`ATH`, see Execution), `CARGO_TARGET_DIR=$ATH/../../target-agent
  CARGO_BUILD_JOBS=4`, one compile at a time (16 GB Mac). Paths below are RF-relative (`src/…`, `tests/…`); `ATH/…` is
  Athenaeum.

## Review focus

1. **A frame with a verdict** (trailed, defocused, no stars): its medians are 0.0 (PC-7). The table must print `-`,
   not `0.00`, and `--sort fwhm` must not rank it best. → T1 tests `a_frame_with_a_verdict_prints_dashes_and_sorts_last`.
2. **Old scripts** pass `--sigma 5`, `--measure-cap 0`, `--no-debayer`: they must still run (exit 0, one stderr warning
   per flag, three lines), not die on "unexpected argument". → T1 `old_flags_still_run_with_a_warning`, T3 the convert
   half.
3. **A bottom-up FITS with `--annotate`**: ellipses must land on the stars after the render's flip. → T2
   `overlay_flips_like_the_legacy_path`, T3 `annotate_draws_on_the_stars`.
4. **An OSC mosaic**: `analyze` measures green sites at native size; the overlay scales native positions onto the
   half-size super-pixel render. → T1 `analyze_measures_a_mosaic_on_its_green_sites`, T3 `annotate_scales_a_mosaic`.
5. **A batch with one unreadable file**: exit 1, the error on stderr with the path, every other row printed (kept:
   `analyze_reports_a_broken_file_and_keeps_the_rest`; T3 adds the `--annotate` twin).

## Decisions this plan takes (review them)

- **P5-1 (5.2): an adapter, not the legacy analyzer.** The overlay must show the stars the numbers were computed on,
  and the legacy analyzer goes in 2.0 anyway; keeping it for overlays would mean doing 5.2 twice. `StarFit` already
  carries the legacy geometry convention (`fwhm_x` = major axis, `theta` = the major axis's angle from image x,
  `src/star_measure/fit.rs:278-287`; `StarMetrics`, `src/analysis/mod.rs:95-110`), so the adapter only converts f64 → f32.
  Every accepted fit is drawn (`pm.stars`, saturated and extended included), as Athenaeum's Analysis star rows store
  them (`ATH/crates/athenaeum-core/src/analysis/analyzer.rs:80-101`).
- **P5-2: flag compatibility.** `--sigma` stays as a hidden alias of `--detection-threshold`, but its value now means
  another detector's σ (the matched filter, default 6.0; the legacy DAOFIND default was 5), so it warns once:
  `warning: --sigma is --detection-threshold since 1.8.0: another detector, default 6.0`. `--measure-cap` (both
  commands) and `analyze --no-debayer` are accepted, hidden from `--help`, ignored, and print one
  `warning: … has no effect since 1.8.0: …` line on stderr. All three go in 2.0. Ranges follow Athenaeum's
  `AnalysisConfig::validate` (`ATH/…/analysis/config.rs:95-107`): threshold `MIN_THRESHOLD..=MAX_THRESHOLD`
  (2–20, `star_measure::seeds`), stars 10–2000, in both commands; `convert --max-stars` (was `1..`) is validated with
  or without `-a`, and the CHANGELOG names it.
- **P5-3: output.** The table follows D-7 (VERDICT replaces TRAIL). CSV / JSON keep the 19 current columns in their
  order, with `FrameSummary`'s values as they are (0.0 stays 0.0: what Athenaeum stores), and append six:
  `channels, stars_refused, frame_verdict, measure_version, detection_threshold, psf_model`. Appending keeps
  position-based scripts working. The columns not in `FrameSummary` follow Athenaeum's `map_measurement`
  (`ATH/…/analysis/analyzer.rs:103-131`): `width`, `height`, `channels` are the `RawFrame`'s (`src/star_measure/raw.rs:50-55`);
  `possibly_trailed` is `pm.trailed` (not `verdict == Trailed`); `stars_measured` is `pm.stars.len()` (Athenaeum's
  star rows).
- **P5-4: `-j` default stays 2** (`DEFAULT_MAX_JOBS_ANALYSIS`) unless T3 measures a peak RSS above 2.25 GB (the legacy
  1.5 GB × 1.5) at `-j 2` on 26 MP frames; then it becomes 1 and the comment says why. T3 rewrites the comment's
  numbers either way (`src/bin/rustafits/convert.rs:24-31`).
- **P5-5: `AnnotationConfig` defaults unchanged** (ecc 0.5 / 0.6). The module's eccentricity reads 5–10 % higher, but
  re-picking the thresholds is one look at the corpus for both this overlay and Athenaeum's Blink (Athenaeum
  `docs/superpowers/open-items.md`, "The Blink ellipse colour thresholds"), not part of the port.
- **P5-6: `rustafits-debug` stays on the legacy pipeline.** It diagnoses that pipeline's stages; its two `analyze*`
  calls (`src/bin/debug.rs:979, :1267`) take `#[allow(deprecated)]`, and its usage text says what it diagnoses. It is
  rewritten or removed with the pipeline at 2.0.
- **P5-7 (only if Q-1 is yes): the FITS reader reads `XBAYROFF` / `YBAYROFF` by Athenaeum's catalog rule.** The
  offsets come from `formats::stored_keys::FitsHeader::get_i32` over the same header bytes, which is what Athenaeum's
  catalog calls (`ATH/…/fits_parser/mod.rs:478`), so the render, the Analysis tab and the catalog agree by
  construction, as the XISF reader has since phase 4 (R-P4-C1). `BAYERPAT` and every other card keep the reader's
  rules. Stacking is not affected: it reads the catalog's pattern (`ATH/…/stacking/run.rs:12698-12715`). Two
  consequences, stated here so nobody later "fixes" them: (a) a repeated card now gives the FIRST value (the catalog's
  `or_insert`, `fits_header_reader.rs:123, :128`), where the reader kept the last; (b) FITS and XISF still differ on a
  fractional offset by design: `get_i32` truncates (`fits_header_reader.rs:160-166`, `0.6` → 0) and the XISF rule
  rounds (`src/formats/xisf.rs:352`, `0.6` → 1), each being Athenaeum's rule for that format.
- **P5-8: branching** continues D-10: create a local `star-measure-port` from `origin/star-measure-port` (4f3710b, the
  second parent of 8b25b85; the local `main` ref is stale, use `origin/main`), fast-forward it to 8b25b85, cut
  `feat/star-measure-phase-5` from it, one PR into `star-measure-port` (squash), then `star-measure-port` → `main` with
  a merge commit, then the 1.8.0 tag on `main`. Athenaeum pins only the tag.

## Questions for the owner (answered 2026-11-05: both yes)

- **Q-1. Fix the quoted `XBAYROFF= '1'` in this phase (T4)?** Today the FITS reader takes the integer prefix of the
  value text, so `'1'` folds as 0, while Athenaeum's catalog reads 1; it also keeps the last of a repeated card and
  matches the name case-sensitively. A mosaic written that way renders and is plate-solved with the wrong Bayer cell,
  and its Analysis row is measured on the wrong sites. The fix moves only such files' renders and Analysis rows
  (MEASURE_VERSION stays 8; a stored row of such a file stays until it is re-analyzed); a repeated card then gives the
  first value, as the catalog does. **Recommendation: yes**, as T4, gated by a corpus scan (G5-e).
- **Q-2. Merge rustafits PR 13 (the release workflow on `gh`, notes from CHANGELOG) before the 1.8.0 tag?**
  **Recommendation: yes**: 1.8.0 is then the first release whose title and notes need no hand edit.

## Execution

- Worktree: Athenaeum `.claude/worktrees/smport-p5` (`ATH`) on branch `star-measure-port-phase-5` from `main`
  b3cbeb85; inside it, `rustafits/` (`RF`) checks out `feat/star-measure-phase-5`. Register it in
  `~/.athenaeum-acc/REGISTRY.tsv`; remove it after the merges (disk rule, `lean-execution` §7).
- Models: implementers `rust-engineer` on sonnet (T1–T5, T7), T6 on sonnet; task reviewers opus; the final
  whole-branch review fable, which also hunts plan-verbatim seams (strings shared between T1 and T3, the flag names).
  Each handoff states one job and a budget of ≤ 120 tool calls; reports ≤ 80 lines in
  `ATH/.superpowers/sdd/2026-11-05-smport-p5/`.
- Order: T1 → T2 → T3 → (T4) → T5 → T6 in rustafits, the gates G5-a…G5-g, the PRs and merges, the release, then T7 in
  Athenaeum. T2 before T3 (T3 consumes T2's API); T5 after T1–T4 so their code is already off the deprecated calls.

---

## Tasks

### T1 (rustafits) — 5.1: `rustafits analyze` on `analysis_mode`

**Files:**
- Modify: `src/bin/rustafits/analyze.rs` (whole module), `src/bin/rustafits/main.rs:20` (the help example's comment)
- Test: `tests/cli.rs:286-341`

**Interfaces:**
- Consumes: `astroimage::star_measure::raw::read_raw_frame(&Path) -> anyhow::Result<RawFrame>`,
  `RawFrame::drop_unmeasured_planes(&mut self)`, `analysis_mode::{FrameOptions, measure_raw_frame, FrameSummary,
  Verdict}`, `star_measure::{PsfModel, MEASURE_VERSION, seeds::{MIN_THRESHOLD, MAX_THRESHOLD}}`.
  `FrameOptions` is `#[non_exhaustive]`: start from `FrameOptions::default()` and assign its public fields; a struct
  literal does not compile outside the crate.
- Produces (T3 reuses them, so they live in a small shared place, `src/bin/rustafits/measure.rs`, `pub(crate)`):
  - `struct Measured { frame_w: usize, frame_h: usize, channels: usize, flip_vertical: bool, pm: PlaneMeasurement, summary: FrameSummary }`
  - `fn measure_file(path: &Path, opts: &FrameOptions, pool: &Arc<ThreadPool>) -> anyhow::Result<Measured>`:
    read, `drop_unmeasured_planes`, measure with `Some(pool)`, summarize. (T3 needs the `RawFrame`'s size and flip only,
    so `Measured` keeps those, not the planes.)
  - `fn budget_parser(s: &str) -> Result<usize, String>` (10–2000) and `fn threshold_parser(s: &str) -> Result<f64, String>`
    (`MIN_THRESHOLD..=MAX_THRESHOLD`), messages `"must be between 10 and 2000"` / `"must be between 2 and 20"`.
  - `fn warn_ignored(flag: &str, why: &str)`: prints `warning: {flag} has no effect since 1.8.0: {why}` to stderr.
    T1 and T3 print these strings through this one function, never a second copy (a plan-verbatim seam).

- [x] **Step 1: rewrite the analyze test and add five (failing).** In `tests/cli.rs` (an integration-test crate: it
  cannot see the binary's `Measured`), the library reference is a helper
  `lib_measure(path: &Path, opts: &FrameOptions) -> (RawFrame, PlaneMeasurement, FrameSummary)` built from
  `read_raw_frame` + `measure_raw_frame(…, None)` + `FrameSummary::from_measurement` (G-6: no result reads a thread
  count, so `None` equals the CLI's pool).
  - `analyze_reports_exactly_what_the_library_measures` (rewrite): for `sharp.fits` / `soft.fits`, assert the CSV's
    `stars_detected`, `stars_measured` (= `pm.stars.len()`), `stars_refused`, `median_fwhm`, `median_hfr`,
    `median_eccentricity`, `median_snr`, `snr_weight`, `psf_signal`, `frame_snr`, `background`, `noise`,
    `trail_r_squared` equal the library's `f64::to_string()`; `median_beta` is empty; `frame_verdict` is empty;
    `measure_version` is `8`; `psf_model` is `moffat4`; the header's first 19 names are the current ones, in order, and
    the last six are P5-3's. Keep the sort, reverse, table and JSON assertions; the table header ends in `VERDICT`.
  - `analyze_options_reach_the_measurement`: `--psf-model auto` gives `median_beta` equal to the library's under
    `PsfModel::Auto`, `Some` or `None` (it is `Some` only for `BetaSource::Fitted`, `analysis_mode.rs:198-200`; if it is
    `None` on this fixture, R5-1 applies); `psf_model` reads `auto`; `--detection-threshold 12` gives the library's
    `stars_detected` at 12.0; `--max-stars 20` gives `stars_measured <= 20` equal to the library's at budget 20.
  - `analyze_measures_a_mosaic_on_its_green_sites`: `write_star_field`'s pixels with `BAYERPAT= 'RGGB'` added (a
    Gaussian star sampled by every site; the measurement reads the green ones) → the CSV equals the library on the same
    file; `width` / `height` are the native 320 / 240, `channels` 1.
  - `a_frame_with_a_verdict_prints_dashes_and_sorts_last`: a flat noise frame (no stars) next to `sharp.fits` → CSV
    `frame_verdict` = `no_stars`, `median_fwhm` = `0`; the table's row for it shows `-` in FWHM, HFR, ECC, SNR, WEIGHT
    and `no stars` in VERDICT; `--sort fwhm` lists `sharp.fits` first with either input order.
  - `old_flags_still_run_with_a_warning`: `analyze sharp.fits --sigma 8 --measure-cap 0 --no-debayer --csv --quiet`
    exits 0, its rows equal `--detection-threshold 8 --csv --quiet`'s, and stderr holds exactly three lines: P5-2's
    `--sigma` warning and the two `has no effect since 1.8.0` lines.
  - `out_of_range_options_are_refused`: `--max-stars 5`, `--detection-threshold 1.5`, `--detection-threshold 25` exit
    2 with the range message; nothing on stdout.
  - `analyze_matches_the_library_on_the_test_frames` (`#[ignore]`, gate G5-f: real 26 MP frames are slow in a test
    run): every CSV column of `tests/{mono,osc,cocoon}.fits` and `tests/test.xisf` equals `lib_measure`'s.
- [x] **Step 2: run them, see them fail.** `cargo test -p rustafits --test cli analyze` → FAIL (the CLI still runs
  `ImageAnalyzer`).
- [x] **Step 3: implement.**
  - `AnalyzeArgs`: `--detection-threshold <SIGMA>` (default 6.0, `threshold_parser`), `--psf-model <moffat4|auto>`
    (a `ValueEnum` mapped to `PsfModel`, default `moffat4`), `--max-stars <N>` (default 500, `budget_parser`). A clap
    `alias` cannot tell that it was used, so `--sigma <SIGMA>` is its own hidden argument (`threshold_parser`,
    `conflicts_with = "detection_threshold"`) that supplies the threshold and prints P5-2's warning once. Hidden
    `--measure-cap <N>` and `--no-debayer` only trigger `warn_ignored` once per run. `--focal-length` / `--pixel-size`
    stay.
  - Rows are `(String, Measured)`. Arcsec = `px × pixel_size / focal_length × 206.265` (the legacy formula,
    `src/analysis/mod.rs:1194-1196`), `None` without optics or with a verdict.
  - Table, CSV, JSON per P5-3; VERDICT texts `trailed`, `defocused`, `no stars`, `-`; `frame_verdict` texts `trailed`,
    `defocused`, `no_stars` (Athenaeum's `FrameVerdict`, snake_case), empty / `null` without.
  - Sort: same keys; a row with a verdict sorts after every row without one, in input order (it replaces today's NaN
    rule); `--reverse` reverses the whole list, as today.
  - The module doc says the numbers are Athenaeum's Analysis numbers (`analysis_mode`, defaults 6.0 σ, Moffat 4, 500).
- [x] **Step 4: pass.** `cargo test -p rustafits --test cli` → PASS; then G-3.
- [x] **Step 5: commit** `feat(cli): analyze on star_measure::analysis_mode, Athenaeum's Analysis numbers (port phase 5.1)`.

### T2 (rustafits) — 5.2, library half: `StarOverlay`

**Files:**
- Modify: `src/annotate.rs` (new types and functions above Tier 1; the three legacy functions become wrappers),
  `src/lib.rs:26-28` (re-exports)
- Test: `src/annotate.rs` test module

**Interfaces:**
- Produces (public, re-exported at the crate root):
  - `pub struct OverlayStar { pub x: f32, pub y: f32, pub fwhm_x: f32, pub fwhm_y: f32, pub theta: f32, pub eccentricity: f32, pub fwhm: f32 }`
    (`Debug, Clone, Copy, PartialEq`): source pixels, file row order, `fwhm_x` the major axis.
  - `pub struct StarOverlay { pub stars: Vec<OverlayStar>, pub width: usize, pub height: usize, pub median_fwhm: f32 }`
    (`Debug, Clone, PartialEq`). Both are plain, all-public structs (callers may build them from their own
    measurements), so their doc says the field set is frozen for 1.x.
  - `StarAnnotation` (`src/annotate.rs:60-61`, derives nothing today) gains `#[derive(Debug, Clone, Copy, PartialEq)]`
    (additive), so the tests can compare it.
  - `impl StarOverlay { pub fn from_measurement(frame_w: usize, frame_h: usize, pm: &PlaneMeasurement, summary: &FrameSummary) -> StarOverlay }`:
    one `OverlayStar` per `pm.stars` entry (`fwhm()` and `eccentricity()` from `StarFit`), `median_fwhm =
    summary.median_fwhm as f32`.
  - `impl From<&AnalysisResult> for StarOverlay`: the legacy fields as they are.
  - `pub fn compute_overlay(overlay: &StarOverlay, output_width: usize, output_height: usize, flip_vertical: bool, config: &AnnotationConfig) -> Vec<StarAnnotation>`
  - `pub fn create_overlay_layer(overlay: &StarOverlay, output_width: usize, output_height: usize, flip_vertical: bool, config: &AnnotationConfig) -> Vec<u8>`
  - `pub fn annotate_image_with(image: &mut ProcessedImage, overlay: &StarOverlay, config: &AnnotationConfig)`
  - `compute_annotations`, `create_annotation_layer`, `annotate_image` keep their signatures and call the new three
    through `StarOverlay::from(result)`. `compute_overlay` keeps today's empty-result guard (`stars` empty or a zero
    size → empty).

- [x] **Step 1: tests (failing).**
  - `legacy_functions_draw_what_they_drew`: on three `make_star` stars of a `dummy_result` (one round; one elongated,
    its `theta` set to 0.6 after `make_star`, which fixes 0.0 at `:494`; one past `ecc_warn`), a source 200×100 drawn
    at 100×50, `compute_annotations` returns, field by field, the values computed by hand from today's rules
    (`x · scale`; `h − 1 − y · scale` and `−θ` with the flip; semi-axes `fwhm · scale · ellipse_scale` clamped to
    `[min_radius, max_radius]`; the colour from `ecc_good` / `ecc_warn`), flip on and off, and equals
    `compute_overlay(&StarOverlay::from(&r), …)`. The drawing primitives do not change, so this pins the refactor.
  - `overlay_flips_like_the_legacy_path`: one star at (10, 5) on a 100×50 source, output 100×50: `flip_vertical =
    true` gives y = 44 and `theta` negated; `false` gives y = 5.
  - `from_measurement_keeps_every_fit`: a `PlaneMeasurement` from `measure_raw_frame` on a `RawFrame` built in the test
    as a struct literal (its fields are public: `RawPlanes::Mono(moffat_field(…))`, `Domain::RawAdu { … }` as
    `analysis_mode`'s tests write it); `analysis_mode`'s own `mono` helper is private to its tests and
    `test_support.rs` stays byte-identical to Athenaeum's → `stars.len() == pm.stars.len()`, each `OverlayStar`'s
    `x, y, fwhm_x, fwhm_y, theta` equal the fit's as f32, `fwhm_x >= fwhm_y`.
  - `median_zero_colours_green` (a pin: `annotate.rs:426-428` already does it): `ColorScheme::Fwhm` with
    `median_fwhm = 0.0` gives `[0, 255, 0]` (a frame whose only fits are saturated).
- [x] **Step 2: fail.** `cargo test -p rustafits annotate` → FAIL (types missing).
- [x] **Step 3: implement** per Interfaces; the drawing primitives stay as they are.
- [x] **Step 4: pass,** then G-3.
- [x] **Step 5: commit** `feat(annotate): StarOverlay, a star list both the measurement and AnalysisResult draw through (port phase 5.2)`.

### T3 (rustafits) — 5.2, CLI half: `convert --annotate` on `analysis_mode`

**Files:**
- Modify: `src/bin/rustafits/convert.rs:24-31` (the jobs comment), `:94-103` (flags), `:336-360` (the verbose line),
  `:381-433` (the annotate path)
- Test: `tests/cli.rs`

**Interfaces:**
- Consumes: T1's `measure_file`, `budget_parser`, `warn_ignored`; T2's `StarOverlay::from_measurement`,
  `annotate_image_with`.

- [x] **Step 1: tests (failing).**
  An ellipse is a ring whose semi-axes clamp to ≥ `min_radius` 6 output px (`src/annotate.rs:49, :112-113`; the
  fixture's FWHM 2.8 is far below it), so the tests look for the ring, not the centre:
  - `annotate_draws_on_the_stars`: `write_star_field` sharp → `rustafits sharp.fits a.png -a` and `plain.png`; the
    annotated PNG holds ≥ 20 pixels exactly `[0, 255, 0]`, the plain one none; the star at file (12, 12) lands at
    output (12, 227) once flipped (bottom-up: no `ROWORDER`), and a green pixel lies within ±1 px of (18, 227).
  - `annotate_scales_a_mosaic`: T1's RGGB field → the output is 160×120 (super-pixel); native (12, 12) → output (6, 6)
    → flipped (6, 113); green pixels lie within ±1 px of (12, 113) and of (6, 107).
  - `annotate_old_flags_still_run`: `-a --measure-cap 0` exits 0 with one warning; `-a --max-stars 5` exits 2.
  - `annotate_reports_a_broken_file_and_keeps_the_rest`: one good and one broken file with `-a` → exit 1, the good
    output written.
- [x] **Step 2: fail.**
- [x] **Step 3: implement.** With `--annotate`: `measure_file` first (budget = `--max-stars`, other options
  default); it returns without the planes, so the render's `read_raw` starts after they are freed and the two reads
  never overlap in memory (the file is read twice: the measurement's reader states the ADU facts, the render's does
  not; the measurement dominates the time); the overlay is
  `StarOverlay::from_measurement(m.frame_w, m.frame_h, &m.pm, &m.summary)`; the verbose line keeps its shape
  (`analysis {ms} ms ({stars_detected} stars, FWHM {median_fwhm:.2})`), with the measurement's own file read inside
  `analysis` and `read` now the render's read only (say so in the line's doc comment). Measure and rewrite the jobs
  comment: wall time and peak RSS (`/usr/bin/time -l`) of `convert -a` and `analyze` over four copies of
  `tests/cocoon.fits` at `-j 1` and `-j 2`, release build, and the same `analyze` runs on `main`'s binary, so the
  CHANGELOG can say how the speed moved; P5-4's bar decides the default.
- [x] **Step 4: pass,** then G-3.
- [x] **Step 5: commit** `feat(cli): --annotate draws the analysis_mode measurement (port phase 5.2)`.

### T4 (rustafits, only if Q-1 is yes) — the FITS reader's Bayer offsets by the catalog rule

**Files:**
- Modify: `src/formats/fits.rs:36-44, :96-125` (the card loop's `XBAYROFF` / `YBAYROFF` arms), and the fold at `:161-163`
  only if its input changes type
- Test: `src/formats/fits.rs` tests, `tests/cli.rs` (`info`)

**Interfaces:**
- Consumes: `formats::stored_keys::FitsHeader::{from_bytes, get_i32}` (crate-private, `src/formats/stored_keys/fits_header_reader.rs:138, :160`).

- [x] **Step 1: tests (failing).** `quoted_bayer_offsets_fold_like_the_catalog`: `BAYERPAT= 'RGGB'`, `XBAYROFF= '1'`
  → pattern GRBG; `YBAYROFF= '1'` alone → GBRG; both → BGGR. `bayer_offsets_follow_the_catalog_rule_on_unusual_cards`:
  a repeated `XBAYROFF` (0 then 1), a lower-case `xbayroff = 1`, and `XBAYROFF= 1.0` give exactly what
  `FitsHeader::from_bytes(…).get_i32("XBAYROFF")` gives — computed in the test, not hand-written, so the test pins
  parity. Plain `XBAYROFF= 1` keeps today's result (a regression line).
- [x] **Step 2: fail.**
- [x] **Step 3: implement:** `read_fits_header` (`:65-130`) reads 2880-byte blocks without keeping them; it now also
  appends each block to a `Vec<u8>` and, after `END`, takes both offsets from `stored_keys::FitsHeader::from_bytes(…)?
  .get_i32(…)` (absent → 0; import it under another name: the reader has its own `FitsHeader` struct). If
  `from_bytes` fails where the reader's own parse succeeded, keep the reader's integer-prefix values and log a
  `tracing::warn!` with the error (never silently). Every other card keeps the reader's own parse.
- [x] **Step 4: pass,** G-3, and G5-b's render identity (none of the four test files carries an offset).
- [x] **Step 5: commit** `fix(fits): XBAYROFF / YBAYROFF read by Athenaeum's catalog rule (quoted, repeated, any case)`.

### T5 (rustafits) — 5.4: `#[deprecated]` on `ImageAnalyzer::analyze*`

**Files:**
- Modify: `src/analysis/mod.rs:499, :513, :532, :548` (the attribute and doc lines; `analyze_batch` calls
  `analyze_impl`, `:572`, and the module's tests call no `analyze*`, so nothing inside the module needs an allow),
  `src/bin/debug.rs:19, :966-979, :1255-1267` and its usage text `:50-62`, `tests/analysis.rs`,
  `tests/fast_detect_real.rs`, and the three examples `cargo test` compiles (no `autoexamples = false`):
  `examples/analyze_batch.rs:29`, `examples/trail_test.rs:9`, `examples/fastdump.rs:9`.

- [x] **Step 1:** add `#[deprecated(since = "1.8.0", note = "use star_measure::analysis_mode (raw::read_raw_frame, measure_raw_frame, FrameSummary); removed in 2.0")]`
  to the four methods; not to `detect_fast*`, `with_*`, `AnalysisResult` or `StarMetrics`.
- [x] **Step 2:** build every G-3 target and list each deprecation warning. Silence them where the legacy pipeline is
  the subject: `#![allow(deprecated)]` at the top of `tests/analysis.rs`, `tests/fast_detect_real.rs`,
  `examples/trail_test.rs` and `examples/fastdump.rs` (each with one comment line: it exercises the legacy analyzer,
  removed in 2.0), and an item-level `#[allow(deprecated)]` on debug.rs's `cmd_pipeline` / `analyze_one_file` (P5-6).
  `examples/analyze_batch.rs` is the library's batch showcase, so it is ported instead: rayon over the paths calling
  `read_raw_frame` → `measure_raw_frame` → `FrameSummary::from_measurement`, printing the same per-file line. Any other
  caller is a defect of T1 / T3, reported, not silenced.
- [x] **Step 3:** `rustafits-debug`'s usage text gains one line: it diagnoses the legacy `ImageAnalyzer` pipeline
  (deprecated since 1.8.0); the measurement Athenaeum and `rustafits analyze` use is `star_measure`.
- [x] **Step 4:** G-3 with zero `deprecated` warnings (`cargo build … 2>&1 | grep -c 'use of deprecated'` → 0 on every
  target).
- [x] **Step 5: commit** `chore(analysis): deprecate ImageAnalyzer::analyze* for star_measure::analysis_mode (port phase 5.4)`.

### T6 (rustafits) — 5.3: docs, CHANGELOG, CLAUDE.md

**Files:** `README.md:16, :76, :103, :130-160, :235-380`; `docs/usage.md`; `docs/annotation.md`;
`docs/athenaeum-integration.md` §2–3 (`:81-360`); `docs/analysis-pipeline.md:1-10`; `CLAUDE.md` (Crate Structure,
Downstream Users); `CHANGELOG.md` `[Unreleased]`; the port plan's §6 phase 5 (status).

- [x] **Step 1: README.** The `analyze` section: what it measures (Athenaeum's Analysis numbers: 6.0 σ, Moffat 4, 500
  stars), the new flags, the VERDICT column, the deprecated flags, and the example table regenerated from the T1
  release build on `tests/mono.fits tests/osc.fits tests/test.xisf --sort fwhm`. The library section leads with
  `star_measure::analysis_mode`; the `ImageAnalyzer` examples move under "Legacy analyzer (deprecated since 1.8.0,
  removed in 2.0)"; fast detection stays where it is; the annotation example uses `StarOverlay` +
  `annotate_image_with`; the dependency snippet's `version = "1.6"` (`README.md:174`) becomes `"1.8"`.
- [x] **Step 2: docs/.** `usage.md`: a new first section for `analysis_mode` (read, measure, summarize, the overlay),
  the rest under the legacy heading. `annotation.md`: `StarOverlay` and the three `*_overlay` / `_with` functions, the
  legacy functions as wrappers, the CLI's `--max-stars` budget. `athenaeum-integration.md` §2: Athenaeum measures
  through `pub use astroimage::star_measure` (stacking, LN, Analysis, Blink, plate-solve gate), not `ImageAnalyzer`;
  §3: the overlay types. `analysis-pipeline.md`: a two-line banner (legacy pipeline, deprecated, see `star_measure`).
- [x] **Step 3: CLAUDE.md.** Crate Structure: `analyze` runs `analysis_mode` and `tests/cli.rs` pins it to
  `measure_raw_frame` + `FrameSummary`; drop "These are no longer the numbers Athenaeum shows". Downstream Users: the
  two stale paragraphs (`ImageAnalyzer::analyze*` "no longer called by the app…", and the port-plan sentence) say
  deprecated since 1.8.0, removed in 2.0, callers left: Athenaeum's probe arm and lab, `rustafits-debug`, the
  `trail_test` / `fastdump` examples and the legacy tests.
- [x] **Step 4: CHANGELOG `[Unreleased]`** (it already holds phase 4): **Changed** — `rustafits analyze` and
  `--annotate` measure with `star_measure::analysis_mode`, so the numbers change (they are Athenaeum's Analysis
  numbers); VERDICT replaces TRAIL; six CSV / JSON columns appended; `--max-stars` 10–2000. **Added** — `StarOverlay`,
  `OverlayStar`, `compute_overlay`, `create_overlay_layer`, `annotate_image_with`; `--detection-threshold`,
  `--psf-model`. **Deprecated** — `ImageAnalyzer::analyze*`; the CLI's `--sigma` (alias), `--measure-cap`,
  `analyze --no-debayer`; `convert --max-stars` below 10 is refused. **Fixed** (if T4) — FITS `XBAYROFF` /
  `YBAYROFF` read by Athenaeum's catalog rule: a quoted value, any case, the first of a repeated card; a fractional value
  still truncates in FITS and rounds in XISF, by design (P5-7).
- [x] **Step 5:** the port plan's phase 5 items 5.1–5.4 ticked with this plan's path; commit
  `docs: phase 5 — analyze and --annotate on analysis_mode, the overlay API, the deprecation (port phase 5.3)`.

### T7 (Athenaeum, after the 1.8.0 tag) — the bump and the deprecation's callers

**Files:** `rustafits` (the submodule pin), `crates/athenaeum-core/src/stacking/measure.rs:1283-1289` (the
`SeedOverride::Full` arm), `crates/athenaeum-core/examples/{star_probe.rs:1024-1039, integrate_probe.rs:395-399}`
and any further caller the build names; `docs/stacking/README.md` (R-OSP-19, phase 5); `docs/superpowers/open-items.md`
(the quoted-`XBAYROFF` item, if T4).

- [ ] **Step 1:** pin `rustafits` at the `v1.8.0` tag (a commit on rustafits `main`).
- [ ] **Step 2:** `#[allow(deprecated)]` on each legacy call, with a comment naming D-6 (removed at rustafits 2.0, with
  `SeedOverride::Full` and the lab's `legacy-analysis` arm). `cargo build -p athenaeum-core --all-targets 2>&1 | grep -c
  'use of deprecated'` → 0.
- [ ] **Step 3:** G5-d.
- [ ] **Step 4:** R-OSP-19 gets a phase-5 paragraph (what changed for Athenaeum: nothing measured; the overlay API; the
  FITS offset fix if T4, with G5-e's result); the open item is deleted if T4 landed (a fixed item is deleted, not
  ticked).
- [ ] **Step 5:** commit `chore: pin rustafits to v1.8.0 (port phase 5), allow the legacy analyzer's last callers`.

---

## Gates (the controller, on the owner's Mac)

| gate | what | bar |
| ---- | ---- | ---- |
| G5-a | G-3 on the branch tip | all pass; zero new warnings; clippy = `main`'s 11 lines |
| G5-b | renders and `info` on `mono.fits`, `osc.fits`, `test.xisf` (`scripts/star_measure_port/` render check, as phase 4's G-1d) | byte-identical to `main` 8b25b85 (6 outputs; `analyze` changes by design) |
| G5-c | `port_digest` through the branch | byte-identical (`a07688d0…`) |
| G5-d | Athenaeum on T7: `cargo test -p athenaeum-core --all-targets`, `ts_contract`, `--no-default-features`, `athenaeum-web` build | pass with no pin edited; `MEASURE_VERSION` 8; zero `deprecated` warnings |
| G5-e (T4 only) | header scan of the star-test corpus FITS for a quoted, repeated or lower-case `XBAYROFF` / `YBAYROFF`, then athenaeum-lab `smport.sh`'s analysis-mode arm, before = b3cbeb85, after = T7 | scan lists the files (expected none); `SMPORT: PASS` with 0 moved, or exactly the scanned files moved |
| G5-f | CLI parity on real frames: T1 adds an `#[ignore]`d `analyze_matches_the_library_on_the_test_frames` to `tests/cli.rs` (the CSV of `mono.fits`, `osc.fits`, `cocoon.fits`, `test.xisf` against `lib_measure` on each); run `cargo test -p rustafits --release --test cli -- --ignored` | every column identical at full precision, `test.xisf` included |
| G5-g | the owner looks at `convert -a` PNGs of `osc.fits` (the half-scale mosaic) and `mono.fits` (bottom-up flip), release build | ellipses sit on the stars |

G5-f is the claim the phase exists for: the CLI prints `analysis_mode`'s numbers on real frames, and phase 4's G-1
proved that Athenaeum's Analysis tab prints the same ones (its rows are `FrameSummary` plus the mapping), so the CLI
prints the Analysis tab's numbers. G5-g is the only check of T3's two geometries off synthetic fields.

## Merges and the release

1. PR `feat/star-measure-phase-5` → `star-measure-port` (squash) after G5-a…c, G5-f and G5-g; PR `star-measure-port` →
   `main` (merge commit, D-10). Pushes and merges on the owner's word.
2. Q-2's PR 13, if approved, merges before the tag.
3. 1.8.0 by the `rustafits-release` skill: the `[Unreleased]` section becomes `[1.8.0]` with phases 4 and 5; tag on
   `main`; verify crates.io `max_version`, the three Linux assets and the tap commit by fact. Tag, push and publish on
   the owner's word.
4. T7 in Athenaeum, G5-d (and G5-e), PR, merge on the owner's word.

## Risks

- **R5-1 The test fixtures may not suit the new path.** `write_star_field` (320×240, Gaussian σ 1.2 / 2.2) was built
  for the legacy detector. If the module refuses it (a verdict, < 20 stars), the implementer reports; the fixture may
  change (amplitude, spacing), the module may not.
- **R5-2 Memory with `-j 2`.** The module's working set is about 9 planes per frame (Athenaeum's
  `MEASURE_PLANES_RESIDENT`), close to the legacy 1.5 GB at 2 jobs on 26 MP; T3's measurement decides whether the
  default stays.
- **R5-3 Users read old CSV by position.** Covered by appending (P5-3); the CHANGELOG names the new columns.
- **R5-4 Deprecation noise downstream.** solvemyastro uses only `detect_fast` (not deprecated); Athenaeum's three legacy
  callers are T7's. A downstream with `-D warnings` would fail on 1.8.0: none is known.
