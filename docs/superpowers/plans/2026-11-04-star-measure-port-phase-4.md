# star_measure port, phase 4: the Analysis-mode entry — task plan

> **Status:** draft r1 (2026-11-04), for review before implementation. It details items 4.1–4.5 of the port plan
> `docs/superpowers/plans/2026-11-03-star-measure-port.md` (section 6), under its decisions D-1…D-10 and gates G-1…G-6.
> Facts were read at rustafits 8ecca88 (v1.7.0, branch `feat/star-measure-phase-4`) and Athenaeum `main` e6344ddd
> (branch `star-measure-port-phase-4`). The four Athenaeum files to move are unchanged since the G-2 freeze sha
> 0b6fb183 (`git log 0b6fb183..e6344ddd -- crates/athenaeum-core/src/analysis/` is empty), so they are copied at
> e6344ddd and the commit messages name that sha (D-9).
>
> **Execution:** one fresh implementer per task, the controller re-gates every task (Athenaeum's `lean-execution`
> skill). Builds: `CARGO_TARGET_DIR=/Volumes/BigMac/Users/astrobureau/Documents/Projects/athenaeum/target-agent
> CARGO_BUILD_JOBS=4`, one compile at a time (16 GB Mac). Every command below runs from the Athenaeum worktree root
> `/Volumes/BigMac/Users/astrobureau/Documents/Projects/athenaeum/.claude/worktrees/smport-p4` (`ATH`), whose
> `rustafits/` (`RF`) is this repository; rustafits is a member of Athenaeum's workspace (Athenaeum `Cargo.toml:3-12`),
> so `-p rustafits` works there. In the text, paths are RF-relative (`src/…`, `tests/…`); `ATH/…/` is Athenaeum's
> `crates/athenaeum-core/src/`; the bare names `input.rs`, `cfa.rs`, `trail.rs`, `analyzer.rs`, `config.rs` are
> Athenaeum's `ATH/…/analysis/` files, and `fits_header_reader.rs`, `stored_header.rs` are under `ATH/…/fits_parser/`.
> In the command blocks, paths are relative to `ATH` (so rustafits files are `rustafits/…`).
>
> **Commits:** author `eg013ra1n <vilen.sharifov@gmail.com>`, trailer `Co-Authored-By: Claude Opus 5.5
> <noreply@anthropic.com>`. Phase 3's cloud commits were authored "Claude" (Athenaeum c05bf332, ca21e19c): check
> `git log -1 --format=%an` after every commit.

## 0. Summary

| what | from (Athenaeum) | to (rustafits) | change beyond paths |
| ---- | ---- | ---- | ---- |
| XISF Bayer offset fold | `ATH/crates/athenaeum-core/src/analysis/input.rs:85-109, :127-128, :240` | `src/formats/xisf.rs` (`parse_xisf_xml`) | the fold moves into the reader (4.1) |
| raw frame input | `analysis/input.rs` | `src/star_measure/raw.rs` | header facts from `ImageConverter::read_header`; the reader folds, `raw` does not |
| CFA green sites | `analysis/cfa.rs` | `src/star_measure/green.rs` | none |
| trail R² | `analysis/trail.rs` | `src/star_measure/trail.rs` | none |
| measurement + frame arithmetic | `analysis/analyzer.rs:24-92, :130-189` | `src/star_measure/analysis_mode.rs` | `AnalysisConfig` → `FrameOptions`; the arithmetic returns a `FrameSummary` |
| `AnalysisConfig`, the row mapping, ids, `config_hash`, `analyzed_at` | `analysis/{config,analyzer}.rs` | stay in Athenaeum | `map_measurement` reads `FrameSummary` |

Order: T1 → T2 → T3 → T4 → T5 in rustafits, then T6 → T7 in Athenaeum, L1 in athenaeum-lab, then the gates (section 4).
**T6 is one Athenaeum commit:** the submodule bump to the rustafits phase-4 tip and the deletion of `input.rs` land
together, never apart (section 1.2).

No number moves: `MEASURE_VERSION` stays 8, `AnalysisConfig::config_hash` is unchanged (so no stored analysis is
re-run), and every Analysis row of the star-test corpus is byte-identical (G-1, section 4).

## 1. Facts (the brief's questions 1–8, with evidence)

### 1.1 Header facts

`input.rs` reads its header facts before the pixels (`input.rs:155-156`) into `HeaderFacts` (`input.rs:65-74`):

```rust
struct HeaderFacts { bitpix: Option<i32>, bzero: f64, bscale: f64, datamax: Option<f64>, egain: Option<f64>, xbayroff: i64, ybayroff: i64 }
```

Nothing else comes from Athenaeum's parsing: width, height, channels, the Bayer pattern and `flip_vertical` come from
`ImageConverter::read_raw`'s `ImageMetadata` (`input.rs:157-159, :236, :256`). Gain is `EGAIN` only.

| fact | Athenaeum today | rustafits `ImageConverter::read_header` (`src/converter.rs:91-93`) | gap → what the plan does |
| ---- | ---- | ---- | ---- |
| FITS `BITPIX` | `FitsHeader::get_i32` (`input.rs:133`): `i64` parse, else D-exponent float, `as i32` (`ATH/…/fits_parser/fits_header_reader.rs:151-158`) | `ImageHeader::get("BITPIX")`, the card's value text (`src/formats/fits.rs:233, :240-275`; `src/types.rs:52-58`, case-insensitive, first match) | numeric rule differs → `raw` copies Athenaeum's rule as `fits_i32` |
| FITS `BZERO`, `BSCALE` (defaults 0, 1) | `get_f64` (`input.rs:134-135`), Fortran `D` exponent allowed (`fits_header_reader.rs:248-255`) | value text | → `fits_f64` = `parse_fits_f64` verbatim; defaults unchanged |
| FITS `DATAMAX`, `EGAIN` | `get_f64` (`input.rs:136-137`) | value text | → `fits_f64` |
| FITS `XBAYROFF`, `YBAYROFF` | not read: 0 (`input.rs:138-141`), the reader folds (`src/formats/fits.rs:161-163`) | `header.bayer_pattern` is the folded pattern | none; still not read |
| XISF `BITPIX`, `BZERO`, `BSCALE` | fixed `None`, 0, 1 (`input.rs:122-124`) | an XISF's `<FITSKeyword>`s are visible through `get` (`src/formats/xisf.rs:627-638`); two corpus XISF (`singles/Light_Pane 4…`, `Light_Pane 7…`) carry `BITPIX 16`, `BZERO 32768` | → `raw` keeps the fixed values for XISF and never reads these keys there |
| XISF `DATAMAX`, `EGAIN` | `<FITSKeyword>` text through `extract_xisf_header` + `parse_stored_header_keys` (`input.rs:118-120`; `ATH/…/fits_parser/mod.rs:189-193`; `stored_header.rs:72-85, :91-122`), parsed by `num` (`input.rs:111-114`): quotes trimmed, `f64` | `get` on the `<FITSKeyword>` (value unquoted, `''` → `'`, `xisf.rs:627-638`) | → `num` verbatim on `get` |
| XISF `XBAYROFF`, `YBAYROFF` | `offset(num(..))` (`input.rs:76-79, :127-128`): `f64`, finite, `round`, `i64`; absent 0; folded at `input.rs:240` | not read, not folded at 8ecca88 (`xisf.rs:312-321`) | → 4.1: the XISF reader folds with this rule; `raw` stops folding (section 1.2) |
| XISF or FITS | extension only (`input.rs:155`, `ATH/…/fits_writer/mod.rs:60-65`) | `header.format` (`src/formats/mod.rs:29-42, :44-58`: extension, else the `XISF0100` signature) | → `raw` uses `header.format`, the decision the pixel reader took; differs only for XISF bytes without a `.xisf` name, which fail today (Athenaeum's FITS header parse) and read after |

**BITPIX 16 / BZERO 0 (R-PC-1a).** The rule itself lives in the reader (`src/formats/fits.rs:308-330`: more than 1 %
negative stored samples → read with the 32768 offset, else negatives clamp to 0). `raw` only picks the full scale
for `(BITPIX 16, BZERO 0, BSCALE 1)` (`input.rs:172-176, :273-283`), so it needs exactly the three facts above,
parsed as today. Unchanged.

**Residual edge cases** (no corpus file has any; the header-only parity check below found 0 differences over the 147
corpus FITS headers):
- a card with `=` in column 9 but no space in column 10: Athenaeum skips it (`fits_header_reader.rs:101`), rustafits
  takes it (`src/formats/fits.rs:246`), and so does the pixel reader (`src/formats/fits.rs:31`). After the port the
  facts agree with the decode; today they can disagree;
- a non-UTF-8 card: rustafits drops it (`src/formats/fits.rs:90`), Athenaeum reads it lossily;
- an XISF `<FITSKeyword>` written as a start/end pair: Athenaeum reads empty elements only (`fits_parser/mod.rs:231`),
  rustafits both (`xisf.rs:618`);
- an XISF header whose keyword pass hits an XML error: Athenaeum keeps a partial result (`fits_parser/mod.rs:263`),
  `read_header` fails (`xisf.rs:616`). The layout pass bails on XML errors too (`xisf.rs:283`), so `read_raw` most
  likely fails on such a file today;
- the error text for an unreadable header changes (it reaches the Analysis tab as `"{path}: {error}"`,
  `ATH/…/api/analysis.rs:338`).

Header-only parity check (done for this plan, scratch script, not committed): each corpus FITS header parsed the
Athenaeum way and the `parse_card` + `get` way, for BITPIX, BZERO, BSCALE, DATAMAX, EGAIN, XBAYROFF, YBAYROFF and
BAYERPAT; plus odd cards and duplicates. Result: 147 headers, 0 differences, 0 odd cards, 0 duplicates. The G-1 run
covers it anyway: the analysis-mode arm prints `clip`, `fullScale` and `gainEPerAdu` (`ATH/…/examples/star_probe.rs:1241-1246`).

### 1.2 The XISF Bayer offset fold (4.1)

- **Athenaeum today** folds XISF itself, in `read_inner` (`input.rs:240`), with the offsets of `header_facts`'
  XISF arm (`input.rs:127-128`) and the parity rule of `fold_bayer_offsets` (`input.rs:85-109`, `rem_euclid(2) == 1`).
  FITS offsets are 0 there because the reader folds (`input.rs:138-141`; rustafits 1.7.0, `src/formats/fits.rs:161-163`).
- **rustafits after 4.1** folds in `parse_xisf_xml` (`src/formats/xisf.rs:237-324`), right after the pattern scan
  (`xisf.rs:312-321`). Both `read_xisf_image` (`xisf.rs:710`) and `read_xisf_header_info` (`xisf.rs:587`) call it, so
  the pixel read and the header-only read report the same folded pattern, as `layout` guarantees for FITS
  (`src/formats/fits.rs:139-177`). `shift_bayer_pattern` (`src/formats/fits.rs:182-206`) moves to `src/formats/mod.rs`
  as `pub(crate)` and both readers call it.
- **The offsets' rule** is Athenaeum's (`input.rs:76-79, :111-114`): the first `<FITSKeyword>` named `XBAYROFF` /
  `YBAYROFF` (case-insensitive), quotes trimmed, parsed as `f64`, finite, rounded; absent or unparsable is 0; odd is
  `rem_euclid(2) == 1`. It is the rule whose results the Analysis tab shows today, so no Analysis number moves, even
  for a non-integer offset. (The FITS reader's integer-prefix rule, `src/formats/fits.rs:36-44`, stays as it is.)
- **Exactly one fold, every path, every commit:**

  | read | Athenaeum e6344ddd + rustafits 8ecca88 | rustafits T1…T5 commits (Athenaeum's committed pin still 8ecca88) | Athenaeum T6 + the rustafits tip |
  | ---- | ---- | ---- | ---- |
  | Analysis, XISF | Athenaeum's fold: 1 | committed: 1. Athenaeum's working tree, whose submodule is at T1+: 2 (see below) | the reader's fold: 1; `raw` folds nothing |
  | Analysis, FITS | the reader: 1 | 1 | 1 |
  | render / `read_header` / detection, XISF | 0 (the bug 4.1 fixes) | 1 | 1 |

  `raw.rs` never contains an XISF fold: T1 precedes T2, and T2 copies `input.rs` with the fold removed. Athenaeum
  pins no rustafits commit that has T1 while `input.rs` still folds: the bump and the deletion are one commit (T6).
  Between T1 and T6 the Athenaeum worktree builds against the submodule's working copy, so its
  `analysis::input::tests::an_xisf_mosaic_folds_its_bayer_offsets` fails there **by design** (it is the tripwire).
  Athenaeum's suite is therefore not a gate before T6; T1–T5 check that Athenaeum compiles (`cargo check`).
- **Tripwires that read an odd offset end to end:** T1's reader tests (0 folds fail them); T2's
  `an_xisf_mosaic_reads_its_bayer_offsets_folded_once` (`raw` over the reader: 0 or 2 folds read RGGB, not GRBG);
  Athenaeum's `an_xisf_mosaic_folds_its_bayer_offsets`, kept with its body unchanged (T6), through Athenaeum's own XISF
  writer at the pinned submodule; and the gate's positive control (section 4, G-1c).
- **Where else rustafits reads an XISF Bayer pattern, and what 4.1 changes there** (only for an XISF mosaic with an
  odd `XBAYROFF` and/or `YBAYROFF`; files without the keywords, or with even offsets, are unchanged):
  - the render: super-pixel and VNG debayer by `meta.bayer_pattern` (`src/pipeline.rs:36-52, :110-119`). The phase
    becomes right: red and blue no longer exchanged (both offsets odd), green and magenta no longer exchanged (one odd);
  - `read_header` (`xisf.rs:601`) and so `rustafits info` (`src/bin/rustafits/info.rs:75, :126-130`): the folded pattern;
  - `ImageAnalyzer`'s green interpolation (`src/analysis/mod.rs:663-667` `detect_fast` on a path, `:692-695`
    `detect_fast_raw`, `:719-722` `analyze_raw`, `:745-748` `analyze`, `:1820-1823` `prepare_luminance`), hence
    `rustafits analyze` and `rustafits-debug` (`src/bin/debug.rs:219`): with one odd offset the true green sites are
    interpolated (with both odd the green sites are the same and nothing changes);
  - downstream: solvemyastro's `detect_fast(path)` (`ATH/solvemyastro/src/orchestrate.rs:1074`) detects on the
    corrected green sites; Athenaeum's previews, thumbnails and Blink renders (`ATH/…/rustafits_processor/mod.rs:155-187`)
    show the right colours, after a restart for an already viewed file (the preview cache is in memory,
    `ATH/…/cache/mod.rs:14-21`). Athenaeum's calibration and stacking are unaffected: they read `BAYERPAT` /
    `XBAYROFF` from the catalog (`ATH/…/calibration_library/light_cal.rs:726-756`), not from the reader.
  - Corpus: no image carries an odd offset (5 XISF, 147 FITS; only `singles/test.xisf` has the keywords, both 0).
- **CHANGELOG `[Unreleased]`, Fixed** (T1 writes it):

  > - XISF: `XBAYROFF` / `YBAYROFF` (`<FITSKeyword>`) are folded into the Bayer pattern, as they are for FITS since
  >   1.7.0: an odd `XBAYROFF` swaps the 2×2 cell's columns (RGGB → GRBG), an odd `YBAYROFF` its rows (RGGB → GBRG).
  >   A mosaic with an odd offset was debayered with the wrong phase by the super-pixel and the VNG render (red and
  >   blue exchanged, or green and magenta), and with one odd offset `ImageAnalyzer`'s green interpolation
  >   (`analyze*`, `detect_fast*`) took red / blue sites for green. `read_header` and `rustafits info` report the
  >   folded pattern. Files without the keywords, or with even offsets, read and render exactly as before.

### 1.3 `analyzer.rs`'s split

What `map_measurement` computes per frame (`analyzer.rs:130-189`) moves into rustafits as `FrameSummary`, statement
for statement and in the same order. The star rows, the ids, `config_hash`, `analyzed_at`, `quality_score`,
`detection_threshold`, `width`, `height`, `source_channels`, `possibly_trailed` and `measure_version` stay in
Athenaeum: they are row facts, not frame arithmetic.

```rust
/// The frame figures of an Analysis-mode measurement (Athenaeum's spec 4.2, rulings PC-6…PC-9, PC-11, PC-13).
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct FrameSummary {
    /// The matched filter's seeds not flagged `extended` (R-PC-5 / R-PC-8), every layout: `analyzer.rs:90-92`.
    pub stars_detected: usize,
    /// Frame aggregates, 0.0 when no fit stands in the frame statistic (PC-7): `analyzer.rs:146, :168-171`.
    pub median_fwhm: f64,
    pub median_eccentricity: f64,
    pub median_snr: f64,
    pub median_hfr: f64,
    /// background / noise, a ratio (PC-11); 0.0 without noise: `analyzer.rs:133`.
    pub frame_snr: f64,
    /// flux² / (noise² · background) over the in-statistic fits' median flux: `analyzer.rs:134-145`.
    pub snr_weight: f64,
    /// median amplitude / noise: `analyzer.rs:142`.
    pub psf_signal: f64,
    /// The sky median and the noise σ, non-finite written 0.0: `analyzer.rs:131-132`.
    pub background: f64,
    pub noise: f64,
    /// Rayleigh R² of the in-statistic fits' θ (`trail::trail_r_squared`), on the UNZEROED field FWHM: `analyzer.rs:147, :181`.
    pub trail_r_squared: f64,
    /// The fitted β under `PsfModel::Auto` when its source is `Fitted`, else `None` (PC-9): `analyzer.rs:148-151`.
    pub median_beta: Option<f64>,
    /// The refusal reasons' sum; the budget's not-selected / trimmed are not in it: `analyzer.rs:188`.
    pub stars_refused: usize,
    /// Trailed / defocused as the module judged, else no stars when no fit stands in the statistic (PC-13): `analyzer.rs:153-161`.
    pub verdict: Option<Verdict>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Verdict { Trailed, Defocused, NoStars }

impl FrameSummary {
    pub fn from_measurement(pm: &PlaneMeasurement, opts: &FrameOptions) -> FrameSummary;
}
```

The body of `from_measurement` is `analyzer.rs:130-161` verbatim, with `config.psf_model`'s `AnalysisPsfModel::Auto`
read as `opts.psf_model`'s `PsfModel::Auto` (the two enums map one to one, `analyzer.rs:40-43`), followed by the
struct literal built from the expressions of `analyzer.rs:167-188`. `finite_or_zero` and `point_like_seeds`
(`analyzer.rs:79-92`) move with it, private; `point_like_seeds` returns the `usize` and Athenaeum casts `as i64`, as
it does today inside the helper. `finite_median` is `crate::star_measure::frame::finite_median` (it is in `adapter`,
D-3; after phase 4 Athenaeum no longer calls it, and it stays listed: removing a doc-hidden item gains nothing).
`measure_version` is not a field: it is the constant `MEASURE_VERSION`.

```rust
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct FrameOptions {
    /// Matched-filter detection threshold, σ (Athenaeum validates 2–20, `seeds::MIN_THRESHOLD..=MAX_THRESHOLD`).
    pub detection_threshold: f64,
    /// `Moffat4` (β fixed at 4) or `Auto` (the frame's own β).
    pub psf_model: PsfModel,
    /// The measurement budget N: the brightest good stars kept over the field (`Mode::Budget`).
    pub max_stars: usize,
}

impl Default for FrameOptions {
    fn default() -> Self {
        FrameOptions { detection_threshold: DEFAULT_THRESHOLD, psf_model: PsfModel::Moffat4, max_stars: 500 }
    }
}

pub fn measure_raw_frame(frame: &RawFrame, opts: &FrameOptions, pool: Option<&Arc<rayon::ThreadPool>>) -> PlaneMeasurement;
pub fn measure_raw_frame_with(frame: &RawFrame, opts: &FrameOptions, pool: Option<&Arc<rayon::ThreadPool>>, diagnostics: bool) -> PlaneMeasurement;
```

- The defaults are Athenaeum's: `DEFAULT_THRESHOLD` is 6.0 (`src/star_measure/seeds.rs:66`; Athenaeum
  `config.rs:6-8`), 500 (`config.rs:9-11`), Moffat 4 (`config.rs:57`). Note that `PsfModel::default()` is `Auto`
  (`src/star_measure/fit.rs:89-90`), so `FrameOptions::default()` must name `Moffat4`.
- `measure_raw_frame_with`'s body is `analyzer.rs:40-76` verbatim, with `config.detection_threshold` →
  `opts.detection_threshold`, `config.max_stars as usize` → `opts.max_stars`, the `psf_model` match replaced by
  `opts.psf_model`, and `cfa::measure_green_sites` → `super::green::measure_green_sites`. `measure_raw_frame` is
  `analyzer.rs:24-30`.
- `psf_model` is the module's `PsfModel`, not a new enum: the CLI (phase 5, D-7) parses `moffat4|auto` into it, and
  `Fixed(β)` measures as given with `median_beta` `None` (it is not `Auto`).
- **Stays in Athenaeum** (`analysis/analyzer.rs`, rewritten in T6): `frame_options(&AnalysisConfig) -> FrameOptions`
  (the mapping of `analyzer.rs:40-51`), `measure_raw_frame` / `measure_raw_frame_with` taking `&AnalysisConfig` (thin
  wrappers, so `api/analysis.rs`, `star_probe` and the tests keep their calls), `map_measurement` (the `StarMetric` rows
  of `analyzer.rs:107-128` unchanged, the `FrameAnalysis` literal from the summary), `impl From<Verdict> for
  FrameVerdict`, `analyze_raw_frame`, `analyze_frame`; and `config.rs` whole.

### 1.4 Every caller

`grep -rn 'analysis::\(input\|cfa\|trail\|analyzer\)\|analyzer::'` over `ATH/crates` (outside `analysis/`), plus the
symbols (`read_raw_frame`, `measure_green_sites`, `trail_r_squared`, `measure_raw_frame`, `map_measurement`,
`RawPlanes`, `native_seed`, `SUPERPIXEL_CENTRE_OFFSET`, `fold_bayer_offsets`):

| caller | today | after T6 |
| ---- | ---- | ---- |
| `ATH/…/api/analysis.rs:31, :322` (`analyze_frame_set`, the batch; Tauri `commands/analysis.rs:53-63`, web `routes/analysis.rs:119-128`) | `analyzer::analyze_frame` | unchanged |
| `ATH/…/api/analysis.rs:569` (`get_frame_star_metrics`: the Analysis tab's frame and **Blink**'s overlay, `src/hooks/useStarMetricsCache.ts:40` ← `BlinkViewer.tsx`; Tauri `commands/analysis.rs:88-95`, web `routes/analysis.rs:150-157`) | `analyzer::analyze_frame` | unchanged |
| the plate-solve input gate (`ATH/…/plate_solve/service.rs:653-654`) | reads the stored row's `measure_version` against `star_measure::MEASURE_VERSION`; calls no analysis code | unchanged |
| the scanner (`ATH/…/scanner/`) | no call site (it reads headers through `fits_parser`) | — |
| `ATH/…/examples/star_probe.rs:156-158` (the analysis-mode arm, `:1166-1300`) | `analysis::analyzer::{map_measurement, measure_raw_frame_with}`, `analysis::config::…`, `analysis::input::{read_raw_frame, RawPlanes}` | line 158 → `athenaeum_core::star_measure::raw::{read_raw_frame, RawPlanes}`; 156-157 unchanged |
| `ATH/…/ts_export.rs:394-395`, `api/settings.rs:7`, Tauri / web `analysis.rs:10`, web `analysis.rs:175` | `analysis::config::…` | unchanged |
| tests: `analysis/analyzer.rs` (10), `analysis/config.rs` (5), `api/analysis.rs` (3), `db/analysis.rs` (8) | Athenaeum types | unchanged (bodies) |

No command or route changes signature, so the two-backend rule is met with no route edit: both hosts call
`api::analyze_frame_set` (`api/analysis.rs:153`) and `api::get_frame_star_metrics` (`api/analysis.rs:499`), whose
signatures do not change. No other example (`grep` over `ATH/crates/athenaeum-core/examples/`) and no integration test
(`tests/`) uses the moving modules.

### 1.5 Tests

| file | tests | where | why |
| ---- | ---- | ---- | ---- |
| `input.rs` | 11: `bitpix16_bzero0_reads_signed`, `bitpix16_bzero32768_reads_unsigned`, `bitpix16_bzero0_holding_offset_data_reads_with_the_32768_offset`, `a_unit_float_fits_is_scaled_to_adu`, `an_adu_float_fits_keeps_its_values`, `datamax_sets_clip_and_egain_sets_gain`, `an_8bit_fits_reads_in_its_own_adu`, `bitpix32_holding_16bit_data_takes_full_scale_65535`, `datamax_scales_only_a_unit_range_plane`, `a_bayer_frame_reads_as_cfa_with_its_layout`, `a_three_plane_frame_reads_as_rgb`; and `dsc5767_reads_with_the_32768_offset` (`#[ignore]`, `STAR_TEST_DIR`) | → `raw.rs`, bodies unchanged | hand-written FITS (`write_fits`, `input.rs:339-357`) and the local `warn_log` |
| `input.rs` | `xisf_uint16_and_float32_read_at_full_scale_65535`, `an_xisf_mosaic_folds_its_bayer_offsets` | stay in Athenaeum (`analysis/raw_seam.rs`), bodies unchanged | they write XISF with Athenaeum's `fits_writer::{write_xisf_f32, Card, CardValue}` (`input.rs:288`); they hold Athenaeum's writer to rustafits' reader |
| `input.rs` | `bayer_offsets_fold_into_the_pattern` | → `raw.rs` as `a_fits_mosaic_with_an_odd_xbayroff_reads_folded`: its file half (`input.rs:638-656`) unchanged; its table half tested `fold_bayer_offsets`, which goes (decision P4-2) and is covered by `formats::fits::tests::bayer_offsets_fold_into_the_pattern` (`src/formats/fits.rs:594-627`, all four patterns, `(3, -1)` included) and T1's XISF twin | |
| `cfa.rs` | all 10 | → `green.rs`, bodies unchanged | `crate::star_measure`, `test_support::{add_noise, moffat_field, MoffatStar, gaussian_field}`, `geometry::ransac::SplitMix64`: all in rustafits (`src/test_support.rs:11, :140, :158, :173`; `src/geometry/mod.rs:7-9`) |
| `trail.rs` | all 5 | → `trail.rs`, bodies unchanged | pure |
| `analyzer.rs` | all 10 | stay, bodies unchanged | they read `AnalysisConfig`, `FrameAnalysis`, `StarMetric`, `FrameVerdict` (the row mapping) |
| `config.rs` | all 5 | stay | the DB, the sha256 `config_hash` pin (`config.rs:197-199`) |

New tests (the plan's own, not moved): T1 3, T2 5, T4 9 (eight adapted from `analyzer.rs`'s, with `FrameSummary`
for `FrameAnalysis` and the row-only asserts left out; one for the defaults), T6 1 (the defaults' parity in Athenaeum).

**Fixtures.** rustafits' `test_support.rs` already holds every generator the moved and adapted tests call
(`gaussian_field`, `add_noise`, `MoffatStar`, `moffat_field`, `trailed_field`, `donut_field`: `src/test_support.rs:11,
:140, :158, :173, :482, :489`) and the RNG (`src/geometry/ransac.rs`, `#[cfg(test)]`). `input.rs`'s tests need
`tempfile`, which rustafits does not have (`src/formats/xisf.rs:866`): T2 adds `tempfile = "3"` as a
dev-dependency (Athenaeum's `Cargo.lock` already has 3.27.0, `Cargo.lock:7256`). `dsc5767_reads_with_the_32768_offset`
reads `singles/_DSC5767.fit` from the star-test corpus (`~/Pictures/star test`), outside both repos. No Athenaeum
`test_support` fixture becomes unused: the staying analyzer tests still call `trailed_field` and `donut_field`, and
`gaussian_field` / `add_noise` / `moffat_field` / `MoffatStar` have stacking callers.

**Pins.** None in the moved code. The only bit-level asserts are self-relative (`cfa.rs:286-287`: two measurements of
one run compared `to_bits`); `input.rs`'s exact values are exact arithmetic (e.g. `0.8f32 * 65535.0`); the analyzer
tests use tolerances. The one pin of the Analysis layer, `config_hash`'s sha256 (`config.rs:197-199`), stays in
Athenaeum and pins `MEASURE_VERSION` 8.

### 1.6 serde and TypeScript

No moved type crosses Athenaeum's wire or a TS export. `RawFrame` / `RawPlanes` have no serde derive
(`input.rs:29-51`); `PlaneMeasurement` is not stored by the Analysis path; the stored and exported types are
Athenaeum's own and stay: `FrameAnalysis`, `StarMetric`, `FrameVerdict` (`ATH/…/models.rs:1066-1129`) and
`AnalysisConfig`, `AnalysisPsfModel` (`config.rs:22-51`, exported by `ts_export.rs:394-395`). So:
- D-4: the phase-4 types (`RawFrame`, `RawPlanes`, `FrameOptions`, `FrameSummary`, `Verdict`) get no serde derive;
  nothing stores them, and `BayerPattern` has none (`src/types.rs:1-9`). Phase 5 adds one if the CLI wants it;
- D-5: no TS mirror is needed (phase 3 needed `PsfModelTs` / `ScaleEstimatorTs` because those types sat in exported
  configs; nothing here does);
- `ts_contract` must pass with `src/types` byte-identical (T6 acceptance).

### 1.7 Logging

Moved events (fields unchanged):

| event | level, fields | old target | new target |
| ---- | ---- | ---- | ---- |
| "raw frame read failed" (`input.rs:149`) | error; `path`, `error` | `athenaeum_core::analysis::input` | `astroimage::star_measure::raw` |
| "float file above 65535 ADU, full scale follows its maximum" (`input.rs:264`) | warn; `path`, `finite_max` | same | same |
| "BITPIX 16 / BZERO 0 file read with the 32768 offset" (`input.rs:278`) | warn; `path`, `plane_max` | same | same |
| "BITPIX 16 / BZERO 0 full scale chosen" (`input.rs:280`) | debug; `path`, `plane_max`, `full_scale` | same | same |
| "cfa green sites: no noise estimate on the green phase" (`cfa.rs:70-74`) | error; `width`, `height` | `athenaeum_core::analysis::cfa` | `astroimage::star_measure::green` |

`trail.rs` and the frame arithmetic log nothing. What references them:
- Athenaeum's filter has no `analysis` module key and no `astroimage` directive (`ATH/…/logging/config.rs:37-66`); the
  base level is `info` (`:81-87`), so the error / warn events stay visible and the debug one stays hidden, as today;
- no test asserts these targets: `input.rs`'s `warn_log` filters by level (`input.rs:308`) and moves with the tests;
- no log-mcp query and no other code names them (`grep` of the messages over `ATH`);
- docs: `ATH/docs/logging/README.md:179` and the logging spec's dictionary line
  (`ATH/docs/superpowers/specs/2026-07-03-logging-overhaul-design.md:64`) name `analysis::input` as the emitter: T7
  updates both. No field name is new, so the dictionary itself is unchanged.

### 1.8 Gating and build cost

- `star_measure` is always compiled (`src/lib.rs:16`), and so are the new submodules. `raw` needs `ImageConverter`
  (`src/lib.rs:30`) and `formats` (private but always compiled, `src/lib.rs:8-11`); `green` needs `processing::cfa`
  (`src/lib.rs:15`); both need only `anyhow`, `rayon`, `tracing`, already regular dependencies (`Cargo.toml:46-72`).
  No feature is needed, and `cargo build --no-default-features` compiles them.
- `cli` (the default, `Cargo.toml:32-35`) pulls `clap`, `clap_complete`, `indicatif` for the binary only; the library
  does not depend on it.
- R-7: phase 4 adds about 650 non-test lines (`input.rs` 284, `cfa.rs` 124, `trail.rs` 25, ~130 of `analyzer.rs`,
  plus the summary type) to a library that already compiles ~12k the CLI does not call. The CLI calls none of it until
  phase 5, so the release binary (LTO) should not grow measurably (phase 2: +752 bytes for the whole module). G-1g
  records the size.

## 2. Decisions this plan takes (review them)

- **P4-1 Header facts** come from `ImageConverter::read_header` with Athenaeum's numeric rules copied (`fits_f64` =
  `parse_fits_f64`, `fits_i32` = `get_i32`'s rule, `num` for XISF); the XISF arm keeps BITPIX/BZERO/BSCALE fixed; the
  format is `header.format`. (Alternative: the reader's own parsed `FitsHeader` values. Closer to the decode, but a
  number could move on a `D`-exponent card; rejected under "no number moves".)
- **P4-2 `raw` folds nothing.** `fold_bayer_offsets`, `offset` and the two offset fields of `HeaderFacts` go; the
  pattern passes through. (Alternative: keep the XISF arm at 0 / 0 as the FITS arm is, `input.rs:138-141`, which keeps
  `bayer_offsets_fold_into_the_pattern` verbatim but leaves a no-op fold and an unused `offset` in a semver-covered
  module.)
- **P4-3 The XISF reader's offset rule is Athenaeum's** (section 1.2), not the FITS reader's integer prefix.
- **P4-4 `#[non_exhaustive]` on `FrameOptions` and `FrameSummary`**: fields can be added in 1.x. `RawFrame`, `RawPlanes`
  and `Verdict` stay exhaustive: Athenaeum's kept tests build `RawFrame` literals (`analyzer.rs:265-272`), `star_probe`
  matches `RawPlanes` exhaustively (`star_probe.rs:1215-1219`), and a new verdict should break Athenaeum's mapping at
  compile time. `analysis_mode` and `raw` are the semver-covered entry (D-3); `green` and `trail` are public at the
  module's own status.
- **P4-5 `tempfile` dev-dependency** rather than a temp-dir shim, so the moved bodies stay unchanged.
- **P4-6 Athenaeum keeps thin wrappers** (`measure_raw_frame(&AnalysisConfig)`, …), so `api/analysis.rs`, `star_probe`
  and the ten analyzer tests keep their calls; it re-exports `RawFrame` / `RawPlanes` from `analysis::analyzer` (its
  public functions take them, and the kept tests reach them through `super::*`).
- **P4-7 No copy script.** Phase 2's `port_star_measure.py` regenerated 12k lines; here three small files are copied by
  hand and reviewed with `git diff --no-index --word-diff` against the originals (each task lists the allowed
  differences).

## 3. Tasks

### T1 (rustafits) — 4.1: the XISF reader folds `XBAYROFF` / `YBAYROFF`

**Files:** `src/formats/mod.rs`, `src/formats/fits.rs`, `src/formats/xisf.rs`, `CHANGELOG.md`.

**Approach.**
- Move `shift_bayer_pattern` (`src/formats/fits.rs:179-206`, with its doc) to `src/formats/mod.rs` as `pub(crate)`;
  `layout` calls `super::shift_bayer_pattern`. Nothing else in `fits.rs` changes.
- In `parse_xisf_xml`'s event loop (`xisf.rs:261-286`) add one arm for `Event::Empty` / `Event::Start` named
  `FITSKeyword`: keep the first `value` of a `name` equal (ignoring case) to `XBAYROFF`, and of `YBAYROFF`. After the
  pattern scan (`xisf.rs:312-321`): `info.bayer_pattern = shift_bayer_pattern(info.bayer_pattern, odd(x), odd(y))`
  with `fn bayer_offset(v: Option<&str>) -> i64` = Athenaeum's `offset(num(..))` (`input.rs:76-79, :111-114`: trim,
  trim `'`, trim, `parse::<f64>()`, finite, `round() as i64`, else 0) and odd = `rem_euclid(2) == 1`. Comment the
  rule's origin (the rule the Analysis tab applied to XISF before the port, so its numbers do not move). Attribute
  values through `normalized_value(XmlVersion::Implicit1_0)`, as `parse_image_attrs` does (`xisf.rs:189-191`).
- The pattern scan itself (`xml.contains`, `xisf.rs:312-321`) does not change.
- CHANGELOG `[Unreleased]` → `### Fixed`: the entry of section 1.2.

**Tests** (in `xisf.rs`'s test module, on `write_xisf_header_and_attachments`, `xisf.rs:878-903`; a UInt16 mosaic,
`<FITSKeyword>`s inside the `<Image>` element):
- `xisf_bayer_offsets_fold_into_the_pattern`: the table of `src/formats/fits.rs:594-627` (four patterns × offsets None,
  (0,0), (1,0), (0,1), (1,1), (2,4), (3,-1)), asserted on `read_xisf_image` AND `read_xisf_header_info`;
- `xisf_bayer_offset_values_read_as_athenaeum_reads_them`: RGGB with `XBAYROFF` = `1`, `1.0`, `'1'`, `-1`, `3`, `0.6`
  → GRBG; `2`, `x`, absent → RGGB;
- `xisf_red_star_stays_red_for_every_bayer_offset`: `src/formats/fits.rs:632-667` on an XISF UInt16 mosaic, super-pixel
  and VNG.

**Acceptance.**
```bash
cargo test -p rustafits --lib formats::          # all pass, 3 more than at 8ecca88
cargo test -p rustafits --test analysis --test cli   # test.xisf has XBAYROFF / YBAYROFF 0: unchanged
cargo check -p athenaeum-core --all-targets      # Athenaeum compiles against T1 (its XISF tripwire test fails: section 1.2)
```
Commit: `fix(formats): XISF folds XBAYROFF / YBAYROFF into the Bayer pattern (port plan 4.1)`.

### T2 (rustafits) — 4.2: `star_measure::raw`

**Files:** `src/star_measure/raw.rs` (new), `src/star_measure/mod.rs` (`pub mod raw;`), `Cargo.toml`
(`[dev-dependencies] tempfile = "3"`), `src/formats/xisf.rs:866` (the comment "no `tempfile` dependency in this
crate" becomes false: say the XISF tests keep their own scratch paths), `CHANGELOG.md`.

**Approach.** Copy `ATH/crates/athenaeum-core/src/analysis/input.rs` at e6344ddd. The only changes:
- module doc: one provenance paragraph (copied from Athenaeum's `analysis/input.rs` at e6344ddd, port plan phase 4; the
  header facts come from `read_header`; the readers fold the Bayer offsets; the rulings cited, PC-… / R-PC-…, are
  Athenaeum's); the `RawPlanes` doc sentence on the fold (`input.rs:26-27`) becomes "both rustafits readers fold them
  (R-PC-4)"; intra-doc paths to this crate's (`star_measure::Domain` → `super::Domain`);
- imports: `use crate::star_measure::Domain;` and `use crate::{BayerPattern, FileFormat, ImageConverter, ImageHeader,
  PixelData};` replace `astroimage::…`, `crate::fits_parser…`, `crate::fits_writer::OutputFormat`,
  `crate::models::FileFormat` (`input.rs:14-19`);
- `header_facts(path, xisf)` (`input.rs:116-144`) → `fn header_facts(header: &ImageHeader) -> HeaderFacts`: the XISF
  arm (`header.format == FileFormat::Xisf`) keeps `bitpix: None, bzero: 0.0, bscale: 1.0` and reads `DATAMAX` /
  `EGAIN` with `num`, now over `header.get(key)`; the FITS arm reads `BITPIX` with `fits_i32`, `BZERO` / `BSCALE`
  (defaults 0.0 / 1.0), `DATAMAX`, `EGAIN` with `fits_f64`. `fits_f64` is `parse_fits_f64`
  (`ATH/…/fits_header_reader.rs:248-255`) and `fits_i32` is `get_i32`'s closure (`:151-158`), both verbatim with a
  comment naming their source;
- `read_inner` (`input.rs:154-156`): `let header = ImageConverter::read_header(path).map_err(|e| anyhow::anyhow!("{}:
  {e:#}", path.display()))?;`, `let xisf = header.format == FileFormat::Xisf;`, `let facts = header_facts(&header);`;
- the fold goes (P4-2): `offset`, `fold_bayer_offsets`, `HeaderFacts::{xbayroff, ybayroff}` are deleted and
  `input.rs:238-241` reads `pattern => RawPlanes::Cfa { mosaic: …, pattern }`;
- clippy: a lint the copied code trips is added to an `#![allow]` at the file top, never fixed (G-3).
Everything else (`RawPlanes`, `RawFrame`, `drop_unmeasured_planes`, `read_raw_frame`, the domain rules, the two
full-scale helpers and their logs) is verbatim.

**Tests.**
- moved, bodies unchanged: the 11 + 1 ignored of section 1.5, with `use crate::fits_writer::…` dropped (only the
  staying XISF tests used it) and `astroimage::` → `crate::`;
- new: `a_fits_mosaic_with_an_odd_xbayroff_reads_folded` (`input.rs:638-656` unchanged as its body);
  `an_xisf_mosaic_reads_its_bayer_offsets_folded_once` (a hand-written monolithic UInt16 XISF built as the UInt16 half
  of `input.rs:506-522` does, with `<FITSKeyword name="BAYERPAT" value="'RGGB'"/>` and `XBAYROFF` 1 → GRBG, and with
  `YBAYROFF` 1 as well → BGGR); `xisf_planes_read_at_full_scale_65535` (that UInt16 file, asserts of `input.rs:524-525`;
  plus a hand-written Float32 XISF with `DATAMAX` 0.9 and `EGAIN` 1.5: sample 1.0 → 65535, full scale 65535, gain 1.5,
  clip 0.9 × 65535 within 1e-6, the asserts of `input.rs:501-504`); `an_xisfs_fits_layout_keywords_are_ignored` (the
  Float32 file with `BITPIX` 16, `BZERO` 32768, `BSCALE` 1 keywords, as the corpus' `Light_Pane` XISF carry, reads
  exactly like the one without them); `fits_numbers_read_as_athenaeum_reads_them` (`fits_f64`: `3.2768D4` → 32768,
  `1.5d0` → 1.5, empty → None; `fits_i32`: `16`, `16.0` → 16, `-32` → -32). The hand-written XISF need no namespace
  attribute.

**Acceptance.**
```bash
cargo test -p rustafits --lib star_measure::raw::     # 16 passed, 1 ignored
STAR_TEST_DIR="$HOME/Pictures/star test" cargo test -p rustafits --lib star_measure::raw:: -- --ignored   # dsc5767 passes
git diff --no-index --word-diff crates/athenaeum-core/src/analysis/input.rs rustafits/src/star_measure/raw.rs
#   only the changes listed above
cargo build -p rustafits --no-default-features
cargo check -p athenaeum-core --all-targets
```
`Cargo.toml`'s dev-dependency changes Athenaeum's `Cargo.lock` (rustafits is a workspace member): leave it
uncommitted; T6 commits it. Commit: `feat(star_measure): raw, the Analysis tab's raw frame input, from Athenaeum
e6344ddd (port plan 4.2)`.

### T3 (rustafits) — 4.3: `star_measure::green` and `star_measure::trail`

**Files:** `src/star_measure/green.rs`, `src/star_measure/trail.rs` (new), `src/star_measure/mod.rs`, `CHANGELOG.md`.

**Approach.** Copy `cfa.rs` → `green.rs` and `trail.rs` → `trail.rs` at e6344ddd. Only: `astroimage::processing::cfa`
→ `crate::processing::cfa` (`cfa.rs:15`), `astroimage::BayerPattern` → `crate::BayerPattern` (`cfa.rs:16`), the
module doc's `input::RawPlanes` → `raw::RawPlanes` (`cfa.rs:10`), one provenance line each. `crate::star_measure::…`
paths are already right. `pub mod green; pub mod trail;` in `mod.rs` (beside `pub mod frame;`, alphabetical).

**Tests.** `cfa.rs`'s 10 and `trail.rs`'s 5, bodies unchanged (paths as above).

**Acceptance.**
```bash
cargo test -p rustafits --lib -- star_measure::green:: star_measure::trail::   # 15 passed
git diff --no-index --word-diff crates/athenaeum-core/src/analysis/cfa.rs rustafits/src/star_measure/green.rs     # paths and provenance only
git diff --no-index --word-diff crates/athenaeum-core/src/analysis/trail.rs rustafits/src/star_measure/trail.rs   # provenance only
cargo check -p athenaeum-core --all-targets
```
Commit: `feat(star_measure): green and trail from Athenaeum e6344ddd (port plan 4.3)`.

### T4 (rustafits) — 4.4: `star_measure::analysis_mode`

**Files:** `src/star_measure/analysis_mode.rs` (new), `src/star_measure/mod.rs`, `CHANGELOG.md`.

**Approach.** Section 1.3 exactly: `FrameOptions` + `Default`, `measure_raw_frame`, `measure_raw_frame_with`,
`Verdict`, `FrameSummary::from_measurement`, the private `finite_or_zero` and `point_like_seeds`. Module doc: what it
is (Athenaeum's Analysis tab measurement, Plan C C5), that it and `raw` are the stable entry (D-3), that the row
mapping stays in Athenaeum, provenance (`analyzer.rs` at e6344ddd). `pub mod analysis_mode;` in `mod.rs`; the module
doc's provenance paragraph (`src/star_measure/mod.rs:10-13`) gains one sentence on phase 4's four submodules.

**Tests** (new, adapted from `analyzer.rs`: `fn summary(f: &RawFrame, o: &FrameOptions) -> FrameSummary` replaces
`analyze_raw_frame(&f, &cfg(), None, "h")`, `FrameOptions::default()` replaces `AnalysisConfig::default()`,
`PsfModel::Auto` replaces `AnalysisPsfModel::Auto`, `Verdict` replaces `FrameVerdict`, the fixtures and assertions are
the originals' minus those on row-only fields):
`frame_options_default_is_the_analysis_tab_default` (6.0 = `DEFAULT_THRESHOLD`, `Moffat4`, 500);
`a_synthetic_mono_frame_summarizes_every_field` (`analyzer.rs:297-323`); `stars_detected_leaves_knots_out_in_every_layout`
(`:328-382`); `a_frame_without_fits_summarizes_to_zero_medians` (`:384-410`);
`a_trailed_frame_is_summarized_with_no_measured_star` (`:412-435`); `stars_refused_is_the_refusal_reasons_alone`
(`:448-481`); `the_summary_records_the_frame_verdict` (`:496-515`); `moffat4_summarizes_no_median_beta_and_auto_the_fitted_one`
(`:517-529`); `an_osc_mosaic_is_measured_on_green_and_a_three_plane_file_on_g` (`:540-632`, positions from
`pm.stars`). `a_non_finite_star_snr_or_hfr_maps_to_none` is row-only (PC-8) and has no twin.

**Acceptance.**
```bash
cargo test -p rustafits --lib -- star_measure::analysis_mode::    # 9 passed
cargo test -p rustafits --features serde                          # aarch64: 556 passed, 0 failed, 3 ignored
cargo test -p rustafits --features debug-pipeline
cargo build -p rustafits --no-default-features
cargo clippy -p rustafits --all-targets --message-format=short 2>&1 | grep -E "src/(star_measure|formats)/"   # empty
cargo check -p athenaeum-core --all-targets
```
556 = 513 at 8ecca88 (port plan l.337: the code is 0025614's, `git diff --stat 0025614 8ecca88` touches no source)
+ 26 moved + 17 new; 3 ignored = 2 + `dsc5767`. Commit: `feat(star_measure): analysis_mode, FrameOptions /
measure_raw_frame / FrameSummary (port plan 4.4)`.

### T5 (rustafits) — docs

**Files:** `CLAUDE.md`, `src/formats/xisf.rs` (doc comment only), `scripts/star_measure_port/mac_check.sh`,
`docs/superpowers/plans/2026-11-03-star-measure-port.md`, `CHANGELOG.md`.

- `CLAUDE.md`: pipeline step 1 (FITS and XISF fold the offsets); the PRODUCTION CONTRACT paragraph names
  `star_measure::raw::read_raw_frame` (this crate since phase 4) instead of Athenaeum's `analysis::input` /
  `analysis::analyzer::analyze_frame`; the module-map `star_measure/` line lists `raw`, `green`, `trail`,
  `analysis_mode` (the entry); Downstream Users: Athenaeum builds with `features = ["serde"]` (stale since phase 3) and
  its Analysis tab, Blink and plate-solve gate measure through `analysis_mode`. The rest of that section is phase 5.3's.
- `xisf.rs:680-681`: the same contract wording (`star_measure::raw`).
- `mac_check.sh`: `BRANCH=${MAC_CHECK_BRANCH:-star-measure-port-phase-3}` (and its usage text), so phase 4 reuses it.
- the port plan: phase 4's items point to this plan; results are added by the controller after the gates.
- `CHANGELOG.md` `[Unreleased]` `### Added` (T2–T4 wrote their lines; T5 checks it reads as one entry):

  > - `star_measure::raw`: `read_raw_frame` reads a FITS / XISF file as a `RawFrame` (mono, a Bayer mosaic or RGB
  >   planes, in the file's ADU) with the facts the measurement needs: full scale, the clip from `DATAMAX`, the gain
  >   from `EGAIN`.
  > - `star_measure::green` (a mosaic measured on its green sites) and `star_measure::trail` (trail R²).
  > - `star_measure::analysis_mode`: `FrameOptions` (6.0 σ, Moffat 4, 500 stars), `measure_raw_frame` and
  >   `FrameSummary`: the measurement and the frame figures Athenaeum's Analysis tab shows, copied from Athenaeum,
  >   which now calls them. The stable entry of `star_measure`.

**Acceptance.**
```bash
grep -n "analysis::input\|analysis::analyzer" rustafits/CLAUDE.md rustafits/src/formats/xisf.rs   # empty
bash -n rustafits/scripts/star_measure_port/mac_check.sh
git -C rustafits diff --stat HEAD~1      # only the files above
```
Commit: `docs: star_measure's Analysis-mode entry (port plan phase 4)`.

### T6 (Athenaeum) — 4.5: the switch, in ONE commit

**Files:** `rustafits` (submodule → the T5 commit), `Cargo.lock`, `crates/athenaeum-core/src/analysis/{mod,analyzer}.rs`,
`crates/athenaeum-core/src/analysis/raw_seam.rs` (new), `crates/athenaeum-core/examples/star_probe.rs`; deleted
`crates/athenaeum-core/src/analysis/{input,cfa,trail}.rs`.

**Approach.**
- `analysis/mod.rs`: `pub mod config; pub mod analyzer; #[cfg(test)] mod raw_seam;`.
- `analysis/analyzer.rs`: imports `super::config::{AnalysisConfig, AnalysisPsfModel}`,
  `crate::models::{FrameAnalysis, FrameVerdict, StarMetric}`, `crate::star_measure::analysis_mode::{self, FrameOptions,
  FrameSummary, Verdict}`, `crate::star_measure::raw::read_raw_frame`, `pub use crate::star_measure::raw::{RawFrame,
  RawPlanes};`, `crate::star_measure::{PlaneMeasurement, PsfModel, MEASURE_VERSION}`.
  - `pub fn frame_options(config: &AnalysisConfig) -> FrameOptions`: `FrameOptions::default()`, then
    `detection_threshold`, `psf_model` (the match of `analyzer.rs:40-43`) and `max_stars` (`config.max_stars as usize`)
    assigned (non-exhaustive struct);
  - `measure_raw_frame` / `measure_raw_frame_with` keep their signatures and call `analysis_mode::…` with
    `&frame_options(config)`;
  - `impl From<Verdict> for FrameVerdict` (one arm per variant);
  - `map_measurement` keeps its signature: the `StarMetric` rows verbatim (`analyzer.rs:107-128`); `let s =
    FrameSummary::from_measurement(pm, &frame_options(config));`; the `FrameAnalysis` literal of `analyzer.rs:163-190`
    with each computed field read from `s` (`stars_detected: s.stars_detected as i64`, `stars_refused:
    Some(s.stars_refused as i64)`, `frame_verdict: s.verdict.map(FrameVerdict::from)`, …), `possibly_trailed:
    pm.trailed`, `measure_version: Some(MEASURE_VERSION)`, the row facts as before;
  - `analyze_raw_frame`, `analyze_frame` and the test module unchanged; the module doc says where the arithmetic lives.
- `analysis/raw_seam.rs` (test-only): `xisf_uint16_and_float32_read_at_full_scale_65535` and
  `an_xisf_mosaic_folds_its_bayer_offsets` with bodies unchanged, the two helpers they call (`mono`, `raw`,
  `input.rs:365-380`) and their imports (`crate::fits_writer::{write_xisf_f32, Card, CardValue}`,
  `crate::star_measure::raw::{read_raw_frame, RawFrame, RawPlanes}`, `crate::star_measure::Domain`,
  `astroimage::BayerPattern`). The fold test's doc comment (`input.rs:610-611`) is rewritten: the rustafits XISF reader
  folds since phase 4, and this test holds the one fold across Athenaeum's writer, the reader and `raw` (zero or two
  folds read RGGB). New: `the_default_analysis_config_measures_with_the_default_frame_options`
  (`analyzer::frame_options(&AnalysisConfig::default()) == FrameOptions::default()`).
- `examples/star_probe.rs:158` → `use athenaeum_core::star_measure::raw::{read_raw_frame, RawPlanes};`.
- The submodule: `git -C rustafits checkout <T5 sha>`, `git add rustafits Cargo.lock`. The bump and the deletion of
  `input.rs` are the same commit (section 1.2); never commit the bump alone.

**Tests.** `analyzer.rs`'s 10 and `config.rs`'s 5 unchanged; `raw_seam.rs` 2 moved + 1 new.

**Acceptance.**
```bash
cargo test -p athenaeum-core --lib analysis       # every test passes; analysis/ itself: 18 passed (was 44 + 1 ignored)
cargo test -p athenaeum-core --test ts_contract && git status --porcelain src/types   # pass, empty
cargo check -p athenaeum-core --no-default-features
cargo check -p athenaeum-web
grep -rn 'analysis::\(input\|cfa\|trail\)\|fold_bayer_offsets' crates     # empty
cargo test -p athenaeum-core --all-targets        # aarch64: 4549 passed, 0 failed, 15 ignored, no pin edited
```
4549 = 4575 at e6344ddd (port plan l.338, measured at ca21e19c; `git diff --stat ca21e19c e6344ddd` touches no source,
and its submodule move 0025614 → 8ecca88 touches none either) − 12
(`input.rs`, moved) − 10 (`cfa.rs`) − 5 (`trail.rs`) + 1 (the parity test); 15 ignored = 16 − `dsc5767`.
Commit: `refactor(analysis): measure through rustafits' star_measure::analysis_mode (port plan 4.5)`, naming the
submodule sha.

### T7 (Athenaeum) — docs

**Files:** `docs/logging/README.md`, `docs/superpowers/specs/2026-07-03-logging-overhaul-design.md`,
`docs/stacking/README.md`, `docs/superpowers/open-items.md`.

- logging README l.179 and the spec's l.64: the emitter is rustafits' `star_measure::raw::read_raw_frame` (target
  `astroimage::star_measure::raw`), the CFA error is `star_measure::green`'s; fields unchanged;
- `docs/stacking/README.md`: under R-OSP-19, a phase-4 sub-entry (the four files, what stays, `FrameSummary`, the
  fold in the reader, the gates and their results); R-PC-4's evidence line (l.4414-4415) gains "phase 4: the XISF fold
  is the rustafits XISF reader's, `star_measure::raw` folds nothing";
- open-items: an owed smoke, "an XISF mosaic with an odd XBAYROFF / YBAYROFF (none in the star-test corpus): its
  preview shows the right colours, and its Analysis row equals the pre-phase-4 one" (deleted once it passes);
- `CLAUDE.md` needs no change: no rule, command surface or path it names changed.

**Acceptance:** `grep -rn "analysis::input\|analysis/input.rs" docs/logging docs/superpowers/specs/2026-07-03-logging-overhaul-design.md`
returns only historical mentions marked as such. Commit: `docs: the Analysis-mode entry lives in rustafits (R-OSP-19,
port phase 4)`.

### L1 (athenaeum-lab, the controller) — the harness takes the phase

`harness/acceptance/smport.sh` hard-codes its arms and its results folder. Three defaults, so phase 3's call is
unchanged: `ARMS=(${=SMPORT_ARMS:-measure analysis-mode})`,
`RES="$LAB/star-lab/results/star-measure-port${SMPORT_TAG:+/$SMPORT_TAG}"`, and the build target
`${SMPORT_TARGET:-$WORK/target}`. Usage text updated. Acceptance: `zsh -n harness/acceptance/smport.sh`; with no
variable set it prints the same usage. Committed in athenaeum-lab.

## 4. Gates (the controller, on the owner's Mac, after T7)

Work dir `~/.athenaeum-acc/star-measure-port-p4/` with a `REGISTRY.tsv` row (job, the controller, start, end = the day
the results are committed); it holds the worktrees, the harness target, the JSONL and the fold copies, and is removed
once the results are committed (`git worktree remove`, then the folder).

- **G-1a Suites, exact aarch64 pins** (`mac_check.sh`, T5): `MAC_CHECK_BRANCH=star-measure-port-phase-4
  rustafits/scripts/star_measure_port/mac_check.sh <athenaeum checkout> ~/.athenaeum-acc/star-measure-port-p4`. It
  fetches `main` and the branch from the checkout's `origin` (push the branch there first, or set `MAC_CHECK_REMOTE`;
  `main` must be e6344ddd) and makes the two worktrees `<work>/main` and `<work>/star-measure-port-phase-4`. Bars: the
  digest is byte-identical `main` against the branch (the module did not change); `cargo test -p rustafits --features
  serde` 556 / 0 / 3; `cargo test -p athenaeum-core --all-targets` 4549 / 0 / 15, no pin edited.
- **G-1b The corpus, analysis-mode arm** (the brief's item 9): `star_probe` (release) built at Athenaeum `main` e6344ddd
  and at the phase-4 branch, on every image of `~/Pictures/star test/{m42,m101,singles,trails}` (77 + 30 + 38 + 7 = 152),
  `measure <img> --paths analysis-mode --diagnostics`. From the athenaeum-lab checkout, after L1:
  ```bash
  cd /Volumes/BigMac/Users/astrobureau/Documents/Projects/athenaeum-lab
  SMPORT_ARMS=analysis-mode SMPORT_TAG=phase-4 zsh harness/acceptance/smport.sh \
      ~/.athenaeum-acc/star-measure-port-p4/main ~/.athenaeum-acc/star-measure-port-p4/star-measure-port-phase-4 \
      ~/.athenaeum-acc/star-measure-port-p4
  ```
  Bar: `SMPORT: PASS (152 image-arms: analysis-mode 152; identical 152, moved 0; … failed on both sides alike 0)`.
  Every JSONL line equal once `…Ms` and `analyzed_at` are dropped: star lines, refusal lines, the frame line with its
  `product` row (verdicts, `stars_refused`, `measure_version` included) and the options line (`clip`, `fullScale`,
  `gainEPerAdu`: the header facts). Phase 3's analysis-mode arm read all 152 with no failure. The Measure arm is not
  re-run: phase 4 changes nothing it reads (the FITS decode is untouched, and the XISF singles fail that arm on both
  sides by construction, R-OSP-19). Any moved row is a STOP.
- **Odd offsets in the corpus:** none (section 1.2). Had one existed, the analysis-mode arm would still not move: both
  sides fold once (Athenaeum's fold before, the reader's after); only the render changes on purpose. G-1c tests exactly
  that.
- **G-1c Positive control of the one fold.** Three copies of `singles/test.xisf` (the corpus copy is identical to
  `tests/test.xisf`; it carries `<FITSKeyword name="XBAYROFF" value="0" …/>` and `YBAYROFF` 0) with the header patched
  in place, one byte each, so the XML length and the attachment offsets do not move: X 1; Y 1; X 1 and Y 1. Put them in
  `<work>/fold-corpus/singles/` (empty `m42`, `m101`, `trails` beside it) and run
  `SMPORT_CORPUS=<work>/fold-corpus SMPORT_ARMS=analysis-mode SMPORT_TAG=phase-4-fold zsh harness/acceptance/smport.sh …`.
  Bars: PASS, 3 identical (one fold on each side). Against the unpatched file's after-JSONL (G-1b): the X-only and
  Y-only copies differ (their green sites move: GRBG / GBRG put green where RGGB has red and blue), and the X-and-Y copy
  is equal (BGGR and RGGB share their green sites, and the analysis-mode arm reads green only). Then `rustafits info
  --json` on the copies: RGGB before (8ecca88), GRBG / GBRG / BGGR after.
- **G-1d rustafits renders unchanged** (G-3): build `rustafits` (debug) in both worktrees, one after the other into the
  harness target, copy each binary out; run `rustafits <f> <out>.jpg`, `rustafits info --json` and `rustafits analyze
  --json` on `tests/mono.fits`, `tests/osc.fits`, `tests/test.xisf`. Bar: every output byte-identical (offsets 0
  everywhere). The fold copies' renders differ (expected).
- **G-1e `MEASURE_VERSION` 8:** `grep -n "pub const MEASURE_VERSION: u32 = 8;" rustafits/src/star_measure/mod.rs`, and
  `config_hash_folds_measure_version_and_ignores_batch_concurrency` passes (G-1a), so no stored analysis goes stale.
- **G-1f `ts_contract`, headless, web:** in T6's acceptance, re-run on the final commit.
- **G-1g Release binary size** (R-7): `cargo build --release -p rustafits` in both worktrees, one after the other; note
  the size change in the port plan.
- **G-3 rustafits:** `cargo test -p rustafits` (default), `--features serde`, `--features debug-pipeline`, `cargo
  build -p rustafits --no-default-features`, clippy with no new warning in `src/{star_measure,formats}/` (T4).

Results: `star-lab/results/star-measure-port/phase-4/` and `…/phase-4-fold/` in athenaeum-lab (committed), the counts in
the port plan's phase 4 and in R-OSP-19.

## 5. Risks

- **R-P4-1 Zero or two folds on XISF.** Mitigated by T1 before T2, `raw` folding nothing, T6 as one commit, three
  tripwire tests and G-1c.
- **R-P4-2 Header-fact semantics.** The copied numeric rules keep every standard and `D`-exponent value; the residual
  card-layout cases (section 1.1) appear in no corpus file. G-1b's options line proves the facts on the corpus.
- **R-P4-3 An XISF whose keyword pass errors** now fails `read_raw_frame` (section 1.1). Low; recorded.
- **R-P4-4 Intended render change downstream** (odd-offset XISF): Athenaeum's previews, solvemyastro's detection, the
  CLI. A Fixed entry; no corpus file; the open-items smoke.
- **R-P4-5 An XISF whose `<ColorFilterArray>` already includes its offset** would be folded twice by the reader's rule.
  This is the rule Athenaeum's Analysis already applies to XISF (R-PC-4), so the Analysis numbers keep their meaning;
  the corpus cannot tell (its one CFA XISF has offsets 0). Owner question Q-2.
- **R-P4-6 Semver** (R-3): `analysis_mode` / `raw` become the covered entry; `#[non_exhaustive]` on the two option /
  summary structs (P4-4).
- **R-P4-7 The squash merge** into `star-measure-port` (D-10) orphans the T1–T5 commits Athenaeum's T6 pins. Athenaeum
  re-pins to the squash commit, and later to `main`'s merge, each time with `git diff --stat <old> <new>` empty
  (identical trees need no gate re-run), as phase 3 did.
- **R-P4-8 The Athenaeum worktree double-folds between T1 and T6** (its submodule is the rustafits working copy):
  Athenaeum's tests are not run as a gate before T6.
- **R-P4-9 Disk and memory:** the gate builds a release `star_probe` twice and runs the suites; one compile at a time,
  the work dir registered and removed (CLAUDE.md disk rule).
- **R-6 Platform pins** (port plan): the counts above are aarch64's; elsewhere `assert_bits_pinned` tolerances apply,
  and Windows failures are checked against open-items' 2026-10-18 entry first.

## 6. Open questions (each with a recommended answer)

- **Q-1 When does phase 4 reach rustafits `main`?** Recommended: merge `star-measure-port` into `main` (merge commit,
  D-10) once phase 4's gates pass, with no tag (D-8: its additions ride with 1.8.0 in phase 5); Athenaeum's PR re-pins
  to that merge and merges on the owner's word. The alternative, holding Athenaeum's phase 4 until phase 5, keeps
  `analysis/` frozen (G-2) for longer.
- **Q-2 The XISF fold rule** (R-P4-5): fold the `<FITSKeyword>` offsets onto whatever pattern the reader found, as
  R-PC-4 does today? Recommended: yes, and if the owner has a real XISF mosaic with an odd offset, check its preview by
  eye (the open-items smoke).
- **Q-3 `#[non_exhaustive]` on `FrameOptions` and `FrameSummary`** (P4-4)? Recommended: yes; Athenaeum builds
  `FrameOptions` from `Default` plus field assignments.
- **Q-4 `raw` folds nothing (P4-2), so `bayer_offsets_fold_into_the_pattern` is split** rather than moved verbatim?
  Recommended: yes; a no-op fold in the semver-covered module is worse than one rewritten test whose table the reader's
  own tests already cover.
- **Q-5 `tempfile` as a rustafits dev-dependency** (P4-5)? Recommended: yes (3.27.0 is already in Athenaeum's lock).
- **Q-6 The lab harness change** (L1) in athenaeum-lab? Recommended: yes, defaults unchanged.
- **Q-7 The XISF offset rule** (P4-3): Athenaeum's `f64` + round, not the FITS reader's integer prefix? Recommended:
  yes; the two differ only on non-integer offsets, where the Analysis tab's numbers must not move.

## 7. Commits, PRs, merges

- rustafits: T1–T5 on `feat/star-measure-phase-4`; PR into `star-measure-port`, squash-merged (D-10). Pushing the
  branch and opening the PR are part of the approved phase; the merge into `main` waits for the owner (Q-1).
- Athenaeum: T6, T7 on `star-measure-port-phase-4`; PR into `main`; T6 pins the T5 commit, re-pinned as R-P4-7 says.
  The merge waits for the owner.
- athenaeum-lab: L1 and the results commit.
- No release (D-8). The G-2 freeze of `analysis/{input,cfa,trail,analyzer}.rs` holds until Athenaeum's phase-4 PR
  merges; a change that cannot wait goes into both copies.
