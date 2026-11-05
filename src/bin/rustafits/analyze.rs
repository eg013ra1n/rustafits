//! `rustafits analyze` — star metrics per frame, for culling subs.
//!
//! The numbers are Athenaeum's Analysis numbers: the frame is measured by `star_measure::analysis_mode` with its
//! defaults (6.0 σ matched-filter threshold, Moffat 4, a budget of 500 stars), so a file reads the same here as in
//! the Analysis tab; only the flags given override. A frame the module judges trailed, defocused or without stars has
//! no medians to show: the table prints `-` for them and `--sort` ranks the frame after every other.

use std::path::PathBuf;

use anyhow::Result;
use astroimage::star_measure::analysis_mode::{FrameOptions, Verdict};
use astroimage::star_measure::{PsfModel, DEFAULT_THRESHOLD, MEASURE_VERSION};
use clap::{Args, ValueEnum};

use crate::batch::{self, file_name};
use crate::convert::DEFAULT_MAX_JOBS_ANALYSIS;
use crate::measure::{
    budget_parser, measure_file, report_file_error, threshold_parser, warn_ignored, Measured, MEASURE_CAP_WHY,
};
use crate::report::{json_num, json_str, print_csv, print_table};

#[derive(Args, Debug)]
pub struct AnalyzeArgs {
    /// Input FITS/XISF files or directories
    #[arg(required = true, value_name = "INPUT")]
    paths: Vec<PathBuf>,

    /// Descend into subdirectories of directory inputs
    #[arg(short, long)]
    recursive: bool,

    /// Sort rows best-first by a metric (default: input order); frames with a verdict go last
    #[arg(short, long, value_enum, value_name = "METRIC")]
    sort: Option<SortKey>,

    /// Reverse the row order (worst-first with --sort)
    #[arg(long)]
    reverse: bool,

    /// JSON: one object per frame with every summary metric
    #[arg(long, conflicts_with = "csv")]
    json: bool,

    /// CSV at full precision
    #[arg(long)]
    csv: bool,

    /// Detection threshold in sigma above the noise, 2-20 [default: 6.0]
    #[arg(long, value_name = "SIGMA", value_parser = threshold_parser)]
    detection_threshold: Option<f64>,

    /// Star profile: Moffat with beta fixed at 4, or the frame's own beta
    #[arg(long, value_enum, value_name = "MODEL", default_value_t = PsfChoice::Moffat4)]
    psf_model: PsfChoice,

    /// Star budget: the brightest good stars measured per frame, 10-2000
    #[arg(long, value_name = "N", value_parser = budget_parser, default_value_t = 500)]
    max_stars: usize,

    /// Focal length in mm: adds FWHM/HFR in arcsec (with --pixel-size)
    #[arg(long, value_name = "MM", requires = "pixel_size")]
    focal_length: Option<f64>,

    /// Pixel size in µm (binned): adds FWHM/HFR in arcsec (with --focal-length)
    #[arg(long, value_name = "UM", requires = "focal_length")]
    pixel_size: Option<f64>,

    /// Frames analyzed concurrently [default: min(2, cores, files)]
    #[arg(short, long, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    jobs: Option<u32>,

    /// Threads for pixel work, shared by all jobs [default: all cores]
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    threads: Option<u32>,

    /// Errors only: no progress bar
    #[arg(long)]
    quiet: bool,

    // Flags of the pre-1.8.0 analyzer, kept so old scripts still run (gone in 2.0). A clap `alias` cannot tell that
    // it was used, and the warning needs that, so `--sigma` is an argument of its own.
    /// Old name of --detection-threshold: another detector since 1.8.0
    #[arg(long, value_name = "SIGMA", value_parser = threshold_parser, conflicts_with = "detection_threshold", hide = true)]
    sigma: Option<f64>,

    /// Ignored since 1.8.0: the star budget is --max-stars
    #[arg(long, value_name = "N", hide = true)]
    measure_cap: Option<usize>,

    /// Ignored since 1.8.0: a mosaic is measured on its green sites
    #[arg(long, hide = true)]
    no_debayer: bool,
}

/// `--psf-model`: the two models the Analysis tab offers.
#[derive(Clone, Copy, Debug, ValueEnum)]
enum PsfChoice {
    /// Moffat with beta fixed at 4 (the Analysis tab's default)
    Moffat4,
    /// The frame's own beta, fitted
    Auto,
}

impl PsfChoice {
    fn model(self) -> PsfModel {
        match self {
            PsfChoice::Moffat4 => PsfModel::Moffat4,
            PsfChoice::Auto => PsfModel::Auto,
        }
    }

    /// The name the flag takes and the `psf_model` column prints.
    fn name(self) -> &'static str {
        match self {
            PsfChoice::Moffat4 => "moffat4",
            PsfChoice::Auto => "auto",
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SortKey {
    /// File name
    Name,
    /// Most stars first
    Stars,
    /// Smallest FWHM first
    Fwhm,
    /// Smallest HFR first
    Hfr,
    /// Roundest stars first
    Ecc,
    /// Highest median star SNR first
    Snr,
    /// Highest SNR weight first
    Weight,
}

/// What the printers need besides the rows: the options the numbers were measured with, and the optics.
#[derive(Clone, Copy)]
struct Context {
    detection_threshold: f64,
    psf: PsfChoice,
    /// (focal length mm, pixel size µm) when both were given.
    optics: Option<(f64, f64)>,
}

type Row = (String, Measured);

pub fn run(args: AnalyzeArgs) -> Result<bool> {
    // Old scripts keep running; each flag is said once, even with --quiet.
    if args.sigma.is_some() {
        eprintln!("warning: --sigma is --detection-threshold since 1.8.0: another detector, default 6.0");
    }
    if args.measure_cap.is_some() {
        warn_ignored("--measure-cap", MEASURE_CAP_WHY);
    }
    if args.no_debayer {
        warn_ignored("--no-debayer", "a mosaic is measured on its green sites");
    }

    let inputs = batch::collect_inputs(&args.paths, args.recursive)?;
    let pool = batch::thread_pool(args.threads)?;
    let workers = batch::workers(args.jobs, DEFAULT_MAX_JOBS_ANALYSIS, inputs.len());
    let bar = batch::progress_bar(inputs.len(), args.quiet);

    let ctx = Context {
        detection_threshold: args.detection_threshold.or(args.sigma).unwrap_or(DEFAULT_THRESHOLD),
        psf: args.psf_model,
        optics: args.focal_length.zip(args.pixel_size),
    };
    let opts = frame_options(&ctx, args.max_stars);

    let results = batch::run_parallel(
        &inputs,
        workers,
        &bar,
        |input| file_name(&input.path),
        |input| measure_file(&input.path, &opts, &pool),
        |input, result| {
            if let Err(e) = result {
                report_file_error(&input.path, e);
            }
        },
    );

    let ok = results.iter().all(Result::is_ok);
    let mut rows: Vec<Row> = inputs
        .iter()
        .zip(results)
        .filter_map(|(input, r)| r.ok().map(|r| (input.path.display().to_string(), r)))
        .collect();
    if let Some(key) = args.sort {
        sort_rows(&mut rows, key);
    }
    if args.reverse {
        rows.reverse();
    }

    if args.json {
        print_json(&rows, &ctx);
    } else if args.csv {
        let (cols, cells) = full_precision(&rows, &ctx);
        print_csv(&cols, &cells);
    } else {
        print_summary(&rows, &ctx);
    }
    Ok(ok)
}

/// `FrameOptions` is `#[non_exhaustive]`: the fields are assigned onto its default.
fn frame_options(ctx: &Context, max_stars: usize) -> FrameOptions {
    let mut opts = FrameOptions::default();
    opts.detection_threshold = ctx.detection_threshold;
    opts.psf_model = ctx.psf.model();
    opts.max_stars = max_stars;
    opts
}

/// Best-first; a frame with a verdict (it has no medians) after every frame without one, in input order.
fn sort_rows(rows: &mut [Row], key: SortKey) {
    let metric = |m: &Measured| -> f64 {
        let s = &m.summary;
        match key {
            SortKey::Name => 0.0,
            SortKey::Stars => -(s.stars_detected as f64),
            SortKey::Fwhm => s.median_fwhm,
            SortKey::Hfr => s.median_hfr,
            SortKey::Ecc => s.median_eccentricity,
            SortKey::Snr => -s.median_snr,
            SortKey::Weight => -s.snr_weight,
        }
    };
    match key {
        SortKey::Name => rows.sort_by(|a, b| a.0.cmp(&b.0)),
        // `sort_by` is stable: two frames with a verdict compare equal and keep their input order.
        _ => rows.sort_by(|a, b| match (a.1.summary.verdict.is_some(), b.1.summary.verdict.is_some()) {
            (false, false) => metric(&a.1).total_cmp(&metric(&b.1)),
            (a_verdict, b_verdict) => a_verdict.cmp(&b_verdict),
        }),
    }
}

/// px -> arcsec at the plate scale (the legacy formula); `None` without optics or for a frame with a verdict.
fn arcsec(px: f64, m: &Measured, optics: Option<(f64, f64)>) -> Option<f64> {
    let (focal_mm, pixel_um) = optics?;
    if m.summary.verdict.is_some() || !(focal_mm > 0.0 && pixel_um > 0.0) {
        return None;
    }
    Some(px * (pixel_um / focal_mm * 206.265))
}

/// The table's VERDICT cell.
fn verdict_cell(v: Option<Verdict>) -> &'static str {
    match v {
        Some(Verdict::Trailed) => "trailed",
        Some(Verdict::Defocused) => "defocused",
        Some(Verdict::NoStars) => "no stars",
        None => "-",
    }
}

/// The `frame_verdict` value (Athenaeum's `FrameVerdict`, snake_case); `None` without a verdict.
fn verdict_key(v: Option<Verdict>) -> Option<&'static str> {
    match v {
        Some(Verdict::Trailed) => Some("trailed"),
        Some(Verdict::Defocused) => Some("defocused"),
        Some(Verdict::NoStars) => Some("no_stars"),
        None => None,
    }
}

fn fixed(v: f64, decimals: usize) -> String {
    if v.is_finite() {
        format!("{v:.decimals$}")
    } else {
        "-".to_string()
    }
}

fn print_summary(rows: &[Row], ctx: &Context) {
    let with_arcsec = ctx.optics.is_some();
    let mut cols = vec!["FILE", "STARS", "FWHM", "HFR"];
    if with_arcsec {
        cols.extend(["FWHM\"", "HFR\""]);
    }
    cols.extend(["ECC", "SNR", "WEIGHT", "BG", "NOISE", "VERDICT"]);
    let cells: Vec<Vec<String>> = rows
        .iter()
        .map(|(file, m)| {
            let s = &m.summary;
            // No fit stands in the frame statistic: the medians are 0.0 (PC-7), which is not a measurement.
            let metric = |text: String| if s.verdict.is_some() { "-".to_string() } else { text };
            let mut row = vec![
                file.clone(),
                s.stars_detected.to_string(),
                metric(fixed(s.median_fwhm, 2)),
                metric(fixed(s.median_hfr, 2)),
            ];
            if with_arcsec {
                row.push(arcsec(s.median_fwhm, m, ctx.optics).map_or("-".into(), |v| fixed(v, 2)));
                row.push(arcsec(s.median_hfr, m, ctx.optics).map_or("-".into(), |v| fixed(v, 2)));
            }
            row.extend([
                metric(fixed(s.median_eccentricity, 3)),
                metric(fixed(s.median_snr, 1)),
                metric(format_weight(s.snr_weight)),
                fixed(s.background, 0),
                fixed(s.noise, 1),
                verdict_cell(s.verdict).into(),
            ]);
            row
        })
        .collect();
    let right: Vec<bool> = cols.iter().map(|c| *c != "FILE" && *c != "VERDICT").collect();
    print_table(&cols, &cells, &right);
}

/// SNR weight spans orders of magnitude between setups; 4 significant digits.
fn format_weight(v: f64) -> String {
    if !v.is_finite() {
        return "-".into();
    }
    let digits = if v == 0.0 { 0 } else { (3 - v.abs().log10().floor() as i32).max(0) as usize };
    format!("{v:.digits$}")
}

/// The columns of the CSV and the keys of the JSON. The first 19 are the ones this command always had, in their
/// order, so a script that reads by position keeps working; the six after them are what `analysis_mode` adds.
const COLUMNS: [&str; 25] = [
    "file", "width", "height", "stars_detected", "stars_measured", "median_fwhm", "median_hfr",
    "median_eccentricity", "median_snr", "median_beta", "snr_weight", "psf_signal", "frame_snr",
    "background", "noise", "trail_r_squared", "possibly_trailed", "median_fwhm_arcsec", "median_hfr_arcsec",
    "channels", "stars_refused", "frame_verdict", "measure_version", "detection_threshold", "psf_model",
];

/// Every scalar metric, unrounded (`{}` of the f64), for CSV. A frame with a verdict keeps its 0.0 medians: that is
/// what Athenaeum stores for it.
fn full_precision(rows: &[Row], ctx: &Context) -> (Vec<&'static str>, Vec<Vec<String>>) {
    let opt = |v: Option<f64>| v.map_or(String::new(), |v| v.to_string());
    let cells = rows
        .iter()
        .map(|(file, m)| {
            let s = &m.summary;
            vec![
                file.clone(),
                m.frame_w.to_string(),
                m.frame_h.to_string(),
                s.stars_detected.to_string(),
                m.pm.stars.len().to_string(),
                s.median_fwhm.to_string(),
                s.median_hfr.to_string(),
                s.median_eccentricity.to_string(),
                s.median_snr.to_string(),
                opt(s.median_beta),
                s.snr_weight.to_string(),
                s.psf_signal.to_string(),
                s.frame_snr.to_string(),
                s.background.to_string(),
                s.noise.to_string(),
                s.trail_r_squared.to_string(),
                m.pm.trailed.to_string(),
                opt(arcsec(s.median_fwhm, m, ctx.optics)),
                opt(arcsec(s.median_hfr, m, ctx.optics)),
                m.channels.to_string(),
                s.stars_refused.to_string(),
                verdict_key(s.verdict).unwrap_or_default().to_string(),
                MEASURE_VERSION.to_string(),
                ctx.detection_threshold.to_string(),
                ctx.psf.name().to_string(),
            ]
        })
        .collect();
    (COLUMNS.to_vec(), cells)
}

fn print_json(rows: &[Row], ctx: &Context) {
    let opt = |v: Option<f64>| v.map_or("null".to_string(), json_num);
    println!("[");
    for (i, (file, m)) in rows.iter().enumerate() {
        let s = &m.summary;
        println!(
            "  {{\"file\": {}, \"width\": {}, \"height\": {}, \"stars_detected\": {}, \"stars_measured\": {}, \
             \"median_fwhm\": {}, \"median_hfr\": {}, \"median_eccentricity\": {}, \"median_snr\": {}, \
             \"median_beta\": {}, \"snr_weight\": {}, \"psf_signal\": {}, \"frame_snr\": {}, \
             \"background\": {}, \"noise\": {}, \"trail_r_squared\": {}, \"possibly_trailed\": {}, \
             \"median_fwhm_arcsec\": {}, \"median_hfr_arcsec\": {}, \"channels\": {}, \"stars_refused\": {}, \
             \"frame_verdict\": {}, \"measure_version\": {}, \"detection_threshold\": {}, \"psf_model\": {}}}{}",
            json_str(file),
            m.frame_w,
            m.frame_h,
            s.stars_detected,
            m.pm.stars.len(),
            json_num(s.median_fwhm),
            json_num(s.median_hfr),
            json_num(s.median_eccentricity),
            json_num(s.median_snr),
            opt(s.median_beta),
            json_num(s.snr_weight),
            json_num(s.psf_signal),
            json_num(s.frame_snr),
            json_num(s.background),
            json_num(s.noise),
            json_num(s.trail_r_squared),
            m.pm.trailed,
            opt(arcsec(s.median_fwhm, m, ctx.optics)),
            opt(arcsec(s.median_hfr, m, ctx.optics)),
            m.channels,
            s.stars_refused,
            verdict_key(s.verdict).map_or("null".to_string(), json_str),
            MEASURE_VERSION,
            json_num(ctx.detection_threshold),
            json_str(ctx.psf.name()),
            if i + 1 < rows.len() { "," } else { "" }
        );
    }
    println!("]");
}
