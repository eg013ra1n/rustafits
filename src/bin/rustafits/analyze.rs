//! `rustafits analyze` — star metrics per frame, for culling subs.
//!
//! The analyzer runs with the library defaults (what Athenaeum uses), so the
//! numbers match Athenaeum's for the same file; only the flags given override.

use std::path::PathBuf;

use anyhow::Result;
use astroimage::{AnalysisResult, ImageAnalyzer};
use clap::{Args, ValueEnum};

use crate::batch::{self, file_name};
use crate::convert::DEFAULT_MAX_JOBS_ANALYSIS;
use crate::report::{json_num, json_str, print_csv, print_table};

#[derive(Args, Debug)]
pub struct AnalyzeArgs {
    /// Input FITS/XISF files or directories
    #[arg(required = true, value_name = "INPUT")]
    paths: Vec<PathBuf>,

    /// Descend into subdirectories of directory inputs
    #[arg(short, long)]
    recursive: bool,

    /// Sort rows best-first by a metric (default: input order)
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

    /// Detection threshold in sigma above background [library default: 5]
    #[arg(long, value_name = "SIGMA")]
    sigma: Option<f32>,

    /// Max stars PSF-fitted per frame, 0 = all [library default: 500]
    #[arg(long, value_name = "N")]
    measure_cap: Option<usize>,

    /// Focal length in mm: adds FWHM/HFR in arcsec (with --pixel-size)
    #[arg(long, value_name = "MM", requires = "pixel_size")]
    focal_length: Option<f64>,

    /// Pixel size in µm (binned): adds FWHM/HFR in arcsec (with --focal-length)
    #[arg(long, value_name = "UM", requires = "focal_length")]
    pixel_size: Option<f64>,

    /// Analyze OSC frames on the raw mosaic (faster, less accurate)
    #[arg(long)]
    no_debayer: bool,

    /// Frames analyzed concurrently [default: min(2, cores, files)]
    #[arg(short, long, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    jobs: Option<u32>,

    /// Threads for pixel work, shared by all jobs [default: all cores]
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    threads: Option<u32>,

    /// Errors only: no progress bar
    #[arg(long)]
    quiet: bool,
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

pub fn run(args: AnalyzeArgs) -> Result<bool> {
    let inputs = batch::collect_inputs(&args.paths, args.recursive)?;
    let pool = batch::thread_pool(args.threads)?;
    let workers = batch::workers(args.jobs, DEFAULT_MAX_JOBS_ANALYSIS, inputs.len());
    let bar = batch::progress_bar(inputs.len(), args.quiet);

    let mut analyzer = ImageAnalyzer::new().with_thread_pool(pool.clone());
    if let Some(s) = args.sigma {
        analyzer = analyzer.with_detection_sigma(s);
    }
    if let Some(n) = args.measure_cap {
        analyzer = analyzer.with_measure_cap(n);
    }
    if let (Some(f), Some(p)) = (args.focal_length, args.pixel_size) {
        analyzer = analyzer.with_optics(f, p);
    }
    if args.no_debayer {
        analyzer = analyzer.without_debayer();
    }

    let results = batch::run_parallel(
        &inputs,
        workers,
        &bar,
        |input| file_name(&input.path),
        |input| analyzer.analyze(&input.path),
        |input, result| {
            if let Err(e) = result {
                eprintln!("error: {}: {e:#}", input.path.display());
            }
        },
    );

    let ok = results.iter().all(Result::is_ok);
    let mut rows: Vec<(String, AnalysisResult)> = inputs
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

    let arcsec = args.focal_length.is_some();
    if args.json {
        print_json(&rows);
    } else if args.csv {
        let (cols, cells) = full_precision(&rows);
        print_csv(&cols, &cells);
    } else {
        print_summary(&rows, arcsec);
    }
    Ok(ok)
}

/// Best-first; NaN (no stars measured) always last.
fn sort_rows(rows: &mut [(String, AnalysisResult)], key: SortKey) {
    let metric = |r: &AnalysisResult| -> f64 {
        match key {
            SortKey::Name => 0.0,
            SortKey::Stars => -(r.stars_detected as f64),
            SortKey::Fwhm => r.median_fwhm as f64,
            SortKey::Hfr => r.median_hfr as f64,
            SortKey::Ecc => r.median_eccentricity as f64,
            SortKey::Snr => -(r.median_snr as f64),
            SortKey::Weight => -(r.snr_weight as f64),
        }
    };
    match key {
        SortKey::Name => rows.sort_by(|a, b| a.0.cmp(&b.0)),
        _ => rows.sort_by(|a, b| {
            let (x, y) = (metric(&a.1), metric(&b.1));
            match (x.is_nan(), y.is_nan()) {
                (false, false) => x.total_cmp(&y),
                (a_nan, b_nan) => a_nan.cmp(&b_nan),
            }
        }),
    }
}

fn fixed(v: f32, decimals: usize) -> String {
    if v.is_finite() {
        format!("{v:.decimals$}")
    } else {
        "-".to_string()
    }
}

fn print_summary(rows: &[(String, AnalysisResult)], arcsec: bool) {
    let mut cols = vec!["FILE", "STARS", "FWHM", "HFR"];
    if arcsec {
        cols.extend(["FWHM\"", "HFR\""]);
    }
    cols.extend(["ECC", "SNR", "WEIGHT", "BG", "NOISE", "TRAIL"]);
    let cells: Vec<Vec<String>> = rows
        .iter()
        .map(|(file, r)| {
            let mut row = vec![file.clone(), r.stars_detected.to_string(), fixed(r.median_fwhm, 2), fixed(r.median_hfr, 2)];
            if arcsec {
                row.push(r.median_fwhm_arcsec.map_or("-".into(), |v| fixed(v, 2)));
                row.push(r.median_hfr_arcsec.map_or("-".into(), |v| fixed(v, 2)));
            }
            row.extend([
                fixed(r.median_eccentricity, 3),
                fixed(r.median_snr, 1),
                format_weight(r.snr_weight),
                fixed(r.background, 0),
                fixed(r.noise, 1),
                if r.possibly_trailed { "yes".into() } else { "-".into() },
            ]);
            row
        })
        .collect();
    let right: Vec<bool> = cols.iter().map(|c| *c != "FILE" && *c != "TRAIL").collect();
    print_table(&cols, &cells, &right);
}

/// SNR weight spans orders of magnitude between setups; 4 significant digits.
fn format_weight(v: f32) -> String {
    if !v.is_finite() {
        return "-".into();
    }
    let digits = if v == 0.0 { 0 } else { (3 - v.abs().log10().floor() as i32).max(0) as usize };
    format!("{v:.digits$}")
}

/// Every scalar metric, unrounded (`{}` of the f32), for CSV.
fn full_precision(rows: &[(String, AnalysisResult)]) -> (Vec<&'static str>, Vec<Vec<String>>) {
    let opt = |v: Option<f32>| v.map_or(String::new(), |v| v.to_string());
    let cols = vec![
        "file", "width", "height", "stars_detected", "stars_measured", "median_fwhm", "median_hfr",
        "median_eccentricity", "median_snr", "median_beta", "snr_weight", "psf_signal", "frame_snr",
        "background", "noise", "trail_r_squared", "possibly_trailed", "median_fwhm_arcsec",
        "median_hfr_arcsec",
    ];
    let cells = rows
        .iter()
        .map(|(file, r)| {
            vec![
                file.clone(),
                r.width.to_string(),
                r.height.to_string(),
                r.stars_detected.to_string(),
                r.stars_measured.to_string(),
                r.median_fwhm.to_string(),
                r.median_hfr.to_string(),
                r.median_eccentricity.to_string(),
                r.median_snr.to_string(),
                opt(r.median_beta),
                r.snr_weight.to_string(),
                r.psf_signal.to_string(),
                r.frame_snr.to_string(),
                r.background.to_string(),
                r.noise.to_string(),
                r.trail_r_squared.to_string(),
                r.possibly_trailed.to_string(),
                opt(r.median_fwhm_arcsec),
                opt(r.median_hfr_arcsec),
            ]
        })
        .collect();
    (cols, cells)
}

fn print_json(rows: &[(String, AnalysisResult)]) {
    let num = |v: f32| json_num(v as f64);
    let opt = |v: Option<f32>| v.map_or("null".to_string(), |v| json_num(v as f64));
    println!("[");
    for (i, (file, r)) in rows.iter().enumerate() {
        println!(
            "  {{\"file\": {}, \"width\": {}, \"height\": {}, \"stars_detected\": {}, \"stars_measured\": {}, \
             \"median_fwhm\": {}, \"median_hfr\": {}, \"median_eccentricity\": {}, \"median_snr\": {}, \
             \"median_beta\": {}, \"snr_weight\": {}, \"psf_signal\": {}, \"frame_snr\": {}, \
             \"background\": {}, \"noise\": {}, \"trail_r_squared\": {}, \"possibly_trailed\": {}, \
             \"median_fwhm_arcsec\": {}, \"median_hfr_arcsec\": {}}}{}",
            json_str(file),
            r.width,
            r.height,
            r.stars_detected,
            r.stars_measured,
            num(r.median_fwhm),
            num(r.median_hfr),
            num(r.median_eccentricity),
            num(r.median_snr),
            opt(r.median_beta),
            num(r.snr_weight),
            num(r.psf_signal),
            num(r.frame_snr),
            num(r.background),
            num(r.noise),
            num(r.trail_r_squared),
            r.possibly_trailed,
            opt(r.median_fwhm_arcsec),
            opt(r.median_hfr_arcsec),
            if i + 1 < rows.len() { "," } else { "" }
        );
    }
    println!("]");
}
