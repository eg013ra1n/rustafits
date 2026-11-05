//! Batch showcase: measure many frames in parallel with `star_measure::analysis_mode` and print one line per file.
//!
//!     cargo run -p rustafits --example analyze_batch --release -- tests/mono.fits tests/osc.fits
//!
//! With no paths it reads every `.fits` in `tests/no_trails`. Two files are measured at a time (a large frame is
//! hundreds of MB); each measurement runs on one shared pool, as `rustafits analyze` does.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use anyhow::{Context, Result};
use astroimage::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
use astroimage::star_measure::raw::read_raw_frame;
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};

/// Files measured at the same time.
const CONCURRENT_FILES: usize = 2;

/// One file's row: the summary plus the facts the summary does not carry.
struct Row {
    summary: FrameSummary,
    kept: usize,
    w: usize,
    h: usize,
    trailed: bool,
}

fn measure(path: &Path, opts: &FrameOptions, pool: &Arc<ThreadPool>) -> Result<Row> {
    let mut frame = read_raw_frame(path)?;
    frame.drop_unmeasured_planes();
    let pm = measure_raw_frame(&frame, opts, Some(pool));
    Ok(Row {
        summary: FrameSummary::from_measurement(&pm, opts),
        kept: pm.stars.len(),
        w: frame.w,
        h: frame.h,
        trailed: pm.trailed,
    })
}

/// The paths on the command line, else the `.fits` files of `tests/no_trails`, sorted.
fn input_paths() -> Result<Vec<PathBuf>> {
    let mut files: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    if files.is_empty() {
        let dir = "tests/no_trails";
        files = std::fs::read_dir(dir)
            .with_context(|| format!("cannot read {dir}"))?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.extension().is_some_and(|ext| ext == "fits"))
            .collect();
        files.sort();
    }
    Ok(files)
}

fn main() -> Result<ExitCode> {
    let files = input_paths()?;
    let opts = FrameOptions::default();
    let measure_pool = Arc::new(ThreadPoolBuilder::new().build().context("failed to create the measurement pool")?);
    let file_pool = ThreadPoolBuilder::new()
        .num_threads(CONCURRENT_FILES)
        .build()
        .context("failed to create the file pool")?;

    let rows: Vec<Result<Row>> =
        file_pool.install(|| files.par_iter().map(|p| measure(p, &opts, &measure_pool)).collect());

    println!(
        "{:<12} {:>5} {:>5} {:>6} {:>6} {:>6} {:>6} {:>7} {:>8} {:>8} {:>6} {:>6} {:>6} {:>7}",
        "FILE", "DET", "KEPT", "FWHM", "ECC", "SNR", "HFR", "SNR_W", "PSF_SIG", "FR_SNR", "W", "H", "R²", "TRAIL?"
    );
    println!("{}", "-".repeat(124));

    let mut failed = 0usize;
    for (path, row) in files.iter().zip(rows) {
        let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("?");
        // Short name: the frame number after the last '_'.
        let short = match stem.rfind('_') {
            Some(pos) => &stem[pos + 1..],
            None => &stem[..stem.len().min(12)],
        };
        match row {
            Ok(r) => {
                let s = &r.summary;
                println!(
                    "{:<12} {:>5} {:>5} {:>6.2} {:>6.3} {:>6.1} {:>6.2} {:>7.3} {:>8.1} {:>8.1} {:>6} {:>6} {:>6.3} {:>7}",
                    short,
                    s.stars_detected,
                    r.kept,
                    s.median_fwhm,
                    s.median_eccentricity,
                    s.median_snr,
                    s.median_hfr,
                    s.snr_weight,
                    s.psf_signal,
                    s.frame_snr,
                    r.w,
                    r.h,
                    s.trail_r_squared,
                    if r.trailed { "YES" } else { "no" },
                );
            }
            Err(e) => {
                failed += 1;
                println!("{:<12} ERROR", short);
                eprintln!("{}: {e:#}", path.display());
            }
        }
    }
    Ok(if failed == 0 { ExitCode::SUCCESS } else { ExitCode::from(1) })
}
