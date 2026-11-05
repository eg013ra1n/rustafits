//! What `analyze` and `convert --annotate` share: one measurement of a file on `star_measure::analysis_mode`
//! (Athenaeum's Analysis numbers), the range checks of the two options both commands take, and the notice an option
//! kept for old scripts prints. One copy of each, so the two commands cannot drift apart.

use std::path::Path;
use std::sync::Arc;

use anyhow::Result;
use astroimage::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
use astroimage::star_measure::raw::read_raw_frame;
use astroimage::star_measure::seeds::{MAX_THRESHOLD, MIN_THRESHOLD};
use astroimage::star_measure::PlaneMeasurement;
use astroimage::ThreadPool;

/// The star budget's range, as Athenaeum's `AnalysisConfig::validate` has it.
const MIN_BUDGET: usize = 10;
const MAX_BUDGET: usize = 2000;

/// Why `--measure-cap` does nothing, for the one notice both commands print (see [`warn_ignored`]).
pub(crate) const MEASURE_CAP_WHY: &str = "the star budget is --max-stars";

/// One file, measured: what the commands read off it. The planes are not kept (a 26 MP frame is hundreds of MB), only
/// the frame's geometry the star positions need.
pub(crate) struct Measured {
    /// The native size: a mosaic is measured, and its stars positioned, on this grid.
    pub frame_w: usize,
    pub frame_h: usize,
    pub channels: usize,
    /// The file is bottom-up; star positions are in file row order, so a render flips them.
    // Read by `convert --annotate` (a later change); `analyze` prints the table, not a picture.
    #[allow(dead_code)]
    pub flip_vertical: bool,
    pub pm: PlaneMeasurement,
    pub summary: FrameSummary,
}

/// Read `path`, measure it with `opts` on `pool` and summarise it.
pub(crate) fn measure_file(path: &Path, opts: &FrameOptions, pool: &Arc<ThreadPool>) -> Result<Measured> {
    let mut frame = read_raw_frame(path)?;
    frame.drop_unmeasured_planes();
    let pm = measure_raw_frame(&frame, opts, Some(pool));
    let summary = FrameSummary::from_measurement(&pm, opts);
    Ok(Measured {
        frame_w: frame.w,
        frame_h: frame.h,
        channels: frame.channels,
        flip_vertical: frame.flip_vertical,
        pm,
        summary,
    })
}

/// `--max-stars`: the budget, 10 to 2000.
pub(crate) fn budget_parser(s: &str) -> Result<usize, String> {
    let n: usize = s.parse().map_err(|e| format!("{e}"))?;
    if (MIN_BUDGET..=MAX_BUDGET).contains(&n) {
        Ok(n)
    } else {
        Err(format!("must be between {MIN_BUDGET} and {MAX_BUDGET}"))
    }
}

/// `--detection-threshold`: sigma, 2 to 20 (`star_measure::seeds`). NaN is out of range.
pub(crate) fn threshold_parser(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|e| format!("{e}"))?;
    if (MIN_THRESHOLD..=MAX_THRESHOLD).contains(&v) {
        Ok(v)
    } else {
        Err(format!("must be between {MIN_THRESHOLD} and {MAX_THRESHOLD}"))
    }
}

/// The notice for a flag that is still accepted but does nothing. It goes to stderr whatever `--quiet` says: the user
/// asked for something that did not happen.
pub(crate) fn warn_ignored(flag: &str, why: &str) {
    eprintln!("warning: {flag} has no effect since 1.8.0: {why}");
}
