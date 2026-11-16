//! The Analysis tab's measurement: a raw frame measured through the one star path (Plan C, Task C5).
//!
//! [`measure_raw_frame`] runs `star_measure` on the frame (mono, a Bayer mosaic's green sites, or a debayered file's G
//! plane); [`FrameSummary::from_measurement`] reduces the result to the frame figures the Analysis tab shows. Positions
//! stay in FILE row order (unflipped); [`super::raw::RawFrame::flip_vertical`] tells a display to flip them.
//!
//! This module and [`super::raw`] are the stable entry of `star_measure`: a caller reads a file with
//! [`super::raw::read_raw_frame`], measures it with [`measure_raw_frame`] and reads the figures off the
//! [`FrameSummary`].
//!
//! What that promise covers: the function names and what [`FrameOptions`] and [`FrameSummary`] carry, which are
//! `#[non_exhaustive]` so that fields can be added in 1.x. [`PlaneMeasurement`] (the return type of
//! [`measure_raw_frame`]) and [`PsfModel`] (the type of `FrameOptions::psf_model`) appear in the signatures but keep
//! the status of `star_measure`'s other public items, outside the stable entry.
//!
//! A list entry sits beside it ([`measure_raw_frame_and_list`], forced-list zero point): the same measurement, then a
//! caller's list of native-pixel positions fitted on the same plane with the measurement's own β, field width and mask,
//! and an aperture correction `c` measured over the list's own bright members. The list adds to the measurement and
//! never moves it: every [`StarFit`], [`FrameSummary`] and [`PlaneMeasurement`] field but the wall-clock `timings` /
//! `sub_steps` is bit-identical with and without a list. [`ListRequest`], [`ListFit`] and [`ListMeasurement`] are
//! `#[non_exhaustive]` like [`FrameOptions`].
//!
//! The row mapping (the per-star rows, the frame row's ids, `config_hash`, `analyzed_at`,
//! `quality_score`, the detection threshold, the frame's size and channel count, `possibly_trailed` and the measurement
//! version) stays in Athenaeum, which builds a [`FrameOptions`] from its own `AnalysisConfig`.
//!
//! **In rustafits** (port plan phase 4, `docs/superpowers/plans/2026-11-04-star-measure-port-phase-4.md`): copied from
//! Athenaeum's `analysis/analyzer.rs` at e6344ddd. `AnalysisConfig` becomes [`FrameOptions`], and what Athenaeum's
//! `map_measurement` computed per frame becomes [`FrameSummary`], statement for statement and in the same order. The
//! rulings the docs cite (`PC-…`, `R-PC-…`) are Athenaeum's.

// Athenaeum's test helper `site` as copied (port plan G-3): the clippy style lint it trips is allowed, not "fixed".
#![allow(clippy::manual_is_multiple_of)]

use std::sync::Arc;

use crate::star_measure::frame::finite_median;
use crate::star_measure::green::{measure_green_sites, measure_green_sites_and_seeds};
use crate::star_measure::list::{measure_list, ListPlane};
use crate::star_measure::raw::{RawFrame, RawPlanes};
use crate::star_measure::trail::trail_r_squared;
use crate::star_measure::{
    self, BetaSource, MeasureOptions, Mode, PixelNoise, PlaneMeasurement, PsfModel, RefusalCounts,
    Seed, SeedSpec, StarFit, StarFlags, DEFAULT_THRESHOLD,
};

/// What the Analysis-mode measurement is asked: the three knobs of Athenaeum's `AnalysisConfig` that reach the
/// measurement. The defaults are the Analysis tab's.
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
    /// 6.0 σ, Moffat 4, 500 stars. (`PsfModel::default()` is `Auto`, so the model is named here.)
    fn default() -> Self {
        FrameOptions {
            detection_threshold: DEFAULT_THRESHOLD,
            psf_model: PsfModel::Moffat4,
            max_stars: 500,
        }
    }
}

/// Measure a raw frame (PC-4, PC-5): mono -> `measure_plane`; a Bayer mosaic -> [`measure_green_sites`]; a
/// debayered file -> `measure_plane` on its G plane (`colour`, mixed pixel noise, green share 0.5).
pub fn measure_raw_frame(
    frame: &RawFrame,
    opts: &FrameOptions,
    pool: Option<&Arc<rayon::ThreadPool>>,
) -> PlaneMeasurement {
    measure_raw_frame_with(frame, opts, pool, false)
}

/// [`measure_raw_frame`] with the module's per-seed diagnostics (`diagnostics`) for the probe's refusal lines; the
/// numbers are the same.
pub fn measure_raw_frame_with(
    frame: &RawFrame,
    opts: &FrameOptions,
    pool: Option<&Arc<rayon::ThreadPool>>,
    diagnostics: bool,
) -> PlaneMeasurement {
    let base = base_options(frame, opts, diagnostics);
    match &frame.planes {
        RawPlanes::Mono(plane) => {
            let opts = MeasureOptions {
                pixel_noise: PixelNoise::White,
                ..base
            };
            star_measure::measure_plane(plane, frame.w, frame.h, &opts, pool)
        }
        RawPlanes::Cfa { mosaic, pattern } => {
            measure_green_sites(mosaic, frame.w, frame.h, *pattern, &base, pool)
        }
        RawPlanes::Rgb(planes) => {
            let opts = rgb_options(base);
            star_measure::measure_plane(&planes[1], frame.w, frame.h, &opts, pool)
        }
    }
}

fn base_options(frame: &RawFrame, opts: &FrameOptions, diagnostics: bool) -> MeasureOptions {
    MeasureOptions {
        domain: frame.domain,
        psf_model: opts.psf_model,
        seeds: SeedSpec::MatchedFilter {
            threshold: opts.detection_threshold,
        },
        mode: Mode::Budget {
            budget: opts.max_stars,
        },
        diagnostics,
        ..Default::default()
    }
}

fn rgb_options(base: MeasureOptions) -> MeasureOptions {
    MeasureOptions {
        colour: true,
        pixel_noise: PixelNoise::Mixed,
        cfa_site_share: Some(0.5),
        ..base
    }
}

/// The default cap of [`ListRequest::max_c_stars`]: the aperture correction is the median over at most this many of the
/// list's kept fits, the brightest by model flux.
pub const LIST_C_STARS: usize = 256;

/// A list of positions to fit beside the Analysis measurement ([`measure_raw_frame_and_list`]).
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ListRequest {
    /// Native pixel positions in FILE row order, pixel centres at integer coordinates (the positions
    /// [`PlaneMeasurement::stars`] carries).
    pub positions: Vec<(f64, f64)>,
    /// A fit is kept when its fitted centre is within this many px of its position. Must be positive.
    pub max_centre_offset_px: f64,
    /// The aperture correction is the median over at most this many kept fits ([`LIST_C_STARS`] from
    /// [`ListRequest::new`]). Must be positive.
    pub max_c_stars: usize,
}

impl ListRequest {
    /// A request for `positions` with the centre tolerance `max_centre_offset_px` and the default `c` cap
    /// ([`LIST_C_STARS`]).
    pub fn new(positions: Vec<(f64, f64)>, max_centre_offset_px: f64) -> Self {
        ListRequest {
            positions,
            max_centre_offset_px,
            max_c_stars: LIST_C_STARS,
        }
    }
}

/// One list position's fit, in request order. A position that was refused (outside the plane, not finite, or a fit the
/// module refused) carries NaN numbers and both flags false.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct ListFit {
    /// The fitted centre, native pixels.
    pub x: f64,
    pub y: f64,
    /// The fit's model flux inside its FWTM (what a star's `signal` is at `c = 1`).
    pub model_flux: f64,
    /// `model_flux` times the list's aperture correction; NaN when that is `None`.
    pub flux: f64,
    /// The star's aperture SNR (NaN without a ring sky).
    pub snr: f64,
    pub flags: StarFlags,
    /// The fit was accepted (no refusal, not a duplicate, not a blend).
    pub accepted: bool,
    /// Accepted, not saturated or near the clip, and centred within the request's tolerance. No SNR cut.
    pub kept: bool,
    /// Distance of the fitted centre from the position, px.
    pub centre_offset_px: f64,
}

/// Why [`measure_raw_frame_and_list`] measured no list on a frame whose measurement stands, when that is the frame's
/// layout and not a verdict of its plane (a trailed or defocused plane is `None` with the plane's own verdict).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ListRefusal {
    /// A debayered three-plane file: the list is fitted on a mono plane or a mosaic's green sites only.
    ColourFrame,
}

impl ListRefusal {
    /// A stable lower-case code the caller can store.
    pub fn code(self) -> &'static str {
        match self {
            ListRefusal::ColourFrame => "colour_frame",
        }
    }
}

/// What the list measured, apart from the frame's own measurement.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ListMeasurement {
    /// One per requested position, in request order; empty when `refusal` is set.
    pub fits: Vec<ListFit>,
    /// The list's aperture correction: the median growth / model-flux ratio over the kept fits that pass Analysis' own
    /// `c` screens (isolation, footprint, clip, SNR 20), at most the request's cap; the masked, phase-matched growth on
    /// a mosaic. `None` under 5 qualifiers.
    pub aperture_correction: Option<f64>,
    /// The qualifiers `aperture_correction` is over; the count (0 to 4) when it is `None`.
    pub aperture_stars: usize,
    /// The standard error of `aperture_correction`, from its ratios' MAD.
    pub aperture_se: Option<f64>,
    /// The β the list was fitted at (the frame's).
    pub beta: f64,
    /// The ring width `F` the list's ring skies and growth were read at; `None` without one.
    pub ring_fwhm: Option<f64>,
    /// The list's seeds not accepted, by reason (the positions refused without a fit are not counted).
    pub refused: RefusalCounts,
    /// Set when no list was measured because of the frame's layout ([`ListRefusal`]); `fits` is then empty.
    pub refusal: Option<ListRefusal>,
}

/// [`measure_raw_frame`], then the caller's list fitted on the same plane. The measurement is exactly
/// `measure_raw_frame(frame, opts, pool)`'s (every number bit-identical but the wall-clock `timings` / `sub_steps`).
/// `list` is called at most once, with that measurement, and returns the request built from it (the catalogue stars
/// projected into the frame, say); it is not called on a trailed or defocused plane or on a three-plane file.
///
/// The second value is
/// - `Some` with empty `fits` and `refusal == Some(ListRefusal::ColourFrame)` for a debayered three-plane file (no list
///   is measured on a colour image; the measurement is the frame's own, as ever);
/// - `None` when the plane is trailed or defocused (the callback is not called): the plane's verdict is the measurement's
///   own (`pm.trailed` / `pm.defocused`, or [`FrameSummary::verdict`]) and is what a caller stores as the reason;
/// - `None` when the callback returns `None` or an empty request, or when the request or the plane's mask cannot be
///   measured (each says why in the log).
///
/// The list is fitted at the measurement's β and field width with its saturation verdict and the plane's mask (a
/// mosaic's green sites); the frame's detected seeds, except one within the tolerance of a list position, keep their
/// `3F` discs out of every list star's ring sky. Cost: one more fit pass over the list and at most `max_c_stars` growth
/// curves beside the budget's.
pub fn measure_raw_frame_and_list<F>(
    frame: &RawFrame,
    opts: &FrameOptions,
    pool: Option<&Arc<rayon::ThreadPool>>,
    list: F,
) -> (PlaneMeasurement, Option<ListMeasurement>)
where
    F: FnOnce(&PlaneMeasurement) -> Option<ListRequest>,
{
    let base = base_options(frame, opts, false);
    // The measurement, as `measure_raw_frame_with(.., false)` makes it, plus what the list needs of it: the plane it ran
    // on, its mask and the seeds that stay in the ring sky.
    let (pm, plane, mask, detected_native): (
        PlaneMeasurement,
        &[f32],
        Option<Vec<bool>>,
        Option<Vec<Seed>>,
    ) = match &frame.planes {
        RawPlanes::Mono(plane) => {
            let o = MeasureOptions {
                pixel_noise: PixelNoise::White,
                ..base
            };
            (
                star_measure::measure_plane(plane, frame.w, frame.h, &o, pool),
                plane,
                None,
                None,
            )
        }
        RawPlanes::Cfa { mosaic, pattern } => {
            let (pm, seeds) =
                measure_green_sites_and_seeds(mosaic, frame.w, frame.h, *pattern, &base, pool);
            let mask = crate::processing::cfa::green_mask(frame.w, frame.h, *pattern);
            (pm, mosaic.as_slice(), Some(mask), Some(seeds))
        }
        RawPlanes::Rgb(planes) => {
            let o = rgb_options(base);
            (
                star_measure::measure_plane(&planes[1], frame.w, frame.h, &o, pool),
                planes[1].as_slice(),
                None,
                None,
            )
        }
    };
    if pm.trailed || pm.defocused {
        tracing::debug!(
            trailed = pm.trailed,
            defocused = pm.defocused,
            "list: the plane is trailed or defocused; no list measured"
        );
        return (pm, None);
    }
    if matches!(frame.planes, RawPlanes::Rgb(_)) {
        tracing::debug!("list: a three-plane file; no list measured");
        return (
            pm,
            Some(ListMeasurement {
                fits: Vec::new(),
                aperture_correction: None,
                aperture_stars: 0,
                aperture_se: None,
                beta: f64::NAN,
                ring_fwhm: None,
                refused: RefusalCounts::default(),
                refusal: Some(ListRefusal::ColourFrame),
            }),
        );
    }
    let Some(request) = list(&pm) else {
        return (pm, None);
    };
    let mosaic = matches!(frame.planes, RawPlanes::Cfa { .. });
    let input = ListPlane {
        plane,
        w: frame.w,
        h: frame.h,
        domain: frame.domain,
        mask: mask.as_deref(),
        mosaic,
    };
    let detected: &[Seed] = detected_native.as_deref().unwrap_or(&pm.detected);
    let measured = measure_list(&input, &pm, detected, &request, pool);
    (pm, measured)
}

/// The frame verdict of a measurement (PC-13): trailed / defocused as the module judged, else no stars when no fit stands
/// in the frame statistic. A frame with measured stars has none.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Verdict {
    Trailed,
    Defocused,
    NoStars,
}

/// The frame figures of an Analysis-mode measurement (Athenaeum's spec 4.2, rulings PC-6…PC-9, PC-11, PC-13).
///
/// Medians come from the module's frame aggregates and read 0.0 when the frame has no fit in its statistic (PC-7);
/// `snr_weight`, `psf_signal` and `trail_r_squared` read the same fits. `median_beta` is the frame's fitted beta under
/// `Auto` and `None` otherwise (PC-9). The medians, `snr_weight`, `psf_signal`, `background`, `noise` and
/// `trail_r_squared` are written 0.0 when non-finite; `frame_snr` (a ratio of the finite `background` and `noise`) and
/// `median_beta` are not wrapped, as in Athenaeum's `map_measurement`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct FrameSummary {
    /// The matched filter's seeds not flagged `extended` (R-PC-5 / R-PC-8), every layout.
    pub stars_detected: usize,
    /// Frame aggregates, 0.0 when no fit stands in the frame statistic (PC-7).
    pub median_fwhm: f64,
    pub median_eccentricity: f64,
    pub median_snr: f64,
    pub median_hfr: f64,
    /// background / noise, a ratio (PC-11); 0.0 without noise.
    pub frame_snr: f64,
    /// flux² / (noise² · background) over the in-statistic fits' median flux.
    pub snr_weight: f64,
    /// median amplitude / noise.
    pub psf_signal: f64,
    /// The sky median and the noise σ, non-finite written 0.0.
    pub background: f64,
    pub noise: f64,
    /// Rayleigh R² of the in-statistic fits' θ ([`super::trail::trail_r_squared`]), on the UNZEROED field FWHM.
    pub trail_r_squared: f64,
    /// The fitted β under `PsfModel::Auto` when its source is `Fitted`, else `None` (PC-9).
    pub median_beta: Option<f64>,
    /// The refusal reasons' sum; the budget's not-selected / trimmed are not in it.
    pub stars_refused: usize,
    /// Trailed / defocused as the module judged, else no stars when no fit stands in the statistic (PC-13).
    pub verdict: Option<Verdict>,
}

/// A non-finite number written as 0.0, so no figure is NaN.
fn finite_or_zero(x: f64) -> f64 {
    if x.is_finite() {
        x
    } else {
        0.0
    }
}

/// `stars_detected` (R-PC-5 / R-PC-8), one rule for every layout: the matched filter's seeds not flagged `extended`.
/// `seeds` is the detection's seed count on the module's own sources and the handed count on the CFA path.
fn point_like_seeds(pm: &PlaneMeasurement) -> usize {
    pm.seeds.saturating_sub(pm.extended_seeds)
}

impl FrameSummary {
    /// Reduce a measurement to the frame figures (spec 4.2 with PC-6..PC-9). `opts.psf_model` says whether `median_beta`
    /// is the fitted one (`Auto`) or `None`.
    pub fn from_measurement(pm: &PlaneMeasurement, opts: &FrameOptions) -> FrameSummary {
        let stat: Vec<&StarFit> = pm.stars.iter().filter(|f| f.flags.in_statistic()).collect();
        let background = finite_or_zero(pm.sky.median);
        let noise = finite_or_zero(pm.noise.sigma);
        let frame_snr = if noise > 0.0 { background / noise } else { 0.0 };
        let (snr_weight, psf_signal) = if noise > 0.0 && !stat.is_empty() {
            let flux = finite_median(stat.iter().map(|f| f.signal)).unwrap_or(0.0);
            let amp = finite_median(stat.iter().map(|f| f.amplitude)).unwrap_or(0.0);
            let w = if background > 0.0 {
                flux * flux / (noise * noise * background)
            } else {
                0.0
            };
            (finite_or_zero(w), finite_or_zero(amp / noise))
        } else {
            (0.0, 0.0)
        };
        let median_fwhm = pm.frame.fwhm_px.unwrap_or(0.0);
        let thetas: Vec<f64> = stat.iter().map(|f| f.theta).collect();
        let median_beta = match (opts.psf_model, pm.beta_source) {
            (PsfModel::Auto, BetaSource::Fitted { .. }) => Some(pm.beta),
            _ => None,
        };

        let verdict = if pm.trailed {
            Some(Verdict::Trailed)
        } else if pm.defocused {
            Some(Verdict::Defocused)
        } else if stat.is_empty() {
            Some(Verdict::NoStars)
        } else {
            None
        };

        FrameSummary {
            stars_detected: point_like_seeds(pm),
            median_fwhm: finite_or_zero(median_fwhm),
            median_eccentricity: finite_or_zero(pm.frame.eccentricity.unwrap_or(0.0)),
            median_snr: finite_or_zero(pm.frame.snr_median.unwrap_or(0.0)),
            median_hfr: finite_or_zero(pm.frame.hfr_median.unwrap_or(0.0)),
            frame_snr,
            snr_weight,
            psf_signal,
            background,
            noise,
            trail_r_squared: finite_or_zero(trail_r_squared(&thetas, median_fwhm)),
            median_beta,
            stars_refused: pm.refused.total(),
            verdict,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::star_measure::Domain;
    use crate::test_support::{add_noise, moffat_field, MoffatStar};
    use crate::BayerPattern;

    const RAW: Domain = Domain::RawAdu {
        clip: 65535.0,
        full_scale: 65535.0,
        gain_e_per_adu: None,
    };

    fn star(x: f64, y: f64, amp: f64, fwhm: f64) -> MoffatStar {
        // Moffat beta 4: FWHM = 2 alpha sqrt(2^(1/4) - 1).
        let alpha = fwhm / (2.0 * (2f64.powf(0.25) - 1.0).sqrt());
        MoffatStar {
            x,
            y,
            amp,
            alpha_x: alpha,
            alpha_y: alpha,
            theta: 0.0,
        }
    }

    /// `n` stars on a jittered grid of a `w` x `h` plane, amplitudes 500..2000.
    fn grid_stars(n: usize, w: usize, h: usize, fwhm: f64, seed: u64) -> Vec<MoffatStar> {
        let mut rng = crate::geometry::ransac::SplitMix64(seed);
        let cols = ((n as f64 * w as f64 / h as f64).sqrt().ceil() as usize).max(1);
        let rows = n.div_ceil(cols);
        let (cw, ch) = (w as f64 / cols as f64, h as f64 / rows as f64);
        (0..n)
            .map(|i| {
                let (c, r) = ((i % cols) as f64, (i / cols) as f64);
                let x = c * cw + 12.0 + rng.next_f64() * (cw - 24.0);
                let y = r * ch + 12.0 + rng.next_f64() * (ch - 24.0);
                star(x, y, 500.0 + rng.next_f64() * 1500.0, fwhm)
            })
            .collect()
    }

    fn mono_frame(n: usize, w: usize, h: usize, sky: f32, seed: u64) -> RawFrame {
        let stars = grid_stars(n, w, h, 3.5, seed);
        let mut p = moffat_field(w, h, &stars, 4.0, sky);
        add_noise(&mut p, 5.0, seed + 1);
        RawFrame {
            planes: RawPlanes::Mono(p),
            w,
            h,
            domain: RAW,
            flip_vertical: false,
            channels: 1,
        }
    }

    fn mono(p: Vec<f32>, w: usize, h: usize) -> RawFrame {
        RawFrame {
            planes: RawPlanes::Mono(p),
            w,
            h,
            domain: RAW,
            flip_vertical: false,
            channels: 1,
        }
    }

    /// The measurement and its summary, the two steps the Analysis tab takes on a frame.
    fn summary(f: &RawFrame, o: &FrameOptions) -> FrameSummary {
        FrameSummary::from_measurement(&measure_raw_frame(f, o, None), o)
    }

    fn all_finite(a: &FrameSummary) {
        for (name, v) in [
            ("median_fwhm", a.median_fwhm),
            ("median_eccentricity", a.median_eccentricity),
            ("median_snr", a.median_snr),
            ("median_hfr", a.median_hfr),
            ("frame_snr", a.frame_snr),
            ("snr_weight", a.snr_weight),
            ("psf_signal", a.psf_signal),
            ("background", a.background),
            ("noise", a.noise),
            ("trail_r_squared", a.trail_r_squared),
        ] {
            assert!(v.is_finite(), "{name} = {v}");
        }
    }

    #[test]
    fn frame_options_default_is_the_analysis_tab_default() {
        let o = FrameOptions::default();
        assert_eq!(o.detection_threshold, 6.0);
        assert_eq!(o.detection_threshold, DEFAULT_THRESHOLD);
        assert_eq!(o.psf_model, PsfModel::Moffat4);
        assert_eq!(o.max_stars, 500);
    }

    #[test]
    fn a_synthetic_mono_frame_summarizes_every_field() {
        let frame = mono_frame(300, 1024, 768, 1000.0, 7);
        let o = FrameOptions::default();
        let pm = measure_raw_frame(&frame, &o, None);
        let a = FrameSummary::from_measurement(&pm, &o);
        all_finite(&a);
        assert!(
            (a.median_fwhm - 3.5).abs() < 0.07,
            "median_fwhm {}",
            a.median_fwhm
        );
        assert!(a.median_eccentricity < 0.1, "ecc {}", a.median_eccentricity);
        assert!(a.stars_detected >= 285, "detected {}", a.stars_detected);
        assert!(!pm.stars.is_empty());
        assert_eq!(a.frame_snr, a.background / a.noise);
        assert!((a.background - 1000.0).abs() < 10.0, "bg {}", a.background);
        assert!(a.snr_weight > 0.0 && a.psf_signal > 0.0);
        assert!(a.median_beta.is_none());
        assert!(!pm.trailed);
        for s in &pm.stars {
            assert!(s.fwhm_x >= s.fwhm_y && s.eccentricity() >= 0.0 && s.eccentricity() < 1.0);
            assert!(s.x.is_finite() && s.y.is_finite() && s.amplitude > 0.0 && s.signal > 0.0);
        }
    }

    /// R-PC-5 / R-PC-8: `stars_detected` counts the matched filter's point-like seeds. A 3F Gaussian knot is a seed the
    /// module flags `extended`: it is fitted but is not counted, in every layout (mono, a mosaic's green sites, a
    /// debayered file's G plane), so the count with the knot equals the count without it.
    #[test]
    fn stars_detected_leaves_knots_out_in_every_layout() {
        use crate::test_support::gaussian_field;
        let (w, h) = (320, 240);
        let stars = grid_stars(40, w, h, 3.5, 21);
        let knot = (160.3, 120.8);
        let g = gaussian_field(
            w,
            h,
            &[(knot.0, knot.1, 400.0)],
            3.0 * 3.5 / 2.354_820_045,
            0.0,
        );
        let plane = |with_knot: bool, seed: u64| {
            let mut p = moffat_field(w, h, &stars, 4.0, 1000.0);
            if with_knot {
                for (d, b) in p.iter_mut().zip(&g) {
                    *d += b;
                }
            }
            add_noise(&mut p, 6.0, seed);
            p
        };
        let frame = |planes: RawPlanes, channels| RawFrame {
            planes,
            w,
            h,
            domain: RAW,
            flip_vertical: false,
            channels,
        };
        let make = |name: &str, k: bool| match name {
            "mono" => frame(RawPlanes::Mono(plane(k, 5)), 1),
            "cfa" => frame(
                RawPlanes::Cfa {
                    mosaic: plane(k, 5),
                    pattern: BayerPattern::Rggb,
                },
                1,
            ),
            _ => frame(
                RawPlanes::Rgb([plane(false, 3), plane(k, 5), plane(false, 4)]),
                3,
            ),
        };
        let o = FrameOptions::default();
        for name in ["mono", "cfa", "rgb"] {
            let (plain, knotted) = (make(name, false), make(name, true));
            let pm = measure_raw_frame(&knotted, &o, None);
            assert!(
                pm.extended_seeds >= 1,
                "{name}: the knot is flagged extended"
            );
            if name == "cfa" {
                assert!(
                    pm.detected.is_empty(),
                    "{name}: handed seeds leave `detected` empty"
                );
            } else {
                assert_eq!(pm.detected.len(), pm.seeds, "{name}");
            }
            let a = summary(&plain, &o);
            let b = summary(&knotted, &o);
            assert_eq!(b.stars_detected, pm.seeds - pm.extended_seeds, "{name}");
            assert_eq!(
                b.stars_detected, a.stars_detected,
                "{name}: the knot is not counted"
            );
            assert!(b.stars_detected >= 30, "{name}: {}", b.stars_detected);
        }
    }

    #[test]
    fn a_frame_without_fits_summarizes_to_zero_medians() {
        let (w, h) = (400, 300);
        let mut noise = vec![1000.0f32; w * h];
        add_noise(&mut noise, 10.0, 3);
        let donut: Vec<f32> = crate::test_support::donut_field(400, 200, 5)
            .iter()
            .map(|v| v * 65535.0)
            .collect();
        let o = FrameOptions::default();
        for (planes, w, h) in [(noise, w, h), (donut, 400, 200)] {
            let frame = RawFrame {
                planes: RawPlanes::Mono(planes),
                w,
                h,
                domain: RAW,
                flip_vertical: false,
                channels: 1,
            };
            let pm = measure_raw_frame(&frame, &o, None);
            let a = FrameSummary::from_measurement(&pm, &o);
            all_finite(&a);
            assert!(pm.stars.is_empty());
            assert_eq!(
                (
                    a.median_fwhm,
                    a.median_eccentricity,
                    a.median_snr,
                    a.median_hfr
                ),
                (0.0, 0.0, 0.0, 0.0)
            );
        }
    }

    #[test]
    fn a_trailed_frame_is_summarized_with_no_measured_star() {
        let p: Vec<f32> = crate::test_support::trailed_field(600, 450, 9)
            .iter()
            .map(|v| v * 65535.0)
            .collect();
        let frame = RawFrame {
            planes: RawPlanes::Mono(p),
            w: 600,
            h: 450,
            domain: RAW,
            flip_vertical: false,
            channels: 1,
        };
        let o = FrameOptions::default();
        let pm = measure_raw_frame(&frame, &o, None);
        let a = FrameSummary::from_measurement(&pm, &o);
        all_finite(&a);
        assert!(pm.trailed);
        assert!(pm.stars.is_empty());
        assert_eq!(
            (
                a.median_fwhm,
                a.median_eccentricity,
                a.median_snr,
                a.median_hfr
            ),
            (0.0, 0.0, 0.0, 0.0)
        );
        assert_eq!(a.trail_r_squared, 0.0);
    }

    /// The refusal reasons' sum is `stars_refused`; with `max_stars` below the seed count the budget leaves seeds
    /// `not_selected` (and may trim), and neither leaks into it.
    #[test]
    fn stars_refused_is_the_refusal_reasons_alone() {
        let scaled = |p: Vec<f32>| p.iter().map(|v| v * 65535.0).collect::<Vec<f32>>();
        let noise = {
            let mut n = vec![1000.0f32; 400 * 300];
            add_noise(&mut n, 10.0, 3);
            mono(n, 400, 300)
        };
        let trailed = mono(
            scaled(crate::test_support::trailed_field(600, 450, 9)),
            600,
            450,
        );
        let clean = mono_frame(120, 640, 480, 1000.0, 11);
        let o = FrameOptions::default();
        for frame in [&noise, &trailed, &clean] {
            let pm = measure_raw_frame(frame, &o, None);
            let a = FrameSummary::from_measurement(&pm, &o);
            // The sum of the refusal reasons; the budget's not-selected / trimmed are not in it.
            assert_eq!(
                a.stars_refused,
                pm.refused.by_reason.values().sum::<usize>()
            );
        }
        let a = summary(&trailed, &o);
        assert!(
            a.stars_refused > 0,
            "a trailed frame's fits are refused: {:?}",
            a.stars_refused
        );

        let budgeted = FrameOptions { max_stars: 20, ..o };
        let pm = measure_raw_frame(&clean, &budgeted, None);
        assert!(
            pm.refused.not_selected > 0,
            "premise: the budget left seeds unselected: {:?}",
            pm.refused
        );
        let a = FrameSummary::from_measurement(&pm, &budgeted);
        assert_eq!(a.stars_refused, pm.refused.total());
        assert_eq!(
            a.stars_refused,
            pm.refused.by_reason.values().sum::<usize>()
        );
    }

    #[test]
    fn the_summary_records_the_frame_verdict() {
        let (w, h) = (400, 300);
        let mut n = vec![1000.0f32; w * h];
        add_noise(&mut n, 10.0, 3);
        let donut: Vec<f32> = crate::test_support::donut_field(400, 200, 5)
            .iter()
            .map(|v| v * 65535.0)
            .collect();
        let trailed: Vec<f32> = crate::test_support::trailed_field(600, 450, 9)
            .iter()
            .map(|v| v * 65535.0)
            .collect();
        let o = FrameOptions::default();
        let verdict = |f: RawFrame| {
            let pm = measure_raw_frame(&f, &o, None);
            (FrameSummary::from_measurement(&pm, &o).verdict, pm.trailed)
        };
        assert_eq!(
            verdict(mono(trailed, 600, 450)),
            (Some(Verdict::Trailed), true)
        );
        assert_eq!(verdict(mono(donut, 400, 200)).0, Some(Verdict::Defocused));
        assert_eq!(verdict(mono(n, w, h)).0, Some(Verdict::NoStars));
        assert_eq!(verdict(mono_frame(120, 640, 480, 1000.0, 11)).0, None);
    }

    #[test]
    fn moffat4_summarizes_no_median_beta_and_auto_the_fitted_one() {
        let frame = mono_frame(120, 640, 480, 1000.0, 11);
        let m4 = summary(&frame, &FrameOptions::default());
        assert_eq!(m4.median_beta, None);
        let auto = FrameOptions {
            psf_model: PsfModel::Auto,
            ..FrameOptions::default()
        };
        let a = summary(&frame, &auto);
        let b = a.median_beta.expect("Auto summarizes the fitted beta");
        assert!(b.is_finite() && (b - 4.0).abs() < 1.0, "beta {b}");
    }

    fn site(p: BayerPattern, x: usize, y: usize) -> char {
        let (ex, ey) = (x % 2 == 0, y % 2 == 0);
        match (p, ex, ey) {
            (BayerPattern::Rggb, true, true) => 'R',
            (BayerPattern::Rggb, false, false) => 'B',
            _ => 'G',
        }
    }

    #[test]
    fn an_osc_mosaic_is_measured_on_green_and_a_three_plane_file_on_g() {
        // The mosaic: scene stars (flux 1.0 on G) x site gains + sky per site + noise.
        let (w, h) = (256, 256);
        let truth = (100.3, 80.7);
        let mut stars = vec![star(truth.0, truth.1, 3000.0, 4.0)];
        let mut rng = crate::geometry::ransac::SplitMix64(11);
        while stars.len() < 41 {
            let x = 16.0 + rng.next_f64() * (w as f64 - 32.0);
            let y = 16.0 + rng.next_f64() * (h as f64 - 32.0);
            if stars.iter().any(|s| (s.x - x).hypot(s.y - y) < 20.0) {
                continue;
            }
            stars.push(star(
                x,
                y,
                250.0 + rng.next_f64() * 500.0,
                3.6 + rng.next_f64() * 0.8,
            ));
        }
        let scene = moffat_field(w, h, &stars, 4.0, 0.0);
        let mut m: Vec<f32> = (0..w * h)
            .map(|i| {
                let (gain, sky) = match site(BayerPattern::Rggb, i % w, i / w) {
                    'R' => (0.5, 900.0),
                    'B' => (0.4, 800.0),
                    _ => (1.0, 1000.0),
                };
                scene[i] * gain + sky
            })
            .collect();
        add_noise(&mut m, 6.0, 5);
        let frame = RawFrame {
            planes: RawPlanes::Cfa {
                mosaic: m,
                pattern: BayerPattern::Rggb,
            },
            w,
            h,
            domain: RAW,
            flip_vertical: false,
            channels: 1,
        };
        let o = FrameOptions::default();
        let pm = measure_raw_frame(&frame, &o, None);
        let a = FrameSummary::from_measurement(&pm, &o);
        all_finite(&a);
        let best = pm
            .stars
            .iter()
            .min_by(|p, q| {
                (p.x - truth.0)
                    .hypot(p.y - truth.1)
                    .total_cmp(&(q.x - truth.0).hypot(q.y - truth.1))
            })
            .expect("a measured star");
        let d = (best.x - truth.0).hypot(best.y - truth.1);
        assert!(d < 0.1, "bright star {d} px off");
        assert!(a.stars_detected >= pm.stars.len() && a.stars_detected >= 35);

        // Three planes: 40 common stars, one only in G, one only in R.
        let (w, h) = (320, 240);
        let common = grid_stars(40, w, h, 3.5, 21);
        // Two free spots, 25 px from every common star and from each other.
        let mut free: Vec<(f64, f64)> = Vec::new();
        let mut rng = crate::geometry::ransac::SplitMix64(77);
        while free.len() < 2 {
            let (x, y) = (20.0 + rng.next_f64() * 280.0, 20.0 + rng.next_f64() * 200.0);
            if common
                .iter()
                .map(|s| (s.x, s.y))
                .chain(free.iter().copied())
                .all(|(a, b)| (a - x).hypot(b - y) > 25.0)
            {
                free.push((x, y));
            }
        }
        let g_only = star(free[0].0, free[0].1, 1500.0, 3.5);
        let r_only = star(free[1].0, free[1].1, 1500.0, 3.5);
        let plane = |extra: Option<MoffatStar>, seed| {
            let mut s = common.clone();
            s.extend(extra);
            let mut p = moffat_field(w, h, &s, 4.0, 1000.0);
            add_noise(&mut p, 10.0, seed);
            p
        };
        let frame = RawFrame {
            planes: RawPlanes::Rgb([
                plane(Some(r_only), 1),
                plane(Some(g_only), 2),
                plane(None, 3),
            ]),
            w,
            h,
            domain: RAW,
            flip_vertical: false,
            channels: 3,
        };
        let pm = measure_raw_frame(&frame, &o, None);
        let near = |s: &MoffatStar| pm.stars.iter().any(|m| (m.x - s.x).hypot(m.y - s.y) < 0.5);
        assert!(near(&g_only), "the G star is measured");
        assert!(!near(&r_only), "the R star is absent");
    }

    // ---- The list entry (forced-list zero point, R-OSP-19) ----

    /// `n` isolated stars (`(x, y, amplitude)`) at seeded positions at least `sep` px apart and `margin` from the edge.
    fn scattered(
        n: usize,
        w: usize,
        h: usize,
        sep: f64,
        margin: f64,
        amp: (f64, f64),
        seed: u64,
    ) -> Vec<(f64, f64, f64)> {
        let mut rng = crate::geometry::ransac::SplitMix64(seed);
        let mut out: Vec<(f64, f64, f64)> = Vec::new();
        let mut tries = 0;
        while out.len() < n {
            tries += 1;
            assert!(tries < 200_000, "{n} stars do not fit");
            let x = margin + rng.next_f64() * (w as f64 - 2.0 * margin);
            let y = margin + rng.next_f64() * (h as f64 - 2.0 * margin);
            if out.iter().any(|s| (s.0 - x).hypot(s.1 - y) < sep) {
                continue;
            }
            out.push((x, y, amp.0 + rng.next_f64() * (amp.1 - amp.0)));
        }
        out
    }

    fn integrated(w: usize, h: usize, stars: &[(f64, f64, f64)], fwhm: f64) -> Vec<f32> {
        use crate::test_support::{integrated_field, IntegratedProfile};
        integrated_field(
            w,
            h,
            stars,
            IntegratedProfile::Moffat { fwhm, beta: 2.5 },
            5,
            0.0,
        )
    }

    /// An RGGB (or GBRG) mosaic of the scene: R / B gains 0.5 / 0.4, greens `g` on even rows and odd rows, sky 1000 per
    /// site (+ `ped2` on the odd-row greens), noise 4.
    fn mosaic_of(
        scene: &[f32],
        w: usize,
        h: usize,
        pattern: BayerPattern,
        g: (f64, f64),
        ped2: f64,
        seed: u64,
    ) -> Vec<f32> {
        let mut m: Vec<f32> = (0..w * h)
            .map(|i| {
                let (x, y) = (i % w, i / w);
                let green_row = |row: usize| {
                    if row % 2 == 0 {
                        (g.0, 1000.0)
                    } else {
                        (g.1, 1000.0 + ped2)
                    }
                };
                let (gain, sky) = match pattern {
                    BayerPattern::Rggb => match (x % 2, y % 2) {
                        (0, 0) => (0.5, 900.0),
                        (1, 1) => (0.4, 800.0),
                        _ => green_row(y),
                    },
                    _ => match (x % 2, y % 2) {
                        (1, 0) => (0.5, 900.0),
                        (0, 1) => (0.4, 800.0),
                        _ => green_row(y),
                    },
                };
                (f64::from(scene[i]) * gain + sky) as f32
            })
            .collect();
        add_noise(&mut m, 4.0, seed);
        m
    }

    fn raw(planes: RawPlanes, w: usize, h: usize, channels: usize) -> RawFrame {
        RawFrame {
            planes,
            w,
            h,
            domain: RAW,
            flip_vertical: false,
            channels,
        }
    }

    /// Every number of a fit as bits (NaN equals NaN by its bits).
    fn fit_bits(f: &StarFit) -> Vec<u64> {
        let mut v: Vec<u64> = [
            f.x,
            f.y,
            f.background,
            f.amplitude,
            f.fwhm_x,
            f.fwhm_y,
            f.fwtm_x,
            f.fwtm_y,
            f.theta,
            f.beta,
            f.residual,
            f.signal,
            f.area,
            f.peak5,
            f.sky_ring,
            f.sky_ring_sigma,
            f.snr,
            f.ap_sum,
            f.hfr,
        ]
        .iter()
        .map(|x| x.to_bits())
        .collect();
        v.push(u64::from(f.sky_ring_n));
        v.push(u64::from(f.flags.bits()));
        v
    }

    /// The measurement with its wall-clock fields reset, as text: `{:?}` prints a float as its shortest round-trip, so
    /// two equal strings are two equal measurements field by field.
    fn measurement_text(pm: &PlaneMeasurement) -> String {
        let mut pm = pm.clone();
        pm.timings = Default::default();
        pm.sub_steps = Default::default();
        format!("{pm:?}")
    }

    /// B9's rustafits half: the Analysis measurement is bit-identical with and without a list, on a mono frame, an RGGB
    /// and a GBRG mosaic and a 3-plane file; the list is 200 positions, which are fitted (the 3-plane file measures no
    /// list, R4c, and its measurement is the frame's own all the same).
    #[test]
    fn the_analysis_measurement_is_identical_with_and_without_a_list() {
        let (w, h) = (400, 400);
        let stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
        let scene = integrated(w, h, &stars, 3.0);
        let mut mono_plane: Vec<f32> = scene.iter().map(|&v| v + 1000.0).collect();
        add_noise(&mut mono_plane, 4.0, 6);
        let plane = |seed| {
            let mut p: Vec<f32> = scene.iter().map(|&v| v + 1000.0).collect();
            add_noise(&mut p, 4.0, seed);
            p
        };
        let frames = vec![
            ("mono", raw(RawPlanes::Mono(mono_plane), w, h, 1)),
            (
                "rggb",
                raw(
                    RawPlanes::Cfa {
                        mosaic: mosaic_of(&scene, w, h, BayerPattern::Rggb, (1.0, 0.98), 3.0, 7),
                        pattern: BayerPattern::Rggb,
                    },
                    w,
                    h,
                    1,
                ),
            ),
            (
                "gbrg",
                raw(
                    RawPlanes::Cfa {
                        mosaic: mosaic_of(&scene, w, h, BayerPattern::Gbrg, (1.0, 0.98), 3.0, 8),
                        pattern: BayerPattern::Gbrg,
                    },
                    w,
                    h,
                    1,
                ),
            ),
            (
                "rgb",
                raw(RawPlanes::Rgb([plane(1), plane(2), plane(3)]), w, h, 3),
            ),
        ];
        let o = FrameOptions::default();
        let mut rng = crate::geometry::ransac::SplitMix64(99);
        let positions: Vec<(f64, f64)> = (0..200)
            .map(|i| {
                if i < stars.len() {
                    (stars[i].0 + 0.1, stars[i].1 - 0.1)
                } else {
                    (20.0 + rng.next_f64() * 360.0, 20.0 + rng.next_f64() * 360.0)
                }
            })
            .collect();
        for (name, frame) in &frames {
            let plain = measure_raw_frame(frame, &o, None);
            let (with, list) = measure_raw_frame_and_list(frame, &o, None, |_| {
                Some(ListRequest::new(positions.clone(), 1.0))
            });
            let list = list.unwrap_or_else(|| panic!("{name}: the list was not measured"));
            if *name == "rgb" {
                assert_eq!(list.refusal, Some(ListRefusal::ColourFrame));
                assert!(list.fits.is_empty());
            } else {
                assert_eq!(list.fits.len(), 200, "{name}");
                assert!(
                    list.fits.iter().filter(|f| f.accepted).count() >= 20,
                    "{name}: the list is fitted"
                );
            }
            assert_eq!(plain.stars.len(), with.stars.len(), "{name}");
            for (a, b) in plain.stars.iter().zip(&with.stars) {
                assert_eq!(fit_bits(a), fit_bits(b), "{name}: a star moved");
            }
            assert_eq!(measurement_text(&plain), measurement_text(&with), "{name}");
            let (sa, sb) = (
                format!("{:?}", FrameSummary::from_measurement(&plain, &o)),
                format!("{:?}", FrameSummary::from_measurement(&with, &o)),
            );
            assert_eq!(sa, sb, "{name}: the summary moved");
        }
    }

    /// The `test-support` constructors build the list rows a caller's own pure tests need.
    #[test]
    fn the_test_constructors_build_list_rows() {
        let mut off = ListFit::test_kept(10.0, 20.0, 500.0);
        off.kept = false;
        let list = ListMeasurement::test_new(
            vec![
                ListFit::test_kept(1.0, 2.0, 100.0),
                off,
                ListFit::test_refused(),
            ],
            Some(1.2),
            7,
        );
        assert_eq!(list.fits.len(), 3);
        assert!(list.fits[0].accepted && list.fits[0].kept && list.fits[0].flux == 100.0);
        assert!(list.fits[1].accepted && !list.fits[1].kept);
        assert!(!list.fits[2].accepted && list.fits[2].flux.is_nan());
        assert_eq!(
            (list.aperture_correction, list.aperture_stars, list.refusal),
            (Some(1.2), 7, None)
        );
    }

    /// A list's `flux` over the true total flux of each star: `(median ratio, median centre offset)`.
    fn recovery(
        frame: &RawFrame,
        stars: &[(f64, f64, f64)],
        fwhm: f64,
        truth_gain: f64,
    ) -> (f64, f64, ListMeasurement) {
        let positions: Vec<(f64, f64)> = stars.iter().map(|s| (s.0 + 0.1, s.1 - 0.1)).collect();
        let (_, list) = measure_raw_frame_and_list(frame, &FrameOptions::default(), None, |_| {
            Some(ListRequest::new(positions.clone(), 1.0))
        });
        let list = list.expect("the list is measured");
        let alpha = fwhm / (2.0 * (2f64.powf(0.4) - 1.0).sqrt());
        let mut ratios = Vec::new();
        let mut offsets = Vec::new();
        for (lf, s) in list.fits.iter().zip(stars) {
            if lf.kept {
                let total = std::f64::consts::PI * alpha * alpha * s.2 / 1.5 * truth_gain;
                ratios.push(lf.flux / total);
                offsets.push((lf.x - s.0).hypot(lf.y - s.1));
            }
        }
        assert!(
            ratios.len() >= 20,
            "kept {} of {}: {:?}",
            ratios.len(),
            stars.len(),
            list.aperture_stars
        );
        (
            crate::star_measure::tests::median(ratios.iter().copied()),
            crate::star_measure::tests::median(offsets.iter().copied()),
            list,
        )
    }

    #[test]
    fn a_list_on_a_mono_field_recovers_its_fluxes() {
        let (w, h) = (400, 400);
        for fwhm in [2.0, 3.0] {
            let stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
            let mut p: Vec<f32> = integrated(w, h, &stars, fwhm)
                .iter()
                .map(|&v| v + 1000.0)
                .collect();
            add_noise(&mut p, 4.0, 6);
            let (ratio, offset, list) =
                recovery(&raw(RawPlanes::Mono(p), w, h, 1), &stars, fwhm, 1.0);
            assert!(
                (ratio - 1.0).abs() < 0.01,
                "fwhm {fwhm}: median flux / truth {ratio}"
            );
            assert!(offset < 0.2, "fwhm {fwhm}: median centre offset {offset}");
            assert!(
                list.aperture_correction.is_some_and(|c| c > 1.0 && c < 1.5),
                "{:?}",
                list.aperture_correction
            );
        }
    }

    #[test]
    fn a_list_on_a_mosaic_recovers_its_fluxes() {
        let (w, h) = (400, 400);
        for fwhm in [2.0, 3.0] {
            let stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
            let scene = integrated(w, h, &stars, fwhm);
            let m = mosaic_of(&scene, w, h, BayerPattern::Rggb, (1.00, 0.98), 3.0, 7);
            let frame = raw(
                RawPlanes::Cfa {
                    mosaic: m,
                    pattern: BayerPattern::Rggb,
                },
                w,
                h,
                1,
            );
            // The two greens carry gains 1.00 and 0.98: the fit and the growth see their mean.
            let (ratio, offset, _) = recovery(&frame, &stars, fwhm, 0.99);
            assert!(
                (ratio - 1.0).abs() < 0.02,
                "fwhm {fwhm}: median flux / truth {ratio}"
            );
            assert!(offset < 0.2, "fwhm {fwhm}: median centre offset {offset}");
        }
    }

    #[test]
    fn a_list_star_whose_own_detection_coincides_still_qualifies_for_c() {
        let (w, h) = (400, 400);
        let stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
        let mut p: Vec<f32> = integrated(w, h, &stars, 3.0)
            .iter()
            .map(|&v| v + 1000.0)
            .collect();
        add_noise(&mut p, 4.0, 6);
        let frame = raw(RawPlanes::Mono(p), w, h, 1);
        let positions: Vec<(f64, f64)> = stars.iter().map(|s| (s.0, s.1)).collect();
        let (pm, list) = measure_raw_frame_and_list(&frame, &FrameOptions::default(), None, |_| {
            Some(ListRequest::new(positions.clone(), 1.0))
        });
        let detected_at =
            |x: f64, y: f64| pm.detected.iter().any(|d| (d.x - x).hypot(d.y - y) < 0.5);
        assert!(
            positions.iter().all(|&(x, y)| detected_at(x, y)),
            "premise: every list star is also detected"
        );
        let list = list.expect("measured");
        assert!(list.aperture_stars >= 5, "{}", list.aperture_stars);
        assert!(list.aperture_correction.is_some());
    }

    #[test]
    fn the_c_cap_is_the_requests() {
        let (w, h) = (1100, 900);
        let mut stars = Vec::new();
        for j in 0..18 {
            for i in 0..22 {
                stars.push((
                    40.0 + 48.0 * i as f64 + 6.0 * ((i * 7 + j * 3) % 5) as f64,
                    40.0 + 48.0 * j as f64 + 5.0 * ((i * 3 + j * 11) % 5) as f64,
                    4000.0 + 30.0 * ((i + 3 * j) % 17) as f64,
                ));
            }
        }
        assert!(stars.len() >= 300);
        let mut p: Vec<f32> = integrated(w, h, &stars, 2.5)
            .iter()
            .map(|&v| v + 1000.0)
            .collect();
        add_noise(&mut p, 4.0, 6);
        let frame = raw(RawPlanes::Mono(p), w, h, 1);
        let positions: Vec<(f64, f64)> = stars.iter().map(|s| (s.0, s.1)).collect();
        let (pm, list) = measure_raw_frame_and_list(&frame, &FrameOptions::default(), None, |_| {
            Some(ListRequest::new(positions.clone(), 1.0))
        });
        let list = list.expect("measured");
        assert_eq!(list.aperture_stars, 256, "the request's cap");
        assert!(
            pm.aperture_stars <= 64,
            "Analysis' own cap: {}",
            pm.aperture_stars
        );
        let capped = measure_raw_frame_and_list(&frame, &FrameOptions::default(), None, |_| {
            let mut r = ListRequest::new(positions.clone(), 1.0);
            r.max_c_stars = 40;
            Some(r)
        })
        .1
        .expect("measured");
        assert_eq!(capped.aperture_stars, 40);
    }

    /// A star's own flux is its model's times the list's `c`, and `c` reads each qualifier's ring sky, so the proof that
    /// an unlisted neighbour does not lift a ring sky is on `c`: every one of eight list stars has an unlisted companion
    /// 11 F away (inside its ring), and the list's fluxes match the field without the companions. A companion as bright
    /// as the list's stars moves `c` by under 0.5 %; one ten times brighter by at most 2 %, the physical wing of a
    /// star ten times brighter outside the masked 3 F disc (1.4 % measured), not a defect of the mask (M3 ruling).
    #[test]
    fn an_unlisted_neighbour_does_not_lift_the_ring_sky() {
        unlisted_neighbour(3000.0, 0.005);
    }

    #[test]
    fn an_unlisted_neighbour_ten_times_brighter_moves_c_by_its_wings_only() {
        unlisted_neighbour(30000.0, 0.02);
    }

    fn unlisted_neighbour(companion_amp: f64, bar: f64) {
        let (w, h) = (500, 500);
        let fwhm = 3.0;
        let stars = scattered(8, w, h, 100.0, 70.0, (3000.0, 5000.0), 5);
        let companion = |t: &(f64, f64, f64)| {
            (0..16)
                .map(|k| {
                    let a = k as f64 * std::f64::consts::TAU / 16.0;
                    (
                        t.0 + 11.0 * fwhm * a.cos(),
                        t.1 + 11.0 * fwhm * a.sin(),
                        companion_amp,
                    )
                })
                .find(|b| {
                    b.0 > 40.0
                        && b.0 < 460.0
                        && b.1 > 40.0
                        && b.1 < 460.0
                        && stars
                            .iter()
                            .filter(|s| *s != t)
                            .all(|s| (s.0 - b.0).hypot(s.1 - b.1) > 40.0)
                })
                .expect("a free spot")
        };
        let companions: Vec<_> = stars.iter().map(companion).collect();
        let frame_of = |extra: &[(f64, f64, f64)]| {
            let mut all = stars.clone();
            all.extend_from_slice(extra);
            let mut p: Vec<f32> = integrated(w, h, &all, fwhm)
                .iter()
                .map(|&v| v + 1000.0)
                .collect();
            add_noise(&mut p, 4.0, 6);
            raw(RawPlanes::Mono(p), w, h, 1)
        };
        let positions: Vec<(f64, f64)> = stars.iter().map(|s| (s.0, s.1)).collect();
        let run = |frame: &RawFrame| {
            let (pm, list) =
                measure_raw_frame_and_list(frame, &FrameOptions::default(), None, |_| {
                    Some(ListRequest::new(positions.clone(), 1.0))
                });
            (pm, list.expect("measured"))
        };
        let (_, without) = run(&frame_of(&[]));
        let (pm, with) = run(&frame_of(&companions));
        for b in &companions {
            assert!(
                pm.detected
                    .iter()
                    .any(|d| (d.x - b.0).hypot(d.y - b.1) < 1.0),
                "premise: the unlisted companion {b:?} is detected"
            );
        }
        let (ca, cb) = (
            without.aperture_correction.expect("c"),
            with.aperture_correction.expect("c"),
        );
        assert!(without.aperture_stars >= 5 && with.aperture_stars >= 5);
        assert!(
            (cb / ca - 1.0).abs() < bar,
            "c {cb} against {ca} without the companions"
        );
        // (The field's brightest star can read as clipped by the plane-max amplitude clause, and a brighter companion
        // lifts the plane max and clears it: only the stars kept in both runs are compared.)
        let mut compared = 0;
        for (a, b) in without.fits.iter().zip(&with.fits) {
            if a.kept && b.kept {
                compared += 1;
                assert!(
                    (b.flux / a.flux - 1.0).abs() < bar,
                    "flux {} against {}",
                    b.flux,
                    a.flux
                );
            }
        }
        assert!(compared >= 6, "{compared} kept");
    }

    #[test]
    fn a_saturated_list_star_is_not_kept() {
        let (w, h) = (400, 400);
        let mut stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
        stars[0].2 = 90000.0;
        let mut p: Vec<f32> = integrated(w, h, &stars, 3.0)
            .iter()
            .map(|&v| (v + 1000.0).min(65535.0))
            .collect();
        add_noise(&mut p, 4.0, 6);
        p.iter_mut().for_each(|v| *v = v.min(65535.0));
        let positions: Vec<(f64, f64)> = stars.iter().map(|s| (s.0, s.1)).collect();
        let (_, list) = measure_raw_frame_and_list(
            &raw(RawPlanes::Mono(p), w, h, 1),
            &FrameOptions::default(),
            None,
            |_| Some(ListRequest::new(positions.clone(), 1.0)),
        );
        let list = list.expect("measured");
        let sat = list.fits[0];
        assert!(
            sat.accepted && sat.flags.saturated(),
            "the clipped star is accepted and flagged saturated: {sat:?}"
        );
        assert!(!sat.kept, "{sat:?}");
        assert!(list.fits[1..].iter().filter(|f| f.kept).count() >= 20);
    }

    #[test]
    fn a_far_centre_is_not_kept() {
        let (w, h) = (400, 400);
        let stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
        let mut p: Vec<f32> = integrated(w, h, &stars, 3.0)
            .iter()
            .map(|&v| v + 1000.0)
            .collect();
        add_noise(&mut p, 4.0, 6);
        let mut positions: Vec<(f64, f64)> = stars.iter().map(|s| (s.0, s.1)).collect();
        positions[0].0 += 3.0;
        let (_, list) = measure_raw_frame_and_list(
            &raw(RawPlanes::Mono(p), w, h, 1),
            &FrameOptions::default(),
            None,
            |_| Some(ListRequest::new(positions.clone(), 1.0)),
        );
        let list = list.expect("measured");
        let far = list.fits[0];
        assert!(!far.kept, "{far:?}");
        assert!(!far.accepted || far.centre_offset_px > 1.0, "{far:?}");
        assert!(list.fits[1..].iter().filter(|f| f.kept).count() >= 20);
    }

    /// A trailed or defocused plane measures no list and does not call the closure; the plane's verdict is the reason a
    /// caller stores (R5: `frame_verdict`).
    #[test]
    fn a_trailed_or_defocused_plane_measures_no_list_and_carries_its_verdict() {
        let trailed: Vec<f32> = crate::test_support::trailed_field(600, 450, 9)
            .iter()
            .map(|v| v * 65535.0)
            .collect();
        let donuts: Vec<f32> = crate::test_support::donut_field(330, 100, 5)
            .iter()
            .map(|v| v * 65535.0)
            .collect();
        for (name, frame, want) in [
            ("trailed", mono(trailed, 600, 450), Verdict::Trailed),
            ("defocused", mono(donuts, 330, 100), Verdict::Defocused),
        ] {
            let asked = std::cell::Cell::new(false);
            let o = FrameOptions::default();
            let (pm, list) = measure_raw_frame_and_list(&frame, &o, None, |_| {
                asked.set(true);
                Some(ListRequest::new(vec![(100.0, 50.0)], 1.0))
            });
            assert!(list.is_none(), "{name}");
            assert!(
                !asked.get(),
                "{name}: a plane with a verdict does not build a list"
            );
            assert_eq!(
                FrameSummary::from_measurement(&pm, &o).verdict,
                Some(want),
                "{name}"
            );
            assert_eq!(pm.trailed, want == Verdict::Trailed, "{name}");
            assert_eq!(pm.defocused, want == Verdict::Defocused, "{name}");
        }
    }

    /// R4c: a debayered three-plane file keeps its measurement and gets an empty list with the layout's refusal; the
    /// closure is not called.
    #[test]
    fn a_three_plane_file_measures_no_list_and_says_why() {
        let (w, h) = (400, 400);
        let stars = scattered(20, w, h, 50.0, 40.0, (1500.0, 5000.0), 5);
        let scene = integrated(w, h, &stars, 3.0);
        let plane = |seed| {
            let mut p: Vec<f32> = scene.iter().map(|&v| v + 1000.0).collect();
            add_noise(&mut p, 4.0, seed);
            p
        };
        let frame = raw(RawPlanes::Rgb([plane(1), plane(2), plane(3)]), w, h, 3);
        let o = FrameOptions::default();
        let asked = std::cell::Cell::new(false);
        let (pm, list) = measure_raw_frame_and_list(&frame, &o, None, |_| {
            asked.set(true);
            Some(ListRequest::new(vec![(stars[0].0, stars[0].1)], 1.0))
        });
        let list = list.expect("a refusal the caller can store");
        assert!(!asked.get());
        assert!(
            list.fits.is_empty() && list.aperture_correction.is_none() && list.aperture_stars == 0
        );
        assert_eq!(list.refusal, Some(ListRefusal::ColourFrame));
        assert_eq!(list.refusal.map(ListRefusal::code), Some("colour_frame"));
        assert_eq!(
            measurement_text(&pm),
            measurement_text(&measure_raw_frame(&frame, &o, None))
        );
        assert!(!pm.stars.is_empty());
    }

    /// I2: a position 1.25 px off its star (inside the fitter's free-centre bound, outside the kept tolerance) is
    /// accepted and not kept, and `c` is over the kept fits alone.
    #[test]
    fn a_fit_centred_off_its_position_is_accepted_but_not_kept_and_not_in_c() {
        let (w, h) = (400, 400);
        let stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
        let mut p: Vec<f32> = integrated(w, h, &stars, 3.0)
            .iter()
            .map(|&v| v + 1000.0)
            .collect();
        add_noise(&mut p, 4.0, 6);
        let frame = raw(RawPlanes::Mono(p), w, h, 1);
        let mut positions: Vec<(f64, f64)> = stars.iter().map(|s| (s.0, s.1)).collect();
        positions[0].0 += 1.25;
        let (_, list) = measure_raw_frame_and_list(&frame, &FrameOptions::default(), None, |_| {
            Some(ListRequest::new(positions.clone(), 1.0))
        });
        let list = list.expect("measured");
        let off = list.fits[0];
        assert!(off.accepted && !off.kept, "{off:?}");
        assert!(
            off.centre_offset_px > 1.0 && off.centre_offset_px < 1.5,
            "{off:?}"
        );
        let kept = list.fits.iter().filter(|f| f.kept).count();
        assert!(kept >= 20);
        assert!(
            list.aperture_stars <= kept,
            "c is over {} stars, {kept} are kept",
            list.aperture_stars
        );
    }

    /// M2: a position off the plane, not finite or on the far edge is refused with NaN numbers and no panic, and the
    /// other positions are fitted as without it.
    #[test]
    fn a_position_off_the_plane_or_not_finite_is_refused_and_the_rest_unchanged() {
        let (w, h) = (400, 400);
        let stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
        let mut p: Vec<f32> = integrated(w, h, &stars, 3.0)
            .iter()
            .map(|&v| v + 1000.0)
            .collect();
        add_noise(&mut p, 4.0, 6);
        let frame = raw(RawPlanes::Mono(p), w, h, 1);
        let good: Vec<(f64, f64)> = stars.iter().map(|s| (s.0, s.1)).collect();
        let bad = [
            (-5.0, 10.0),
            (f64::NAN, 20.0),
            (30.0, f64::INFINITY),
            (w as f64 + 3.0, 5.0),
            (5.0, h as f64 - 0.5),
        ];
        let run = |positions: Vec<(f64, f64)>| {
            measure_raw_frame_and_list(&frame, &FrameOptions::default(), None, |_| {
                Some(ListRequest::new(positions, 1.0))
            })
            .1
            .expect("measured")
        };
        let base = run(good.clone());
        let mut mixed: Vec<(f64, f64)> = bad.to_vec();
        mixed.extend_from_slice(&good);
        let with = run(mixed);
        assert_eq!(with.fits.len(), bad.len() + good.len());
        for (i, f) in with.fits[..bad.len()].iter().enumerate() {
            assert!(!f.accepted && !f.kept, "{i}: {f:?}");
            assert!(
                f.x.is_nan() && f.y.is_nan() && f.flux.is_nan() && f.model_flux.is_nan(),
                "{i}: {f:?}"
            );
        }
        assert_eq!(with.aperture_stars, base.aperture_stars);
        assert_eq!(
            with.aperture_correction.map(f64::to_bits),
            base.aperture_correction.map(f64::to_bits)
        );
        for (a, b) in base.fits.iter().zip(&with.fits[bad.len()..]) {
            assert_eq!(
                (a.kept, a.x.to_bits(), a.flux.to_bits()),
                (b.kept, b.x.to_bits(), b.flux.to_bits())
            );
        }
    }

    /// R3: the near-clip class is the MEASURED plane's floor. A clip-free plane has none (its list derives none), so a
    /// floor on the measurement is the only way a list star can be near the clip: set one at half the range and every
    /// star with a top above it is flagged, where the unchanged measurement flags none of the stars below the brightest.
    #[test]
    fn the_list_reads_the_planes_clip_floor_not_one_from_its_own_seeds() {
        use crate::star_measure::list::{measure_list, ListPlane};
        let (w, h) = (400, 400);
        let stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
        let mut p: Vec<f32> = integrated(w, h, &stars, 3.0)
            .iter()
            .map(|&v| v + 1000.0)
            .collect();
        add_noise(&mut p, 4.0, 6);
        let o = FrameOptions::default();
        let frame = raw(RawPlanes::Mono(p), w, h, 1);
        let pm = measure_raw_frame(&frame, &o, None);
        assert!(pm.clip_floor.is_none(), "premise: the plane shows no clip");
        let RawPlanes::Mono(plane) = &frame.planes else {
            unreachable!()
        };
        let input = ListPlane {
            plane,
            w,
            h,
            domain: RAW,
            mask: None,
            mosaic: false,
        };
        let req = ListRequest::new(stars.iter().map(|s| (s.0, s.1)).collect(), 1.0);
        let run = |pm: &PlaneMeasurement| {
            measure_list(&input, pm, &pm.detected, &req, None).expect("measured")
        };
        let plain = run(&pm);
        let mut floored = pm.clone();
        floored.clip_floor = Some(0.5);
        floored.clip_evidence = 1;
        let with = run(&floored);
        let amp_max = stars.iter().map(|s| s.2).fold(0.0, f64::max);
        let mut flagged = 0;
        for ((a, b), s) in plain.fits.iter().zip(&with.fits).zip(&stars) {
            if s.2 > 0.75 * amp_max {
                assert!(
                    b.flags.saturated() && !b.kept,
                    "a star at {:.2} of the range: {b:?}",
                    s.2 / amp_max
                );
                flagged += usize::from(!a.flags.saturated());
            }
            if s.2 < 0.4 * amp_max {
                assert!(!b.flags.saturated(), "{b:?}");
            }
        }
        assert!(
            flagged >= 2,
            "the floor flagged {flagged} stars the plain list did not"
        );
    }

    /// R2: a mosaic's detection is a `2s + 0.5` centroid of the green-sum plane, up to about a pixel off the star's true
    /// position, which the list position is not. Every detection moved 1.3 px: each is still its list star's own (the
    /// nearest within max(tolerance, 1 F)), so no list star is a blend of its own detection and `c` is measured.
    #[test]
    fn a_mosaic_stars_own_detection_one_point_three_px_off_is_still_its_own() {
        use crate::star_measure::list::{measure_list, ListPlane};
        let (w, h) = (400, 400);
        let stars = scattered(30, w, h, 50.0, 45.0, (1500.0, 5000.0), 5);
        let scene = integrated(w, h, &stars, 3.0);
        let mosaic = mosaic_of(&scene, w, h, BayerPattern::Rggb, (1.0, 0.98), 3.0, 7);
        let frame = raw(
            RawPlanes::Cfa {
                mosaic: mosaic.clone(),
                pattern: BayerPattern::Rggb,
            },
            w,
            h,
            1,
        );
        let base = base_options(&frame, &FrameOptions::default(), false);
        let (pm, seeds) =
            measure_green_sites_and_seeds(&mosaic, w, h, BayerPattern::Rggb, &base, None);
        let mask = crate::processing::cfa::green_mask(w, h, BayerPattern::Rggb);
        let input = ListPlane {
            plane: &mosaic,
            w,
            h,
            domain: RAW,
            mask: Some(&mask),
            mosaic: true,
        };
        let req = ListRequest::new(stars.iter().map(|s| (s.0, s.1)).collect(), 1.0);
        let moved: Vec<Seed> = seeds.iter().map(|s| Seed { x: s.x + 1.3, ..*s }).collect();
        assert!(
            pm.field_width.fwhm_px > 1.3,
            "premise: F_field {}",
            pm.field_width.fwhm_px
        );
        let list = measure_list(&input, &pm, &moved, &req, None).expect("measured");
        let blends = list
            .refused
            .by_reason
            .get(&crate::star_measure::RefusalReason::Blend)
            .copied()
            .unwrap_or(0);
        assert_eq!(blends, 0, "{:?}", list.refused);
        assert!(list.fits.iter().filter(|f| f.kept).count() >= 20);
        assert!(list.aperture_stars >= 5, "c: {} stars", list.aperture_stars);
        assert!(list.aperture_correction.is_some());
        // Control: the list's own detections, unmoved.
        let plain = measure_list(&input, &pm, &seeds, &req, None).expect("measured");
        assert!(plain.aperture_stars >= 5);
    }

    #[test]
    fn an_empty_request_measures_no_list() {
        let frame = mono_frame(120, 640, 480, 1000.0, 11);
        let o = FrameOptions::default();
        let (_, a) = measure_raw_frame_and_list(&frame, &o, None, |_| {
            Some(ListRequest::new(Vec::new(), 1.0))
        });
        let (_, b) = measure_raw_frame_and_list(&frame, &o, None, |_| None);
        let (_, c) = measure_raw_frame_and_list(&frame, &o, None, |_| {
            Some(ListRequest::new(vec![(50.0, 50.0)], 0.0))
        });
        assert!(a.is_none() && b.is_none() && c.is_none());
    }
}
