//! The Analysis tab's measurement: a raw frame measured through the one star path (Plan C, Task C5).
//!
//! [`measure_raw_frame`] runs `star_measure` on the frame (mono, a Bayer mosaic's green sites, or a debayered file's G
//! plane); [`FrameSummary::from_measurement`] reduces the result to the frame figures the Analysis tab shows. Positions
//! stay in FILE row order (unflipped); [`super::raw::RawFrame::flip_vertical`] tells a display to flip them.
//!
//! This module and [`super::raw`] are the stable entry of `star_measure`: a caller reads a file with
//! [`super::raw::read_raw_frame`], measures it with [`measure_raw_frame`] and reads the figures off the
//! [`FrameSummary`]. The row mapping (the per-star rows, the frame row's ids, `config_hash`, `analyzed_at`,
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
use crate::star_measure::green::measure_green_sites;
use crate::star_measure::raw::{RawFrame, RawPlanes};
use crate::star_measure::trail::trail_r_squared;
use crate::star_measure::{
    self, BetaSource, MeasureOptions, Mode, PixelNoise, PlaneMeasurement, PsfModel, SeedSpec, StarFit,
    DEFAULT_THRESHOLD,
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
    let base = MeasureOptions {
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
    };
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
            let opts = MeasureOptions {
                colour: true,
                pixel_noise: PixelNoise::Mixed,
                cfa_site_share: Some(0.5),
                ..base
            };
            star_measure::measure_plane(&planes[1], frame.w, frame.h, &opts, pool)
        }
    }
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
/// `Auto` and `None` otherwise (PC-9). A non-finite figure is written 0.0, so no figure is NaN.
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
            assert!(
                s.fwhm_x >= s.fwhm_y && s.eccentricity() >= 0.0 && s.eccentricity() < 1.0
            );
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
        let g = gaussian_field(w, h, &[(knot.0, knot.1, 400.0)], 3.0 * 3.5 / 2.354_820_045, 0.0);
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
            assert!(pm.extended_seeds >= 1, "{name}: the knot is flagged extended");
            if name == "cfa" {
                assert!(pm.detected.is_empty(), "{name}: handed seeds leave `detected` empty");
            } else {
                assert_eq!(pm.detected.len(), pm.seeds, "{name}");
            }
            let a = summary(&plain, &o);
            let b = summary(&knotted, &o);
            assert_eq!(b.stars_detected, pm.seeds - pm.extended_seeds, "{name}");
            assert_eq!(b.stars_detected, a.stars_detected, "{name}: the knot is not counted");
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
                (a.median_fwhm, a.median_eccentricity, a.median_snr, a.median_hfr),
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
            (a.median_fwhm, a.median_eccentricity, a.median_snr, a.median_hfr),
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
        let trailed = mono(scaled(crate::test_support::trailed_field(600, 450, 9)), 600, 450);
        let clean = mono_frame(120, 640, 480, 1000.0, 11);
        let o = FrameOptions::default();
        for frame in [&noise, &trailed, &clean] {
            let pm = measure_raw_frame(frame, &o, None);
            let a = FrameSummary::from_measurement(&pm, &o);
            // The sum of the refusal reasons; the budget's not-selected / trimmed are not in it.
            assert_eq!(a.stars_refused, pm.refused.by_reason.values().sum::<usize>());
        }
        let a = summary(&trailed, &o);
        assert!(a.stars_refused > 0, "a trailed frame's fits are refused: {:?}", a.stars_refused);

        let budgeted = FrameOptions { max_stars: 20, ..o };
        let pm = measure_raw_frame(&clean, &budgeted, None);
        assert!(pm.refused.not_selected > 0, "premise: the budget left seeds unselected: {:?}", pm.refused);
        let a = FrameSummary::from_measurement(&pm, &budgeted);
        assert_eq!(a.stars_refused, pm.refused.total());
        assert_eq!(a.stars_refused, pm.refused.by_reason.values().sum::<usize>());
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
        assert_eq!(verdict(mono(trailed, 600, 450)), (Some(Verdict::Trailed), true));
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
            stars.push(star(x, y, 250.0 + rng.next_f64() * 500.0, 3.6 + rng.next_f64() * 0.8));
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
            planes: RawPlanes::Rgb([plane(Some(r_only), 1), plane(Some(g_only), 2), plane(None, 3)]),
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
}
