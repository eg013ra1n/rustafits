//! The CFA green-site path of the one star path (Plan C, Task C3, ruling PC-4).
//!
//! A Bayer mosaic is measured on its green sites. Detection runs on the 2x2 sum of the two green sites of each block
//! (G1 + G2, ruling R-PC-3: the fit reads green only, so R / B sky, chromatic focus and R / B hot pixels cannot move a
//! seed; the sum's noise is white); every seed is mapped to native pixels, and those the matched filter flagged
//! extended (knots, R-PC-8) are handed to the fit flagged: fitted, and kept out of the frame statistic as on the mono
//! path; the fit runs on the native mosaic with the green-site mask; the plane's sky and noise are read from one
//! green phase. `c` is not applied here (Analysis star fluxes on an OSC frame stay uncorrected, PC-4); since the OSC zero point it is
//! measured on the green sites and reported apart (`PlaneMeasurement::reported_aperture_correction`).
//!
//! Positions are in FILE row order (the pattern addresses that order, see `raw::RawPlanes`): the caller flips them
//! for display.
//!
//! **In rustafits** (port plan phase 4, `docs/superpowers/plans/2026-11-04-star-measure-port-phase-4.md`): copied from
//! Athenaeum's `analysis/cfa.rs` at e6344ddd, paths only. The rulings the docs cite (`PC-…`, `R-PC-…`) are Athenaeum's.

// Athenaeum's test helper `site` as copied (port plan G-3): the clippy style lint it trips is allowed, not "fixed".
#![allow(clippy::manual_is_multiple_of)]

use std::sync::Arc;

use crate::processing::cfa::{green_mask, green_phase, green_sum};
use crate::BayerPattern;
use tracing::error;

use crate::star_measure::FieldWidth;
use crate::star_measure::{
    self, domain, seeds, MeasureOptions, NoiseSource, PixelNoise, PlaneMeasurement, Seed, SeedSpec,
    DEFAULT_THRESHOLD,
};

/// Native-pixel offset of a super-pixel seed. The module's pixel-centre convention puts a pixel's centre at its integer
/// coordinate, so super-pixel `i` (native columns `2i`, `2i + 1`) is centred at native `2i + 0.5`: a seed at `s` on the
/// super-pixel plane is at `2s + SUPERPIXEL_CENTRE_OFFSET` natively. Pinned at seed level by
/// `the_super_pixel_offset_is_half_a_pixel_at_seed_level` (the fit's position test cannot: the fitted centre is free).
pub const SUPERPIXEL_CENTRE_OFFSET: f64 = 0.5;

/// A super-pixel-plane seed in native pixels: positions x2 plus the centre offset; the peak is a 2-site sum, so one
/// site's is a half (the flux stays: the fit re-measures it).
pub fn native_seed(s: &Seed) -> Seed {
    Seed {
        x: 2.0 * s.x + SUPERPIXEL_CENTRE_OFFSET,
        y: 2.0 * s.y + SUPERPIXEL_CENTRE_OFFSET,
        peak: s.peak / 2.0,
        flux: s.flux,
    }
}

/// Measure a Bayer mosaic on its green sites (PC-4); `base` carries the domain, PSF model, threshold
/// ([`SeedSpec::MatchedFilter`], else the default), detection cap, fit parameters and mode. `aperture_correction` of the
/// result is `None` (`c` = 1: no fit's flux is corrected); the green-site `c` is reported, not applied, in
/// `reported_aperture_correction` / `_stars` / `_se`. Every seed the green-sum matched filter found is handed to the fit in native pixels with
/// its `extended` flag (R-PC-8); `PlaneMeasurement::detected` stays empty (given seeds), `seeds` / `extended_seeds`
/// count them.
///
/// # Panics
/// When `mosaic.len() != w * h` (the reader guarantees it).
pub fn measure_green_sites(
    mosaic: &[f32],
    w: usize,
    h: usize,
    pattern: BayerPattern,
    base: &MeasureOptions,
    pool: Option<&Arc<rayon::ThreadPool>>,
) -> PlaneMeasurement {
    measure_green(mosaic, w, h, pattern, base, true, false, pool).0
}

/// [`measure_green_sites`] that also hands back the native seeds it handed the fit (every detected seed, `extended`
/// included, in native pixels): what the list entry ([`super::analysis_mode`]) keeps out of
/// a list star's ring sky. The measurement is the same bit for bit. Crate-private.
pub(crate) fn measure_green_sites_and_seeds(
    mosaic: &[f32],
    w: usize,
    h: usize,
    pattern: BayerPattern,
    base: &MeasureOptions,
    pool: Option<&Arc<rayon::ThreadPool>>,
) -> (PlaneMeasurement, Vec<Seed>) {
    measure_green(mosaic, w, h, pattern, base, true, true, pool)
}

/// [`measure_green_sites`]' body. `report_aperture` selects the green-site `c` (reported, never applied); `false` is the
/// path as it was, which `reporting_moves_no_fit` measures the reporting path against. `keep_seeds` returns a copy of
/// the native seeds beside the measurement (empty otherwise); nothing it does reaches a number.
#[allow(clippy::too_many_arguments)]
fn measure_green(
    mosaic: &[f32],
    w: usize,
    h: usize,
    pattern: BayerPattern,
    base: &MeasureOptions,
    report_aperture: bool,
    keep_seeds: bool,
    pool: Option<&Arc<rayon::ThreadPool>>,
) -> (PlaneMeasurement, Vec<Seed>) {
    assert_eq!(mosaic.len(), w * h, "mosaic does not match its geometry");
    let threshold = match base.seeds {
        SeedSpec::MatchedFilter { threshold } => threshold,
        _ => DEFAULT_THRESHOLD,
    };

    // Sky and noise: one green phase, in the working units.
    let (phase, pw, ph) = green_phase(mosaic, w, h, pattern);
    let sky = domain::sky_level(&phase, pw, ph, &base.domain, base.scale_estimator);
    let phase_working = domain::to_working_units(&phase, &base.domain);
    let (noise, _) = domain::noise_model(&phase_working, pw, ph, pool);
    if noise.source == NoiseSource::None {
        error!(
            width = w,
            height = h,
            "cfa green sites: no noise estimate on the green phase"
        );
    }

    // Detection: the G1 + G2 sum, white noise by construction. The `fieldWidthAtFloor` warn fires in `detect`, in
    // super-pixel units; the handed width does not carry the flag.
    let (sum, sw, sh) = green_sum(mosaic, w, h, pattern);
    let sum_working = domain::to_working_units(&sum, &base.domain);
    let (sum_noise, _) = domain::noise_model(&sum_working, sw, sh, pool);
    let (found, width, report) = seeds::detect(
        &sum_working,
        sw,
        sh,
        &sum_noise,
        threshold,
        base.max_stars,
        None,
        pool,
        PixelNoise::White,
    );
    // Knots (ruling A3, amended by R-PC-8) ride with their flag: `SeedSpec::Given` + `given_extended` fits them like the
    // mono path and keeps them out of `c`'s qualifiers and the frame statistic. Then native pixels, `native_seed`.
    let native: Vec<Seed> = found.iter().map(native_seed).collect();
    let kept_seeds = if keep_seeds { native.clone() } else { Vec::new() };
    let given_extended = report.extended;
    let field_width = FieldWidth {
        fwhm_px: 2.0 * width.fwhm_px,
        from_stars: width.from_stars,
    };

    let opts = MeasureOptions {
        domain: base.domain,
        psf_model: base.psf_model,
        max_stars: base.max_stars,
        seeds: SeedSpec::Given(native),
        fit: base.fit,
        noise: Some(noise),
        scale_estimator: base.scale_estimator,
        diagnostics: base.diagnostics,
        mask: Some(green_mask(w, h, pattern)),
        field_width: Some(field_width),
        // Given seeds: nothing is detected, so there is no detection width to tell.
        detect_width_hint: None,
        colour: false,
        pixel_noise: PixelNoise::White,
        cfa_site_share: None,
        mode: base.mode,
        sky: Some(sky),
        given_extended,
    };
    let pm = if report_aperture {
        star_measure::measure_plane_reporting_aperture(mosaic, w, h, &opts, pool)
    } else {
        star_measure::measure_plane(mosaic, w, h, &opts, pool)
    };
    (pm, kept_seeds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::star_measure::{Domain, Mode, PsfModel};
    use crate::test_support::{add_noise, moffat_field, MoffatStar};

    const LAYOUTS: [BayerPattern; 4] = [
        BayerPattern::Rggb,
        BayerPattern::Bggr,
        BayerPattern::Grbg,
        BayerPattern::Gbrg,
    ];

    /// 'R', 'G' or 'B' at `(x, y)`.
    fn site(p: BayerPattern, x: usize, y: usize) -> char {
        let (ex, ey) = (x % 2 == 0, y % 2 == 0);
        match (p, ex, ey) {
            (BayerPattern::Rggb, true, true) | (BayerPattern::Grbg, false, true) => 'R',
            (BayerPattern::Rggb, false, false) | (BayerPattern::Grbg, true, false) => 'B',
            (BayerPattern::Bggr, true, true) | (BayerPattern::Gbrg, false, true) => 'B',
            (BayerPattern::Bggr, false, false) | (BayerPattern::Gbrg, true, false) => 'R',
            _ => 'G',
        }
    }

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

    /// A 256x256-style mosaic: scene stars (flux 1.0 on G) x site gains + sky per site + Gaussian noise.
    fn mosaic(p: BayerPattern, w: usize, h: usize, stars: &[MoffatStar], seed: u64) -> Vec<f32> {
        let scene = moffat_field(w, h, stars, 4.0, 0.0);
        let mut m: Vec<f32> = (0..w * h)
            .map(|i| {
                let (gain, sky) = match site(p, i % w, i / w) {
                    'R' => (0.5, 900.0),
                    'B' => (0.4, 800.0),
                    _ => (1.0, 1000.0),
                };
                scene[i] * gain + sky
            })
            .collect();
        add_noise(&mut m, 6.0, seed);
        m
    }

    fn field_stars(bright: (f64, f64), w: usize, h: usize, seed: u64) -> Vec<MoffatStar> {
        let mut rng = crate::geometry::ransac::SplitMix64(seed);
        let mut stars = vec![star(bright.0, bright.1, 3000.0, 4.0)];
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
        stars
    }

    fn base() -> MeasureOptions {
        MeasureOptions {
            domain: Domain::RawAdu {
                clip: 65535.0,
                full_scale: 65535.0,
                gain_e_per_adu: None,
            },
            psf_model: PsfModel::Moffat4,
            seeds: SeedSpec::MatchedFilter {
                threshold: DEFAULT_THRESHOLD,
            },
            mode: Mode::Budget { budget: 500 },
            ..Default::default()
        }
    }

    #[test]
    fn green_mask_marks_the_green_sites_of_each_layout() {
        let at = |m: &[bool], x: usize, y: usize| m[y * 4 + x];
        for (p, g00) in [
            (BayerPattern::Rggb, false),
            (BayerPattern::Bggr, false),
            (BayerPattern::Grbg, true),
            (BayerPattern::Gbrg, true),
        ] {
            let m = green_mask(4, 6, p);
            assert_eq!(m.len(), 24);
            assert_eq!(at(&m, 0, 0), g00, "{p:?} (0,0)");
            assert_eq!(at(&m, 1, 0), !g00, "{p:?} (1,0)");
            assert_eq!(at(&m, 0, 1), !g00, "{p:?} (0,1)");
            assert_eq!(at(&m, 1, 1), g00, "{p:?} (1,1)");
            assert_eq!(m.iter().filter(|&&b| b).count(), 12);
            for (i, &g) in m.iter().enumerate() {
                assert_eq!(
                    g,
                    site(p, i % 4, i / 4) == 'G',
                    "{p:?} site ({}, {})",
                    i % 4,
                    i / 4
                );
            }
        }
    }

    #[test]
    fn a_star_on_a_mosaic_is_fitted_at_its_native_position() {
        let truth = (100.3, 80.7);
        for p in LAYOUTS {
            let stars = field_stars(truth, 256, 256, 11);
            let m = mosaic(p, 256, 256, &stars, 5);
            let pm = measure_green_sites(&m, 256, 256, p, &base(), None);
            let best = pm
                .stars
                .iter()
                .min_by(|a, b| {
                    (a.x - truth.0)
                        .hypot(a.y - truth.1)
                        .total_cmp(&(b.x - truth.0).hypot(b.y - truth.1))
                })
                .unwrap_or_else(|| panic!("{p:?}: no accepted fit; refused {:?}", pm.refused));
            let d = (best.x - truth.0).hypot(best.y - truth.1);
            assert!(
                d < 0.1,
                "{p:?}: fit at ({}, {}) is {d} px from the star",
                best.x,
                best.y
            );
            let fwhm = 0.5 * (best.fwhm_x + best.fwhm_y);
            assert!((fwhm / 4.0 - 1.0).abs() < 0.02, "{p:?}: FWHM {fwhm}");
        }
    }

    #[test]
    fn the_cfa_sky_and_noise_read_only_green_sites() {
        let p = BayerPattern::Rggb;
        let stars = field_stars((100.3, 80.7), 256, 256, 11);
        let a = mosaic(p, 256, 256, &stars, 5);
        let mut b = a.clone();
        for (i, v) in b.iter_mut().enumerate() {
            if site(p, i % 256, i / 256) != 'G' {
                *v += 5000.0;
            }
        }
        let (ma, mb) = (
            measure_green_sites(&a, 256, 256, p, &base(), None),
            measure_green_sites(&b, 256, 256, p, &base(), None),
        );
        assert_eq!(ma.sky.median.to_bits(), mb.sky.median.to_bits());
        assert_eq!(ma.noise.sigma.to_bits(), mb.noise.sigma.to_bits());
        assert!(!ma.stars.is_empty());
        // Seeds come from the two green sites only (R-PC-3), so R / B never reach a seed, the fit, the sky or the noise:
        // the accepted fits are identical bit for bit.
        assert_eq!(ma.seeds, mb.seeds);
        assert_eq!(ma.field_width, mb.field_width);
        // The 3 F aperture sum (Plan D D2c, PD-30) is NaN on a mosaic by design: its disc holds the absent R / B sites. A NaN
        // is not equal to itself, so it is checked on both sides and left out of the comparison of the rest.
        assert!(ma.stars.iter().chain(&mb.stars).all(|f| f.ap_sum.is_nan()));
        use crate::star_measure::StarFit;
        let without_sum = |s: &[StarFit]| -> Vec<StarFit> {
            s.iter().map(|f| StarFit { ap_sum: 0.0, ..*f }).collect()
        };
        assert_eq!(without_sum(&ma.stars), without_sum(&mb.stars));
    }

    #[test]
    fn a_hot_red_site_is_not_a_star() {
        let p = BayerPattern::Rggb;
        let mut m = mosaic(p, 128, 128, &[], 3);
        let (hx, hy) = (60, 40); // both even: an R site in RGGB
        assert_eq!(site(p, hx, hy), 'R');
        m[hy * 128 + hx] = 60000.0;
        let pm = measure_green_sites(&m, 128, 128, p, &base(), None);
        assert!(
            !pm.stars
                .iter()
                .any(|s| (s.x - hx as f64).hypot(s.y - hy as f64) < 3.0),
            "accepted fits near the hot red site: {:?}",
            pm.stars.iter().map(|s| (s.x, s.y)).collect::<Vec<_>>()
        );
        // A pin: the green-only seeds never see the site, so no seed was even attempted there (not merely refused
        // `NarrowerThanPixel`).
        assert!(
            !native_seeds(&m, 128, 128, p)
                .0
                .iter()
                .any(|s| (s.x - hx as f64).hypot(s.y - hy as f64) < 3.0),
            "a seed at the hot red site"
        );
    }

    #[test]
    fn an_odd_mosaic_leaves_its_last_row_and_column_out_of_detection_only() {
        let (w, h) = (257, 255);
        let (sum, sw, sh) = green_sum(&vec![1.0; w * h], w, h, BayerPattern::Rggb);
        assert_eq!((sw, sh, sum.len()), (128, 127, 128 * 127));
        assert!(sum.iter().all(|&v| v == 2.0));
        // A pin of the stride (`w`, not `2 * sw`) and of the unshifted grid: a wrong stride would shear every row and
        // shift this star's seed. The star sits ~10 px from the odd column (the module's border rule), so the test does
        // not use the last column itself; the fit lands at its native position.
        let p = BayerPattern::Rggb;
        let truth = (246.2, 100.4);
        let mut stars = field_stars((100.3, 80.7), w, h, 11);
        stars.push(star(truth.0, truth.1, 3000.0, 4.0));
        let m = mosaic(p, w, h, &stars, 5);
        let pm = measure_green_sites(&m, w, h, p, &base(), None);
        assert!(
            pm.stars
                .iter()
                .any(|s| (s.x - truth.0).hypot(s.y - truth.1) < 0.3),
            "no fit at the last-column star; refused {:?}",
            pm.refused
        );
    }

    /// R-PC-10b: the green-site Analysis path (`Mode::Budget`, through `measure_plane`) keeps the Blend refusal. A
    /// resolved pair (1.75 F apart, both seeded on the green sum) has a fit refused `Blend` on every layout, and
    /// neither star of it is an accepted fit. Plan D D3 (PD-10): the same mosaic under `Mode::All` (the stacking
    /// Measure's mode) is refused the same way; before D3 it refused none as `Blend`. (The pair's other fit is refused
    /// `RegionTooSmall`: the merged blob outgrows the region.)
    #[test]
    fn the_green_site_analysis_path_still_refuses_a_resolved_pair_blend() {
        use crate::star_measure::RefusalReason;
        let pair = [(180.0, 150.0), (187.0, 150.0)];
        for p in LAYOUTS {
            let mut stars = field_stars((100.3, 80.7), 256, 256, 11);
            stars.retain(|s| {
                pair.iter()
                    .all(|&(px, py)| (s.x - px).hypot(s.y - py) > 25.0)
            });
            for &(x, y) in &pair {
                stars.push(star(x, y, 700.0, 4.0));
            }
            let m = mosaic(p, 256, 256, &stars, 5);
            let blends = |pm: &PlaneMeasurement| {
                pm.refused
                    .by_reason
                    .get(&RefusalReason::Blend)
                    .copied()
                    .unwrap_or(0)
            };
            let near_pair = |pm: &PlaneMeasurement| {
                pm.stars
                    .iter()
                    .filter(|f| (f.x - 183.5).abs() < 8.0 && (f.y - 150.0).abs() < 4.0)
                    .count()
            };
            let budget = measure_green_sites(&m, 256, 256, p, &base(), None);
            assert!(blends(&budget) >= 1, "{p:?}: {:?}", budget.refused);
            assert_eq!(near_pair(&budget), 0, "{p:?}: {:?}", budget.stars);
            let all = measure_green_sites(
                &m,
                256,
                256,
                p,
                &MeasureOptions {
                    mode: Mode::All,
                    ..base()
                },
                None,
            );
            assert!(blends(&all) >= 1, "{p:?}: {:?}", all.refused);
            assert_eq!(near_pair(&all), 0, "{p:?}: {:?}", all.stars);
        }
    }

    /// The detection source of `measure_green_sites`: seeds on the G1 + G2 sum, in native pixels.
    fn native_seeds(m: &[f32], w: usize, h: usize, p: BayerPattern) -> (Vec<Seed>, Vec<bool>) {
        let b = base();
        let (sum, sw, sh) = green_sum(m, w, h, p);
        let working = domain::to_working_units(&sum, &b.domain);
        let (noise, _) = domain::noise_model(&working, sw, sh, None);
        let (found, _, report) = seeds::detect(
            &working,
            sw,
            sh,
            &noise,
            DEFAULT_THRESHOLD,
            b.max_stars,
            None,
            None,
            PixelNoise::White,
        );
        (found.iter().map(native_seed).collect(), report.extended)
    }

    /// I2: the 0.5 offset pinned where it is applied. A star at a block centre (2i + 0.5) or at a block corner
    /// (2i + 1.5) maps to within 0.15 px of its true position; an offset of 0.0 or 1.0 misses by 0.5 px per axis.
    #[test]
    fn the_super_pixel_offset_is_half_a_pixel_at_seed_level() {
        let truths = [
            (40.5, 40.5), // block centres
            (90.5, 50.5),
            (50.5, 90.5),
            (91.5, 91.5), // block corners
            (41.5, 71.5),
            (71.5, 41.5),
        ];
        let stars: Vec<MoffatStar> = truths
            .iter()
            .map(|&(x, y)| star(x, y, 3000.0, 4.0))
            .collect();
        for p in LAYOUTS {
            let m = mosaic(p, 128, 128, &stars, 5);
            let (seeds, _) = native_seeds(&m, 128, 128, p);
            for &(tx, ty) in &truths {
                let near = seeds
                    .iter()
                    .map(|s| (s.x - tx).abs().max((s.y - ty).abs()))
                    .fold(f64::INFINITY, f64::min);
                assert!(
                    near < 0.15,
                    "{p:?}: nearest seed to ({tx}, {ty}) is {near} px off"
                );
            }
        }
    }

    /// I1 / R-PC-8: a knot (an extended blob, ruling A3) is flagged by the matched filter on the green sum, HANDED to the
    /// fit flagged and fitted (as on the mono path), and never enters the frame statistic.
    #[test]
    fn a_knot_does_not_enter_the_cfa_statistic() {
        use crate::test_support::gaussian_field;
        let p = BayerPattern::Rggb;
        let (w, h) = (400, 300);
        let mut stars = Vec::new();
        for j in 0..3 {
            for i in 0..5 {
                stars.push(star(
                    50.0 + 80.0 * f64::from(i),
                    50.0 + 80.0 * f64::from(j) + 3.3,
                    700.0,
                    4.0,
                ));
            }
        }
        let plain = mosaic(p, w, h, &stars, 5);
        // A Gaussian blob of FWHM 12 px (3 F), ~40 px from its neighbours.
        let blob = (90.3, 94.8);
        let g = gaussian_field(w, h, &[(blob.0, blob.1, 300.0)], 12.0 / 2.354_820_045, 0.0);
        let mut knotted = plain.clone();
        for (i, v) in knotted.iter_mut().enumerate() {
            if site(p, i % w, i / w) == 'G' {
                *v += g[i];
            }
        }
        let near = |x: f64, y: f64| (x - blob.0).hypot(y - blob.1) < 8.0;
        // Premise: the matched filter seeds the blob and flags it.
        let (seeds, extended) = native_seeds(&knotted, w, h, p);
        let flagged = seeds
            .iter()
            .zip(&extended)
            .any(|(s, &e)| e && near(s.x, s.y));
        assert!(flagged, "premise: the blob is seeded and flagged extended");
        let (a, b) = (
            measure_green_sites(&plain, w, h, p, &base(), None),
            measure_green_sites(&knotted, w, h, p, &base(), None),
        );
        assert!(!a.stars.is_empty());
        assert!(b.detected.is_empty(), "given seeds leave `detected` empty");
        assert_eq!(b.seeds, seeds.len(), "every seed is handed to the fit");
        assert_eq!(b.extended_seeds, extended.iter().filter(|&&e| e).count());
        assert!(b.extended_seeds >= 1);
        let knot = b
            .stars
            .iter()
            .find(|s| near(s.x, s.y))
            .expect("the knot is handed and fitted");
        assert!(
            knot.flags.contains(crate::star_measure::StarFlags::EXTENDED),
            "{:?}",
            knot.flags
        );
        let (fa, fb) = (a.frame.fwhm_px.unwrap(), b.frame.fwhm_px.unwrap());
        assert!(
            (fb / fa - 1.0).abs() < 0.01,
            "frame FWHM {fb} vs {fa} without the knot"
        );
    }

    /// R-PC-8: a bright broad star (1.3 F at SNR 300, the injector's rung that the significance test flags) among F-width
    /// stars is fitted and flagged `extended`; the frame statistic is within 1 % of the run without it.
    #[test]
    fn a_bright_broad_star_on_a_mosaic_is_fitted_and_flagged() {
        let p = BayerPattern::Rggb;
        let (w, h) = (400, 300);
        let mut stars = Vec::new();
        for j in 0..3 {
            for i in 0..5 {
                stars.push(star(
                    50.0 + 80.0 * f64::from(i),
                    50.0 + 80.0 * f64::from(j) + 3.3,
                    700.0,
                    4.0,
                ));
            }
        }
        let plain = mosaic(p, w, h, &stars, 5);
        let broad = (90.3, 94.8);
        stars.push(star(broad.0, broad.1, 20000.0, 1.3 * 4.0));
        let with = mosaic(p, w, h, &stars, 5);
        let near = |x: f64, y: f64| (x - broad.0).hypot(y - broad.1) < 4.0;
        let (a, b) = (
            measure_green_sites(&plain, w, h, p, &base(), None),
            measure_green_sites(&with, w, h, p, &base(), None),
        );
        let star_b = b
            .stars
            .iter()
            .find(|s| near(s.x, s.y))
            .expect("the broad star is in pm.stars");
        assert!(
            star_b.flags.extended(),
            "premise: the significance test flags it ({:?}, extended seeds {})",
            star_b.flags,
            b.extended_seeds
        );
        let (fa, fb) = (a.frame.fwhm_px.unwrap(), b.frame.fwhm_px.unwrap());
        assert!(
            (fb / fa - 1.0).abs() < 0.01,
            "frame FWHM {fb} vs {fa} without it"
        );
    }

    #[test]
    fn a_cfa_plane_reports_its_aperture_correction() {
        let p = BayerPattern::Rggb;
        let stars = field_stars((100.3, 80.7), 256, 256, 11);
        let m = mosaic(p, 256, 256, &stars, 5);
        let pm = measure_green_sites(&m, 256, 256, p, &base(), None);
        assert!(!pm.stars.is_empty());
        // Applied: nothing (O3). Reported: the green-site `c`.
        assert_eq!(pm.aperture_correction, None);
        assert_eq!(pm.aperture_stars, 0);
        assert_eq!(pm.aperture_se, None);
        assert!(pm.aperture_ratios.is_empty());
        let c = pm
            .reported_aperture_correction
            .unwrap_or_else(|| panic!("no reported c; {} stars", pm.reported_aperture_stars));
        assert!(pm.reported_aperture_stars >= 5, "{}", pm.reported_aperture_stars);
        assert!(c.is_finite() && c > 1.0 && c < 1.5, "c = {c}");
        assert!(pm.reported_aperture_se.is_some_and(|se| se.is_finite() && se >= 0.0));
        // No fit's signal carries it: the model flux in the FWTM, as at c = 1.
        for f in &pm.stars {
            let model = crate::star_measure::photometry::model_flux_in_fwtm(f);
            assert_eq!(f.signal.to_bits(), model.to_bits(), "signal is the bare model flux");
        }
    }

    /// A mosaic of pixel-integrated Moffat stars (β 2.5 wings): `(field, stars (x, y, amplitude))` on an RGGB mosaic with
    /// R / B gains 0.5 / 0.4 and green sites of gain `g1` on even rows and `g2` on odd rows, sky 1000 per site (+ `ped2`
    /// on the odd-row greens), Gaussian noise 4.
    fn integrated_mosaic(
        fwhm: f64,
        gains: (f64, f64),
        ped2: f64,
        seed: u64,
    ) -> (Vec<f32>, Vec<(f64, f64, f64)>) {
        integrated_mosaic_beta(fwhm, 2.5, gains, ped2, seed)
    }

    /// [`integrated_mosaic`] with the stars' wing exponent `beta` stated.
    fn integrated_mosaic_beta(
        fwhm: f64,
        beta: f64,
        (g1, g2): (f64, f64),
        ped2: f64,
        seed: u64,
    ) -> (Vec<f32>, Vec<(f64, f64, f64)>) {
        use crate::test_support::{integrated_field, IntegratedProfile};
        let (w, h) = (400usize, 400usize);
        let mut rng = crate::geometry::ransac::SplitMix64(seed);
        let mut stars: Vec<(f64, f64, f64)> = Vec::new();
        let mut attempts = 0;
        while stars.len() < 30 && attempts < 20000 {
            attempts += 1;
            let x = 45.0 + rng.next_f64() * (w as f64 - 90.0);
            let y = 45.0 + rng.next_f64() * (h as f64 - 90.0);
            if stars.iter().any(|s| (s.0 - x).hypot(s.1 - y) < 50.0) {
                continue;
            }
            stars.push((x, y, 1500.0 + rng.next_f64() * 3500.0));
        }
        let scene = integrated_field(
            w,
            h,
            &stars,
            IntegratedProfile::Moffat { fwhm, beta },
            5,
            0.0,
        );
        let mut m: Vec<f32> = (0..w * h)
            .map(|i| {
                let (x, y) = (i % w, i / w);
                let (gain, sky) = match site(BayerPattern::Rggb, x, y) {
                    'R' => (0.5, 900.0),
                    'B' => (0.4, 800.0),
                    _ if y % 2 == 0 => (g1, 1000.0),
                    _ => (g2, 1000.0 + ped2),
                };
                (scene[i] as f64 * gain + sky) as f32
            })
            .collect();
        add_noise(&mut m, 4.0, seed + 1);
        (m, stars)
    }

    /// A round Moffat's analytic flux inside radius `r` (peak `amp`, FWHM `fwhm`, exponent `beta`).
    fn moffat_flux_within(amp: f64, fwhm: f64, beta: f64, r: f64) -> f64 {
        let alpha = fwhm / (2.0 * (2f64.powf(1.0 / beta) - 1.0).sqrt());
        std::f64::consts::PI * alpha * alpha * amp / (beta - 1.0)
            * (1.0 - (1.0 + r * r / (alpha * alpha)).powf(1.0 - beta))
    }

    /// One star's growth flux at `5F` on a mosaic's green sites, twice: `(old, new)`. `old` is task 1's (the green sum
    /// over the constant lattice density 0.5), `new` the phase-matched growth the reported `c` is made of. `f_ring` is
    /// the plane's `F`, B_i the fit's own ring sky.
    fn green_growth_old_new(
        m: &[f32],
        w: usize,
        h: usize,
        p: BayerPattern,
        f: &crate::star_measure::StarFit,
        f_ring: f64,
    ) -> Option<(f64, f64)> {
        let mask = green_mask(w, h, p);
        let mut scratch = Vec::new();
        let new = crate::star_measure::photometry::growth_curve_in(
            m,
            Some(&mask),
            None,
            w,
            h,
            f.x,
            f.y,
            f_ring,
            f.sky_ring,
            Some(f),
            &mut scratch,
        )?;
        let r2 = (crate::star_measure::photometry::GROWTH_RADIUS_PER_FWHM * f_ring).powi(2);
        let sum: f64 = (0..w * h)
            .filter(|&i| {
                mask[i] && ((i % w) as f64 - f.x).powi(2) + ((i / w) as f64 - f.y).powi(2) <= r2
            })
            .map(|i| m[i] as f64 - f.sky_ring)
            .sum();
        Some((sum / 0.5, new))
    }

    /// 16 pixel-integrated Moffat stars (peak `amp` over a sky of 1000 per site, noise 3) on an RGGB mosaic, one at each
    /// of the checkerboard's 4 x 4 sub-pixel phases `((i + 0.5) / 4, (j + 0.5) / 4)`, 60 px apart (an even pitch keeps
    /// the parity of the integer part). Returns the mosaic and the stars `(x, y, peak)`.
    fn phase_grid_mosaic(fwhm: f64, beta: f64, amp: f64, seed: u64) -> (Vec<f32>, Vec<(f64, f64, f64)>) {
        use crate::test_support::{integrated_field, IntegratedProfile};
        let (w, h) = (300usize, 300usize);
        let mut stars = Vec::new();
        for j in 0..4 {
            for i in 0..4 {
                // The peaks differ by a few percent each, so no two pixels share the plane's maximum.
                let k = (j * 4 + i) as f64;
                stars.push((
                    42.0 + 60.0 * i as f64 + (i as f64 + 0.5) / 4.0,
                    42.0 + 60.0 * j as f64 + (j as f64 + 0.5) / 4.0,
                    amp * (1.0 - 0.012 * k),
                ));
            }
        }
        let scene = integrated_field(
            w,
            h,
            &stars,
            IntegratedProfile::Moffat { fwhm, beta },
            5,
            0.0,
        );
        let mut m: Vec<f32> = scene.iter().map(|&v| (v as f64 + 1000.0) as f32).collect();
        add_noise(&mut m, 3.0, seed);
        (m, stars)
    }

    /// OSC zero point, spec §5.1 as amended by task 1b: the green sum of a sharp star depends on where its centre falls
    /// on the checkerboard; the phase-matched growth (the green sum scaled by the star's OWN fitted model's all-sites
    /// over green-sites sum) does not. Pixel-integrated Moffat stars at FWHM 1.25 to 3 px on a 4 x 4 phase grid, measured
    /// through the mosaic path, Moffat 4 fits, `c_i = growth / model flux in the FWTM` per star.
    ///
    /// The per-star phase spread of `c_i` (largest `|c_i / median − 1|` over the 16 phases) is what `c` is made of, and
    /// it falls from task 1's 10 % at FWHM 1.25 to under 3 % (and from 6 % to under 3 % at 1.4): old against new is
    /// printed. The growth against the analytic flux inside `5F` is not the bar for wings the model has: the model is
    /// fitted on the green sites alone, so `S_d / S_mg` is near 1 and the growth is the MODEL's own total, which a
    /// point-sampled fit to an undersampled star carries a phase error of (±6 % at FWHM 1.25, in `signal` itself and not
    /// in `c`). For wings of β 2.5 that error is what is left: within 2 % per phase at FWHM 1.25, median within 1 %.
    #[test]
    fn phase_matched_growth_removes_the_checkerboard_phase_from_c() {
        let p = BayerPattern::Rggb;
        let spread = |v: &[f64]| {
            let c = crate::star_measure::tests::median(v.iter().copied());
            v.iter().fold(0.0f64, |a, x| a.max((x / c - 1.0).abs()))
        };
        for beta in [4.0, 2.5] {
            for fwhm in [1.25, 1.4, 1.6, 2.0, 3.0] {
                let amp = 8000.0;
                let (m, stars) = phase_grid_mosaic(fwhm, beta, amp, 17);
                let pm = measure_green_sites(&m, 300, 300, p, &base(), None);
                let f_ring = pm.frame.fwhm_px.unwrap_or_else(|| {
                    panic!(
                        "β {beta} FWHM {fwhm}: no frame FWHM; seeds {}, fits {}, refused {:?}, frame {:?}",
                        pm.seeds,
                        pm.stars.len(),
                        pm.refused,
                        pm.frame,
                    )
                });
                let r_ap = crate::star_measure::photometry::GROWTH_RADIUS_PER_FWHM * f_ring;
                let (mut old_errs, mut new_errs) = (Vec::new(), Vec::new());
                let (mut cs_old, mut cs_new): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
                for &(sx, sy, peak) in &stars {
                    let Some((f, (old, new))) = pm
                        .stars
                        .iter()
                        .find(|f| (f.x - sx).hypot(f.y - sy) < 0.7)
                        .and_then(|f| {
                            green_growth_old_new(&m, 300, 300, p, f, f_ring).map(|g| (f, g))
                        })
                    else {
                        continue;
                    };
                    let truth = moffat_flux_within(peak, fwhm, beta, r_ap);
                    let model = crate::star_measure::photometry::model_flux_in_fwtm(f);
                    old_errs.push(old / truth - 1.0);
                    new_errs.push(new / truth - 1.0);
                    cs_old.push(old / model);
                    cs_new.push(new / model);
                }
                assert!(
                    new_errs.len() >= 14,
                    "β {beta} FWHM {fwhm}: only {} of 16 stars fitted and grown",
                    new_errs.len()
                );
                let max_abs = |v: &[f64]| v.iter().fold(0.0f64, |a, e| a.max(e.abs()));
                let med = crate::star_measure::tests::median(new_errs.iter().copied());
                let (so, sn) = (spread(&cs_old), spread(&cs_new));
                eprintln!(
                    "phase grid β {beta} FWHM {fwhm} (F {f_ring:.3}, {} phases): c phase spread old (0.5) {so:.4} -> new {sn:.4}; growth/analytic 5F max |error| old {:.4} -> new {:.4}, median new {med:+.5}",
                    new_errs.len(),
                    max_abs(&old_errs),
                    max_abs(&new_errs)
                );
                assert!(sn < 0.03, "β {beta} FWHM {fwhm}: c phase spread {sn:.4} (bar 0.03)");
                if fwhm <= 1.4 {
                    assert!(
                        sn < 0.5 * so,
                        "β {beta} FWHM {fwhm}: c phase spread new {sn:.4} is not under half of the old {so:.4}"
                    );
                } else {
                    assert!(
                        sn < so + 0.005,
                        "β {beta} FWHM {fwhm}: c phase spread new {sn:.4} vs old {so:.4}"
                    );
                }
                assert!(
                    med.abs() < 0.01,
                    "β {beta} FWHM {fwhm}: median growth / analytic flux {med:+.4} (bar 0.01)"
                );
                if beta < 3.0 && fwhm <= 1.25 {
                    assert!(
                        max_abs(&new_errs) < 0.02,
                        "β {beta} FWHM {fwhm}: per-phase growth / analytic flux |error| {:.4} (bar 0.02)",
                        max_abs(&new_errs)
                    );
                }
            }
        }
    }

    /// OSC zero point, O2 / O3 with task 1b: on pixel-integrated Moffat (β 2.5) stars fitted with Moffat 4, `signal ·
    /// reported c` is the star's true flux: the median over stars within 1 %, at FWHM 1.3, 1.6, 2 and 3 px, and with the
    /// two green sites at gains 1.00 / 0.98 and a +3 ADU pedestal on the second (the truth then is the mean green gain's,
    /// 0.99, of the scene). At FWHM 1.3 the per-star spread of the growth over the truth is printed old (task 1's constant
    /// 0.5) against new.
    #[test]
    fn mosaic_corrected_flux_matches_the_truth() {
        let p = BayerPattern::Rggb;
        let beta = 2.5;
        for (fwhm, gains, ped2) in [
            (1.3, (1.0, 1.0), 0.0),
            (1.6, (1.0, 1.0), 0.0),
            (2.0, (1.0, 1.0), 0.0),
            (3.0, (1.0, 1.0), 0.0),
            (2.0, (1.0, 0.98), 3.0),
            (3.0, (1.0, 0.98), 3.0),
        ] {
            let (m, stars) = integrated_mosaic(fwhm, gains, ped2, 31);
            let pm = measure_green_sites(&m, 400, 400, p, &base(), None);
            let c = pm.reported_aperture_correction.unwrap_or_else(|| {
                panic!("FWHM {fwhm}: no c; stars {}", pm.reported_aperture_stars)
            });
            let alpha = fwhm / (2.0 * (2f64.powf(1.0 / beta) - 1.0).sqrt());
            let mean_gain = 0.5 * (gains.0 + gains.1);
            let f_ring = pm.frame.fwhm_px.expect("frame FWHM");
            let (mut ratios, mut old_g, mut new_g): (Vec<f64>, Vec<f64>, Vec<f64>) =
                (Vec::new(), Vec::new(), Vec::new());
            for f in &pm.stars {
                let Some(t) = stars
                    .iter()
                    .find(|t| (t.0 - f.x).hypot(t.1 - f.y) < 1.5)
                else {
                    continue;
                };
                let total = std::f64::consts::PI * alpha * alpha * t.2 / (beta - 1.0) * mean_gain;
                ratios.push(f.signal * c / total);
                if let Some((old, new)) = green_growth_old_new(&m, 400, 400, p, f, f_ring) {
                    let model = crate::star_measure::photometry::model_flux_in_fwtm(f);
                    old_g.push(f.signal * (old / model) / total);
                    new_g.push(f.signal * (new / model) / total);
                }
            }
            assert!(ratios.len() >= 10, "FWHM {fwhm}: {} stars matched", ratios.len());
            let med = crate::star_measure::tests::median(ratios.iter().copied());
            let mad = |v: &[f64]| {
                let c = crate::star_measure::tests::median(v.iter().copied());
                crate::star_measure::tests::median(v.iter().map(|x| (x - c).abs()))
            };
            eprintln!(
                "mosaic truth, FWHM {fwhm} gains {gains:?} ped {ped2}: c {c:.4} ({} stars, se {:?}), median signal*c/true {med:.4} over {} stars; per-star growth/true MAD old {:.4} new {:.4} (medians {:.4} / {:.4}, {} stars)",
                pm.reported_aperture_stars,
                pm.reported_aperture_se,
                ratios.len(),
                mad(&old_g),
                mad(&new_g),
                crate::star_measure::tests::median(old_g.iter().copied()),
                crate::star_measure::tests::median(new_g.iter().copied()),
                new_g.len()
            );
            assert!(
                (med - 1.0).abs() < 0.01,
                "FWHM {fwhm} gains {gains:?} pedestal {ped2}: median signal*c / true flux = {med:.4}"
            );
        }
    }

    /// Every field the old path produced, as text (`Debug` prints a float's exact digits and NaN as NaN).
    fn fits_text(pm: &PlaneMeasurement) -> Vec<String> {
        pm.stars.iter().map(|f| format!("{f:?}")).collect()
    }

    fn assert_same_but_the_report(a: &PlaneMeasurement, b: &PlaneMeasurement, what: &str) {
        assert_eq!(fits_text(a), fits_text(b), "{what}: every StarFit");
        assert_eq!(format!("{:?}", a.frame), format!("{:?}", b.frame), "{what}: frame aggregates");
        assert_eq!(format!("{:?}", a.noise), format!("{:?}", b.noise), "{what}: noise");
        assert_eq!(format!("{:?}", a.sky), format!("{:?}", b.sky), "{what}: sky");
        assert_eq!(a.beta.to_bits(), b.beta.to_bits(), "{what}: beta");
        assert_eq!(a.seeds, b.seeds, "{what}: seeds");
        assert_eq!(a.refused, b.refused, "{what}: refusals");
        assert_eq!(a.aperture_correction, b.aperture_correction, "{what}: c");
        assert_eq!(a.aperture_stars, b.aperture_stars, "{what}: c stars");
        assert_eq!(a.aperture_se, b.aperture_se, "{what}: c se");
        assert_eq!(a.aperture_ratios.len(), b.aperture_ratios.len(), "{what}: ratios");
        assert_eq!(a.point_source_screened, b.point_source_screened, "{what}: screens");
        assert_eq!(a.aperture_below_snr, b.aperture_below_snr, "{what}: snr screen");
        assert_eq!(a.aperture_rescued, b.aperture_rescued, "{what}: rescued");
        assert_eq!(a.ring_sky_missing, b.ring_sky_missing, "{what}: ring sky");
        assert_eq!(a.clip_floor.map(f64::to_bits), b.clip_floor.map(f64::to_bits), "{what}: clip");
        assert_eq!(format!("{:?}", a.budget), format!("{:?}", b.budget), "{what}: budget");
    }

    /// OSC zero point, O3: reporting `c` moves nothing. One mosaic through the reporting path and through the path as it
    /// was (`FromPlane`): every `StarFit` field (`ap_sum` NaN on both), the frame aggregates and every other output
    /// equal; the same for a mono plane and a debayered 3-plane file, where the reporting path is not selected.
    #[test]
    fn reporting_moves_no_fit() {
        use crate::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
        use crate::star_measure::raw::{RawFrame, RawPlanes};
        let p = BayerPattern::Rggb;
        let (m, _) = integrated_mosaic(2.4, (1.0, 0.98), 3.0, 77);
        let new = measure_green_sites(&m, 400, 400, p, &base(), None);
        let old = measure_green(&m, 400, 400, p, &base(), false, false, None).0;
        assert!(!new.stars.is_empty());
        assert!(new.stars.iter().chain(&old.stars).all(|f| f.ap_sum.is_nan()));
        assert_same_but_the_report(&new, &old, "mosaic");
        assert!(new.reported_aperture_correction.is_some());
        assert_eq!(
            (old.reported_aperture_correction, old.reported_aperture_stars, old.reported_aperture_se),
            (None, 0, None),
            "the path as it was reports nothing"
        );
        let fo = FrameOptions::default();
        assert_eq!(
            format!("{:?}", FrameSummary::from_measurement(&new, &fo)),
            format!("{:?}", FrameSummary::from_measurement(&old, &fo)),
            "FrameSummary"
        );

        // A mono plane and a debayered 3-plane file: not selected, nothing reported, the outputs are `measure_plane`'s.
        let stars = field_stars((100.3, 80.7), 256, 256, 11);
        let plane = {
            let mut d = moffat_field(256, 256, &stars, 4.0, 1000.0);
            add_noise(&mut d, 6.0, 9);
            d
        };
        for (what, planes, channels) in [
            ("mono", RawPlanes::Mono(plane.clone()), 1),
            (
                "rgb",
                RawPlanes::Rgb([plane.clone(), plane.clone(), plane.clone()]),
                3,
            ),
        ] {
            let frame = RawFrame {
                planes,
                w: 256,
                h: 256,
                domain: base().domain,
                flip_vertical: false,
                channels,
            };
            let pm = measure_raw_frame(&frame, &fo, None);
            assert!(!pm.stars.is_empty(), "{what}");
            assert_eq!(
                (pm.reported_aperture_correction, pm.reported_aperture_stars, pm.reported_aperture_se),
                (None, 0, None),
                "{what}: nothing reported"
            );
            assert!(pm.aperture_correction.is_some(), "{what}: its own c is applied, as before");
        }
    }
}
