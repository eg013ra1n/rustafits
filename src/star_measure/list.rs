//! The list entry's machinery (forced-list zero point, R-OSP-19): the fits of a caller's list of positions on a plane
//! this module has already measured. The public face is the list entry of [`super::analysis_mode`] and its
//! three types; everything here is crate-private.
//!
//! The list is fitted with the budget path's halves ([`super::measure_budget`]) on the same plane: the saturation plane
//! from the measurement's own noise and sky, the verdict pass over the list's seeds at the measurement's field width,
//! then [`fit::fit_after_passes`] at the measurement's β with the plane's mask. [`super::fit_stars_at`] is not used: it
//! runs no saturation verdict and masks no detected neighbour in the ring sky, and a list star needs both.

use std::collections::HashMap;
use std::sync::Arc;

use crate::star_measure::analysis_mode::{ListFit, ListMeasurement, ListRequest};
use crate::star_measure::fit::{self, Correction, Neighbours};
use crate::star_measure::{
    domain, install_in_span, mask_mismatch, photometry, prepare_passes, saturation_plane,
    ApertureNoise, Passes, Domain, FitParams, PlaneMeasurement, PsfModel, Seed, SnrNoise, StarFit, StarFlags,
};

/// The plane a list is fitted on, as the measurement read it: the pixels, the domain, the mask of the pixels that carry
/// the plane's own colour (a mosaic's green sites), whether the growth is masked and whether it is one plane of a
/// debayered colour image.
pub(crate) struct ListPlane<'a> {
    pub(crate) plane: &'a [f32],
    pub(crate) w: usize,
    pub(crate) h: usize,
    pub(crate) domain: Domain,
    pub(crate) mask: Option<&'a [bool]>,
    pub(crate) mosaic: bool,
    pub(crate) colour: bool,
}

/// The seeds of `detected` that are not within `tol` px of a seed of `listed`: a list star's own detection sits at
/// distance about 0 and would fail `c`'s isolation test, so no list star would ever qualify for `c`. A cell grid of side
/// `tol` over `listed`, so the cost is linear.
fn without_coinciding(detected: &[Seed], listed: &[Seed], tol: f64) -> Vec<Seed> {
    let cell = |x: f64, y: f64| ((x / tol).floor() as i64, (y / tol).floor() as i64);
    let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (j, s) in listed.iter().enumerate() {
        grid.entry(cell(s.x, s.y)).or_default().push(j);
    }
    detected
        .iter()
        .filter(|d| {
            let (cx, cy) = cell(d.x, d.y);
            !(cx - 1..=cx + 1).any(|gx| {
                (cy - 1..=cy + 1).any(|gy| {
                    grid.get(&(gx, gy)).is_some_and(|at| {
                        at.iter()
                            .any(|&j| (listed[j].x - d.x).hypot(listed[j].y - d.y) <= tol)
                    })
                })
            })
        })
        .copied()
        .collect()
}

/// A list fit that was not made, or was refused: NaN numbers, both flags false.
fn refused_fit() -> ListFit {
    ListFit {
        x: f64::NAN,
        y: f64::NAN,
        model_flux: f64::NAN,
        flux: f64::NAN,
        snr: f64::NAN,
        flags: StarFlags(0),
        accepted: false,
        kept: false,
        centre_offset_px: f64::NAN,
    }
}

/// Fit `req`'s positions on `input` (`pm` is the plane's own measurement, `detected` its detected seeds in native
/// pixels: the matched filter's on a mono or G plane, the green-sum filter's on a mosaic) and measure `c` over the
/// kept fits. `None` (always with a `warn!` or `debug!` that says why) when the request cannot be measured: no
/// position, no usable tolerance or `c` cap, a plane or mask that does not match its geometry.
pub(crate) fn measure_list(
    input: &ListPlane<'_>,
    pm: &PlaneMeasurement,
    detected: &[Seed],
    req: &ListRequest,
    pool: Option<&Arc<rayon::ThreadPool>>,
) -> Option<ListMeasurement> {
    let (w, h) = (input.w, input.h);
    if req.positions.is_empty() {
        tracing::debug!("list: no position asked; nothing measured");
        return None;
    }
    let tol = req.max_centre_offset_px;
    if !(tol.is_finite() && tol > 0.0) || req.max_c_stars == 0 {
        tracing::warn!(
            max_centre_offset_px = tol,
            max_c_stars = req.max_c_stars,
            count = req.positions.len(),
            "list: the centre tolerance and the c cap must be positive; nothing measured"
        );
        return None;
    }
    if input.plane.len() != w * h {
        tracing::warn!(
            plane_len = input.plane.len(),
            width = w,
            height = h,
            "list: the plane does not match its geometry; nothing measured"
        );
        return None;
    }
    if mask_mismatch(input.mask, w, h) {
        return None;
    }
    let mask = input.mask;
    let working = domain::to_working_units(input.plane, &input.domain);
    let working: &[f32] = &working;
    let saturation = saturation_plane(working, mask, &pm.noise, &pm.sky, &input.domain, input.colour);

    // The seeds: a position on the plane (a refused one keeps its NaN row), the working plane at its nearest present
    // pixel less the sky median as its peak and flux, which the verdict pass reads for its plateau tests.
    let mut seeds: Vec<Seed> = Vec::with_capacity(req.positions.len());
    let mut request_of: Vec<usize> = Vec::with_capacity(req.positions.len());
    for (r, &(x, y)) in req.positions.iter().enumerate() {
        if !(x.is_finite() && y.is_finite() && x >= 0.0 && y >= 0.0 && x <= (w - 1) as f64 && y <= (h - 1) as f64) {
            continue;
        }
        let i = y.round() as usize * w + x.round() as usize;
        let present = working[i].is_finite() && mask.map_or(true, |m| m[i]);
        let peak = if present {
            (f64::from(working[i]) - saturation.bg).max(0.0)
        } else {
            0.0
        };
        seeds.push(Seed {
            x,
            y,
            peak,
            flux: peak,
        });
        request_of.push(r);
    }
    let outside = req.positions.len() - seeds.len();
    if outside > 0 {
        tracing::debug!(
            count = outside,
            width = w,
            height = h,
            "list: positions outside the plane or not finite; refused without a fit"
        );
    }
    let mut fits: Vec<ListFit> = (0..req.positions.len()).map(|_| refused_fit()).collect();
    if seeds.is_empty() {
        return Some(ListMeasurement {
            fits,
            aperture_correction: None,
            aperture_stars: 0,
            aperture_se: None,
            beta: pm.beta,
            ring_fwhm: None,
            refused: Default::default(),
        });
    }

    // The ring sky and `c`'s isolation read every detected seed except a list star's own detection, then the list.
    let mut neighbour_seeds = without_coinciding(detected, &seeds, tol);
    let first = neighbour_seeds.len();
    neighbour_seeds.extend_from_slice(&seeds);
    let of_fit: Vec<usize> = (first..first + seeds.len()).collect();

    let fit_params = FitParams::default();
    let passes = prepare_passes(
        working,
        w,
        h,
        &fit_params,
        mask,
        pool,
        &seeds,
        pm.field_width,
        saturation,
    );
    let Passes {
        fit_field,
        passes,
        saturation,
        ..
    } = passes;
    let rho = if pm.snr_noise == SnrNoise::Read {
        pm.detection
            .as_ref()
            .and_then(|d| d.noise_correlation.as_ref())
            .map(|c| c.rho())
    } else {
        None
    };
    let ring_fwhm = pm.frame.fwhm_px;
    let fit_call = || {
        fit::fit_after_passes(
            working,
            w,
            h,
            &seeds,
            PsfModel::Fixed(pm.beta),
            &fit_params,
            mask,
            fit_field,
            Some(&passes),
            Correction {
                saturation: Some(&saturation),
                ring_fwhm,
                gain: pm.snr_gain,
                noise: ApertureNoise {
                    rho,
                    shot_spread: pm.shot_spread,
                },
                neighbours: Some(Neighbours {
                    seeds: &neighbour_seeds,
                    of_fit: &of_fit,
                }),
                blend_refusal: true,
                report_aperture: true,
                masked_growth: input.mosaic,
                max_c_stars: req.max_c_stars,
                warn_short_c: false,
                ..Default::default()
            },
            None,
            false,
            false,
            false,
            &[],
        )
    };
    let outcome = install_in_span(pool, fit_call);

    // Kept: accepted, not saturated (the verdict pass's saturated and near-clip classes), centred within the tolerance.
    let mut kept_fits: Vec<(usize, StarFit)> = Vec::new();
    let mut accepted_unsaturated = 0usize;
    for (f, &fit_seed) in outcome.fits.iter().zip(&outcome.fit_seeds) {
        let r = request_of[fit_seed];
        let (px, py) = req.positions[r];
        let offset = (f.x - px).hypot(f.y - py);
        let kept = !f.flags.saturated() && offset.is_finite() && offset <= tol;
        accepted_unsaturated += usize::from(!f.flags.saturated());
        fits[r] = ListFit {
            x: f.x,
            y: f.y,
            model_flux: photometry::model_flux_in_fwtm(f),
            flux: f64::NAN,
            snr: f.snr,
            flags: f.flags,
            accepted: true,
            kept,
            centre_offset_px: offset,
        };
        if kept {
            kept_fits.push((fit_seed, *f));
        }
    }

    // `c` over the kept fits. The fit call measured it over every accepted unsaturated fit, which is the same set unless
    // a fit was centred off its position; then it is measured again over the kept ones alone.
    let (c, stars, se) = if kept_fits.len() == accepted_unsaturated {
        let r = outcome.reported_aperture;
        (r.correction, r.stars, r.se)
    } else {
        let (c, stars, se, _) = photometry::aperture_correction(
            working,
            mask,
            w,
            h,
            &kept_fits,
            &neighbour_seeds,
            Some(&of_fit),
            fit_field.fwhm_px,
            outcome.ring_fwhm,
            Some(&saturation),
            None,
            None,
            None,
            None,
            false,
            input.mosaic,
            req.max_c_stars,
            false,
        );
        (c, stars, se)
    };
    if let Some(c) = c {
        for lf in fits.iter_mut().filter(|lf| lf.accepted) {
            lf.flux = lf.model_flux * c;
        }
    }
    let kept = fits.iter().filter(|lf| lf.kept).count();
    tracing::debug!(
        count = req.positions.len(),
        kept,
        aperture_correction = c,
        aperture_stars = stars,
        "list measured"
    );
    Some(ListMeasurement {
        fits,
        aperture_correction: c,
        aperture_stars: stars,
        aperture_se: se,
        beta: pm.beta,
        ring_fwhm: outcome.ring_fwhm,
        refused: outcome.refused,
    })
}
