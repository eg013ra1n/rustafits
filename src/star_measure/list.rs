//! The list entry's machinery (forced-list zero point, R-OSP-19): the fits and the aperture sums of a caller's list of
//! positions on a plane this module has already measured. The public face is the list entry of [`super::analysis_mode`]
//! and its types; everything here is crate-private.
//!
//! Both fit passes use the budget path's halves ([`super::measure_budget`]) on the same plane: the saturation plane from
//! the measurement's own noise and sky, the verdict pass over the pass's seeds at the measurement's field width, then
//! [`fit::fit_after_passes`] at the measurement's β with the plane's mask. The near-clip class reads the MEASURED plane's
//! clip floor ([`PlaneMeasurement::clip_floor`]), never one re-derived from the list's seeds (a list holds few or no
//! clipped stars, so its own evidence would miss a floor the plane shows). [`super::fit_stars_at`] is not used: it runs
//! no saturation verdict and masks no detected neighbour in the ring sky, and a list star needs both. No aperture
//! correction is measured ([`super::FluxCorrection::Given`] of 1, so no growth curve runs; the ring sky and the SNR pass
//! still do). The sums are [`photometry::aperture_sums`], on the same working plane.

use std::borrow::Cow;
use std::sync::Arc;
use std::time::Instant;

use rayon::prelude::*;

use crate::geometry::kdtree::KdTree2;
use crate::star_measure::analysis_mode::{
    ListAperture, ListFit, ListMeasurement, ListPlacement, ListRequest, ListTimings,
};
use crate::star_measure::fit::{self, Correction, NeighbourSplit, Neighbours};
use crate::star_measure::photometry::SaturationPlane;
use crate::star_measure::{
    domain, install_in_span, mask_mismatch, photometry, prepare_passes, saturation_plane,
    ApertureNoise, Domain, FitParams, FluxCorrection, Passes, PlaneMeasurement, PsfModel,
    RefusalCounts, Seed, SnrNoise, StarFlags,
};

/// The plane a list is fitted on, as the measurement read it: the pixels, the domain and the mask of the pixels that
/// carry the plane's own colour (a mosaic's green sites; `None` = mono). A debayered colour file has no list (the entry
/// refuses it), so the plane is never one plane of a colour image.
pub(crate) struct ListPlane<'a> {
    pub(crate) plane: &'a [f32],
    pub(crate) w: usize,
    pub(crate) h: usize,
    pub(crate) domain: Domain,
    pub(crate) mask: Option<&'a [bool]>,
}

/// Which of `detected` is a list star's own detection: for every seed of `listed`, the NEAREST detected seed within
/// `radius` px (ties to the lower index). A list star's own detection sits at distance about 0 and would be read as a
/// neighbour of the star itself, so with `blend_refusal` on every list star would be refused as a `Blend`; `radius` is
/// larger than the kept tolerance because a mosaic's detection is a `2s + 0.5` centroid of the green-sum plane, up to
/// about a pixel off the star's true position, which the list position (a catalogue's) is not. Only the nearest is
/// taken: a second detection near a list star is a neighbour. The candidates come from `tree` (over `detected`, asked a
/// hair wider than `radius`), and the distance and the tie are judged on the candidates as ever: `hypot`, `<= radius`,
/// the lower index, so the answer does not depend on the tree's own rounding of the boundary.
fn own_detections(
    detected: &[Seed],
    tree: &KdTree2,
    listed: &[Seed],
    radius: f64,
) -> Vec<bool> {
    let mut own = vec![false; detected.len()];
    let wide = radius * (1.0 + 1e-9) + 1e-9;
    for l in listed {
        let mut nearest: Option<(f64, usize)> = None;
        for (j, _) in tree.within(l.x, l.y, wide) {
            let d = (detected[j].x - l.x).hypot(detected[j].y - l.y);
            if d <= radius && nearest.map_or(true, |(best, bj)| d < best || (d == best && j < bj)) {
                nearest = Some((d, j));
            }
        }
        if let Some((_, j)) = nearest {
            own[j] = true;
        }
    }
    own
}

/// A list fit that was not made, or was refused: NaN numbers, both flags false.
pub(crate) fn refused_fit() -> ListFit {
    ListFit {
        x: f64::NAN,
        y: f64::NAN,
        model_flux: f64::NAN,
        snr: f64::NAN,
        flags: StarFlags(0),
        accepted: false,
        kept: false,
        centre_offset_px: f64::NAN,
    }
}

/// One fit pass: a row per requested position, in request order.
pub(crate) struct Pass {
    pub(crate) fits: Vec<ListFit>,
    pub(crate) refused: RefusalCounts,
    pub(crate) ring_fwhm: Option<f64>,
    /// The fit call's own aperture correction and the stars it was measured over: `Some(1.0)` over 0 on every pass (the
    /// given unit correction; no growth curve ran). Read by one test only.
    #[cfg(test)]
    pub(crate) aperture_correction: Option<f64>,
    #[cfg(test)]
    pub(crate) aperture_stars: usize,
}

/// What both passes need of the plane and the detected seeds before the list is known: the working copy, its saturation
/// plane (a scan of every pixel) and a tree over the detected seeds. None of it reads the request, so
/// [`measure_raw_frame_and_list`](super::analysis_mode::measure_raw_frame_and_list) builds it on the pool while the
/// caller's callback builds the list.
pub(crate) struct Prebuilt<'a> {
    working: Cow<'a, [f32]>,
    saturation: SaturationPlane,
    tree: KdTree2,
}

impl<'a> Prebuilt<'a> {
    /// `None` when the plane or its mask does not match its geometry (nothing is scanned; [`PassContext::new`] says
    /// why when the list is measured).
    pub(crate) fn new(
        input: &'a ListPlane<'a>,
        pm: &PlaneMeasurement,
        detected: &[Seed],
    ) -> Option<Self> {
        let (w, h) = (input.w, input.h);
        if w == 0
            || h == 0
            || input.plane.len() != w * h
            || input.mask.is_some_and(|m| m.len() != w * h)
        {
            return None;
        }
        let (working, (saturation, tree)) = {
            let working = domain::to_working_units(input.plane, &input.domain);
            let both = rayon::join(
                || {
                    saturation_plane(
                        &working,
                        input.mask,
                        &pm.noise,
                        &pm.sky,
                        &input.domain,
                        false,
                    )
                },
                || {
                    let points: Vec<(f64, f64)> = detected.iter().map(|s| (s.x, s.y)).collect();
                    KdTree2::build(&points)
                },
            );
            (working, both)
        };
        Some(Prebuilt {
            working,
            saturation,
            tree,
        })
    }
}

/// What both passes share: the validated request's tolerance, the working copy of the plane and its saturation plane.
pub(crate) struct PassContext<'a> {
    input: &'a ListPlane<'a>,
    pm: &'a PlaneMeasurement,
    detected: &'a [Seed],
    tol: f64,
    working: Cow<'a, [f32]>,
    saturation: SaturationPlane,
    tree: KdTree2,
    pool: Option<&'a Arc<rayon::ThreadPool>>,
}

impl<'a> PassContext<'a> {
    /// `None` (always with a `warn!` or `debug!` that says why) when nothing can be measured: an empty plane, no usable
    /// tolerance, a plane or mask that does not match its geometry. `pre` is [`Prebuilt::new`]'s (built here, on this
    /// thread's pool, when the caller did not).
    pub(crate) fn new(
        input: &'a ListPlane<'a>,
        pm: &'a PlaneMeasurement,
        detected: &'a [Seed],
        tol: f64,
        pool: Option<&'a Arc<rayon::ThreadPool>>,
        pre: Option<Prebuilt<'a>>,
    ) -> Option<Self> {
        let (w, h) = (input.w, input.h);
        if w == 0 || h == 0 {
            tracing::debug!(
                width = w,
                height = h,
                "list: an empty plane; nothing measured"
            );
            return None;
        }
        if !(tol.is_finite() && tol > 0.0) {
            tracing::warn!(
                max_centre_offset_px = tol,
                "list: the centre tolerance must be positive; nothing measured"
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
        let Prebuilt {
            working,
            saturation,
            tree,
        } = match pre {
            Some(p) => p,
            None => install_in_span(pool, || Prebuilt::new(input, pm, detected))?,
        };
        Some(PassContext {
            input,
            pm,
            detected,
            tol,
            working,
            saturation,
            tree,
            pool,
        })
    }

    /// The plane's near-clip floor in working units: `bg + clip_floor·range` of the saturation plane, `None` when the
    /// plane shows no clipped evidence.
    pub(crate) fn clip_level(&self) -> Option<f64> {
        self.pm
            .clip_floor
            .map(|r| self.saturation.bg + r * self.saturation.range)
    }

    /// Fit `positions` (seeds, in request order): a position off the plane or not finite keeps its refused row. Kept:
    /// accepted, not saturated (the verdict pass's saturated and near-clip classes), centred within the tolerance of its
    /// seed.
    pub(crate) fn fit(&self, positions: &[(f64, f64)]) -> Pass {
        let (w, h) = (self.input.w, self.input.h);
        let mask = self.input.mask;
        let working: &[f32] = &self.working;
        let pm = self.pm;
        let tol = self.tol;
        let saturation = self.saturation;
        let mut fits: Vec<ListFit> = (0..positions.len()).map(|_| refused_fit()).collect();

        // The seeds: a position on the plane (a refused one keeps its NaN row), the working plane at its nearest present
        // pixel less the sky median as its peak and flux, which the verdict pass reads for its plateau tests.
        let mut seeds: Vec<Seed> = Vec::with_capacity(positions.len());
        let mut request_of: Vec<usize> = Vec::with_capacity(positions.len());
        for (r, &(x, y)) in positions.iter().enumerate() {
            if !(x.is_finite()
                && y.is_finite()
                && x >= 0.0
                && y >= 0.0
                && x <= (w - 1) as f64
                && y <= (h - 1) as f64)
            {
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
        let outside = positions.len() - seeds.len();
        if outside > 0 {
            tracing::debug!(
                count = outside,
                width = w,
                height = h,
                "list: positions outside the plane or not finite; refused without a fit"
            );
        }
        if seeds.is_empty() {
            return Pass {
                fits,
                refused: RefusalCounts::default(),
                ring_fwhm: None,
                #[cfg(test)]
                aperture_correction: None,
                #[cfg(test)]
                aperture_stars: 0,
            };
        }

        // The ring sky reads every detected seed except a list star's own detection (the nearest one within
        // max(tolerance, 1 F_field), R2), then the list.
        let own_radius = if pm.field_width.fwhm_px.is_finite() {
            tol.max(pm.field_width.fwhm_px)
        } else {
            tol
        };
        let own = own_detections(self.detected, &self.tree, &seeds, own_radius);
        let mut neighbour_seeds: Vec<Seed> = self
            .detected
            .iter()
            .zip(&own)
            .filter(|(_, &o)| !o)
            .map(|(d, _)| *d)
            .collect();
        let first = neighbour_seeds.len();
        neighbour_seeds.extend_from_slice(&seeds);
        let of_fit: Vec<usize> = (first..first + seeds.len()).collect();

        let fit_params = FitParams::default();
        let (bg, range) = (saturation.bg, saturation.range);
        let passes = prepare_passes(
            working,
            w,
            h,
            &fit_params,
            mask,
            self.pool,
            &seeds,
            pm.field_width,
            saturation,
            // R3: the near-clip class is the plane's own floor, `(floor - bg) / range` in range units on the measurement.
            Some((pm.clip_floor.map(|r| bg + r * range), pm.clip_evidence)),
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
                    // The list measures no `c`: the unit correction is given, so no growth curve runs (the ring sky and
                    // the SNR pass still do; the caller reads the SNR).
                    flux: FluxCorrection::Given(1.0),
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
                    split: Some(NeighbourSplit {
                        detected: self.detected,
                        tree: &self.tree,
                        own: &own,
                    }),
                    blend_refusal: true,
                    report_aperture: false,
                    ..Default::default()
                },
                None,
                false,
                false,
                false,
                &[],
            )
        };
        let outcome = install_in_span(self.pool, fit_call);
        for (f, &fit_seed) in outcome.fits.iter().zip(&outcome.fit_seeds) {
            let r = request_of[fit_seed];
            let (px, py) = positions[r];
            let offset = (f.x - px).hypot(f.y - py);
            fits[r] = ListFit {
                x: f.x,
                y: f.y,
                model_flux: photometry::model_flux_in_fwtm(f),
                snr: f.snr,
                flags: f.flags,
                accepted: true,
                kept: !f.flags.saturated() && offset.is_finite() && offset <= tol,
                centre_offset_px: offset,
            };
        }
        Pass {
            fits,
            refused: outcome.refused,
            ring_fwhm: outcome.ring_fwhm,
            #[cfg(test)]
            aperture_correction: outcome.aperture_correction,
            #[cfg(test)]
            aperture_stars: outcome.aperture_stars,
        }
    }

    /// The aperture sums of `placement`, one per aperture position, in parallel on the pool.
    pub(crate) fn apertures(&self, placement: &ListPlacement) -> Vec<ListAperture> {
        let (w, h) = (self.input.w, self.input.h);
        let working: &[f32] = &self.working;
        let mask = self.input.mask;
        // The positions are in the list's key order, scattered over the plane; each sum reads a window of rows. Taken in
        // row order, the stars a worker meets in a row run share those rows' cache lines; every position is summed on its
        // own, so the order it is taken in moves no number, and the sums are put back in the request's order.
        let positions = &placement.aperture_positions;
        let mut order: Vec<usize> = (0..positions.len()).collect();
        order.sort_unstable_by(|&a, &b| {
            positions[a].1.total_cmp(&positions[b].1).then(positions[a].0.total_cmp(&positions[b].0)).then(a.cmp(&b))
        });
        let body = || {
            order
                .par_iter()
                .map_init(photometry::ApertureScratch::default, |scratch, &i| {
                    (
                        i,
                        photometry::aperture_sums_in(
                            scratch,
                            working,
                            mask,
                            w,
                            h,
                            positions[i],
                            placement.radius_px,
                            &placement.extra_radii_px,
                            placement.annulus_start_px,
                            placement.annulus_area_px,
                            placement.peak_radius_px,
                        ),
                    )
                })
                .collect::<Vec<(usize, ListAperture)>>()
        };
        let mut slots: Vec<Option<ListAperture>> = (0..positions.len()).map(|_| None).collect();
        for (i, sums) in install_in_span(self.pool, body) {
            slots[i] = Some(sums);
        }
        slots
            .into_iter()
            .map(|sums| sums.expect("every position was summed once"))
            .collect()
    }
}

/// Fit `req`'s positions on `input` (`pm` is the plane's own measurement, `detected` its detected seeds in native
/// pixels: the matched filter's on a mono plane, the green-sum filter's on a mosaic), hand the fits to `place`, fit the
/// placement's refit positions and sum its aperture positions. `None` (always with a `warn!` or `debug!` that says why)
/// when the request cannot be measured: no position, no usable tolerance, a plane or mask that does not match its
/// geometry. A `place` that returns `None`, or an unusable placement (logged), ends the list with the first fits alone.
pub(crate) fn measure_list<P>(
    input: &ListPlane<'_>,
    pm: &PlaneMeasurement,
    detected: &[Seed],
    req: &ListRequest,
    place: P,
    pool: Option<&Arc<rayon::ThreadPool>>,
    pre: Option<Prebuilt<'_>>,
) -> Option<ListMeasurement>
where
    P: FnOnce(&[ListFit]) -> Option<ListPlacement>,
{
    if req.fit_positions.is_empty() {
        tracing::debug!("list: no position asked; nothing measured");
        return None;
    }
    let ctx = PassContext::new(input, pm, detected, req.max_centre_offset_px, pool, pre)?;
    let clip_level = ctx.clip_level();
    let mut timings = ListTimings::default();

    let t = Instant::now();
    let first = ctx.fit(&req.fit_positions);
    timings.fit_ms = t.elapsed().as_millis() as u64;
    let mut out = ListMeasurement {
        fits: first.fits,
        refits: Vec::new(),
        apertures: Vec::new(),
        clip_level,
        beta: pm.beta,
        ring_fwhm: first.ring_fwhm,
        refused: first.refused,
        refit_refused: RefusalCounts::default(),
        timings,
        refusal: None,
    };
    let Some(placement) = place(&out.fits) else {
        tracing::debug!(
            count = req.fit_positions.len(),
            "list: no placement; the first fits alone"
        );
        return Some(out);
    };
    if !placement.is_usable() {
        tracing::warn!(
            radius_px = placement.radius_px,
            extra_radii = placement.extra_radii_px.len(),
            annulus_start_px = placement.annulus_start_px,
            annulus_area_px = placement.annulus_area_px,
            peak_radius_px = placement.peak_radius_px,
            "list: the placement's radii and annulus must be positive and finite; the first fits alone"
        );
        return Some(out);
    }

    let t = Instant::now();
    let second = ctx.fit(&placement.refit_positions);
    out.timings.refit_ms = t.elapsed().as_millis() as u64;
    out.refits = second.fits;
    out.refit_refused = second.refused;

    let t = Instant::now();
    out.apertures = ctx.apertures(&placement);
    out.timings.aperture_ms = t.elapsed().as_millis() as u64;

    tracing::debug!(
        count = req.fit_positions.len(),
        kept = out.fits.iter().filter(|f| f.kept).count(),
        refits = out.refits.iter().filter(|f| f.kept).count(),
        apertures = out.apertures.iter().filter(|a| a.on_frame).count(),
        aperture_ms = out.timings.aperture_ms,
        "list measured"
    );
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// The own-detection rule as the cell grid wrote it (cells of side `radius`, the 3 x 3 around a list star), kept as
    /// the reference of the tree version.
    fn own_detections_grid(detected: &[Seed], listed: &[Seed], radius: f64) -> Vec<bool> {
        let cell = |x: f64, y: f64| ((x / radius).floor() as i64, (y / radius).floor() as i64);
        let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (j, d) in detected.iter().enumerate() {
            grid.entry(cell(d.x, d.y)).or_default().push(j);
        }
        let mut own = vec![false; detected.len()];
        for l in listed {
            let (cx, cy) = cell(l.x, l.y);
            let mut nearest: Option<(f64, usize)> = None;
            for gx in cx - 1..=cx + 1 {
                for gy in cy - 1..=cy + 1 {
                    for &j in grid.get(&(gx, gy)).into_iter().flatten() {
                        let d = (detected[j].x - l.x).hypot(detected[j].y - l.y);
                        if d <= radius
                            && nearest.map_or(true, |(best, bj)| d < best || (d == best && j < bj))
                        {
                            nearest = Some((d, j));
                        }
                    }
                }
            }
            if let Some((_, j)) = nearest {
                own[j] = true;
            }
        }
        own
    }

    fn seed(x: f64, y: f64) -> Seed {
        Seed {
            x,
            y,
            peak: 1.0,
            flux: 1.0,
        }
    }

    /// The tree finds the same own detections as the cell grid: on a half-pixel lattice (many equal distances, points
    /// exactly on the radius), with repeated points and a list star on top of several detections.
    #[test]
    fn the_own_detections_of_the_tree_are_the_grids() {
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for (radius, step) in [(2.0, 0.5), (3.0, 0.5), (1.0, 0.25), (2.5, 1.0)] {
            let lattice = |r: u64| (r % 80) as f64 * step;
            let detected: Vec<Seed> = (0..4000)
                .map(|_| {
                    let r = next();
                    seed(lattice(r), lattice(r >> 20))
                })
                .collect();
            let listed: Vec<Seed> = (0..600)
                .map(|_| {
                    let r = next();
                    seed(lattice(r) + 0.25 * ((r >> 40) % 3) as f64, lattice(r >> 20))
                })
                .collect();
            let points: Vec<(f64, f64)> = detected.iter().map(|s| (s.x, s.y)).collect();
            let tree = KdTree2::build(&points);
            let got = own_detections(&detected, &tree, &listed, radius);
            let want = own_detections_grid(&detected, &listed, radius);
            assert_eq!(got, want, "radius {radius}, step {step}");
            assert!(want.iter().any(|&o| o) && want.iter().any(|&o| !o));
        }
    }
}
