/// PSF metrics: windowed moments, optional 2D Gaussian fit, HFR.

use super::detection::DetectedStar;
use super::fitting::{fit_gaussian_2d, PixelSample};

/// Fully measured star with all metrics.
pub(crate) struct MeasuredStar {
    pub x: f32,
    pub y: f32,
    pub peak: f32,
    pub flux: f32,
    pub fwhm_x: f32,
    pub fwhm_y: f32,
    pub fwhm: f32,
    pub eccentricity: f32,
    pub hfr: f32,
    /// Per-star SNR — set later by snr::compute_star_snr.
    pub snr: f32,
    /// PSF position angle in radians, counter-clockwise from +X axis.
    pub theta: f32,
}

const FWHM_FACTOR: f32 = 2.3548;

/// Measure PSF metrics for all detected stars.
///
/// `data`: single-channel f32 image.
/// `stars`: detected star candidates from detection stage.
/// `background`: global background level.
/// `bg_map`: optional per-pixel background map.
/// `use_gaussian_fit`: if true, use full 2D Gaussian fit instead of windowed moments.
pub(crate) fn measure_stars(
    data: &[f32],
    width: usize,
    height: usize,
    stars: &[DetectedStar],
    background: f32,
    bg_map: Option<&[f32]>,
    use_gaussian_fit: bool,
    green_mask: Option<&[bool]>,
) -> Vec<MeasuredStar> {
    use rayon::prelude::*;

    stars
        .par_iter()
        .filter_map(|star| {
            measure_single_star(data, width, height, star, background, bg_map, use_gaussian_fit, green_mask)
        })
        .collect()
}

fn measure_single_star(
    data: &[f32],
    width: usize,
    height: usize,
    star: &DetectedStar,
    background: f32,
    bg_map: Option<&[f32]>,
    use_gaussian_fit: bool,
    green_mask: Option<&[bool]>,
) -> Option<MeasuredStar> {
    // Estimate sigma from the star's area: area ≈ π*(2σ)² for a Gaussian at low_threshold
    let estimated_sigma = (star.area as f32 / std::f32::consts::PI).sqrt() * 0.5;
    let estimated_sigma = estimated_sigma.max(1.0).min(20.0);

    // Stamp radius: capture >99.9% of Gaussian flux
    let stamp_radius = (4.0 * estimated_sigma).ceil() as usize;
    let stamp_radius = stamp_radius.max(8).min(50);

    let cx = star.x.round() as i32;
    let cy = star.y.round() as i32;

    // Check bounds
    if cx - stamp_radius as i32 <= 0
        || cy - stamp_radius as i32 <= 0
        || cx + stamp_radius as i32 >= width as i32 - 1
        || cy + stamp_radius as i32 >= height as i32 - 1
    {
        return None;
    }

    let x0 = (cx - stamp_radius as i32) as usize;
    let y0 = (cy - stamp_radius as i32) as usize;
    let x1 = (cx + stamp_radius as i32) as usize;
    let y1 = (cy + stamp_radius as i32) as usize;

    // Extract background-subtracted stamp
    let stamp_w = x1 - x0 + 1;
    let stamp_h = y1 - y0 + 1;
    let mut stamp = Vec::with_capacity(stamp_w * stamp_h);
    for sy in y0..=y1 {
        for sx in x0..=x1 {
            let bg = bg_map.map_or(background, |m| m[sy * width + sx]);
            stamp.push(data[sy * width + sx] - bg);
        }
    }

    // Centroid relative to stamp origin
    let rel_cx = star.x - x0 as f32;
    let rel_cy = star.y - y0 as f32;

    // Build stamp-local green mask: true at real green CFA positions
    let stamp_mask: Option<Vec<bool>> = green_mask.map(|gm| {
        let mut mask = Vec::with_capacity(stamp_w * stamp_h);
        for sy in y0..=y1 {
            for sx in x0..=x1 {
                mask.push(gm[sy * width + sx]);
            }
        }
        mask
    });
    let stamp_mask_ref = stamp_mask.as_deref();

    // Robust sigma estimate used by all paths (HFR, moments, Gaussian fit)
    // (uses bilinear interpolation — fine to use all pixels for initialization)
    let robust_sigma = estimate_sigma_halfmax(&stamp, stamp_w, rel_cx, rel_cy);
    let hfr_radius = 5.0_f32.max(4.0 * robust_sigma);

    // Compute HFR within the star vicinity (not the entire stamp)
    let hfr = compute_hfr(&stamp, stamp_w, rel_cx, rel_cy, hfr_radius, stamp_mask_ref);

    if use_gaussian_fit {
        // Full 2D Gaussian fit
        measure_with_gaussian_fit(&stamp, stamp_w, stamp_h, rel_cx, rel_cy, robust_sigma, star, hfr, x0, y0, stamp_mask_ref)
    } else {
        // Windowed moments
        measure_with_moments(&stamp, stamp_w, stamp_h, rel_cx, rel_cy, robust_sigma, star, hfr, x0, y0, stamp_mask_ref)
    }
}

/// Compute Half-Flux Radius within a given radius of the centroid.
/// HFR = Σ(w_i * d_i) / Σ(w_i) where w_i = max(0, pixel), d_i = distance from centroid.
/// When `stamp_mask` is `Some`, only pixels at green CFA positions are used.
fn compute_hfr(stamp: &[f32], stamp_w: usize, cx: f32, cy: f32, radius: f32, stamp_mask: Option<&[bool]>) -> f32 {
    let stamp_h = stamp.len() / stamp_w;
    let radius_sq = (radius as f64) * (radius as f64);
    let mut sum_wd = 0.0_f64;
    let mut sum_w = 0.0_f64;

    for sy in 0..stamp_h {
        for sx in 0..stamp_w {
            if let Some(mask) = stamp_mask {
                if !mask[sy * stamp_w + sx] {
                    continue;
                }
            }
            let dx = sx as f64 - cx as f64;
            let dy = sy as f64 - cy as f64;
            let dist_sq = dx * dx + dy * dy;
            if dist_sq > radius_sq {
                continue;
            }
            let val = stamp[sy * stamp_w + sx].max(0.0) as f64;
            if val > 0.0 {
                sum_wd += val * dist_sq.sqrt();
                sum_w += val;
            }
        }
    }

    if sum_w > 0.0 {
        (sum_wd / sum_w) as f32
    } else {
        0.0
    }
}

/// Measure using intensity-weighted second-order moments.
/// Uses I²-weighting for centroid refinement, then I-weighted second moments
/// for FWHM/eccentricity. Two iterations to converge centroid.
fn measure_with_moments(
    stamp: &[f32],
    stamp_w: usize,
    stamp_h: usize,
    init_cx: f32,
    init_cy: f32,
    robust_sigma: f32,
    star: &DetectedStar,
    hfr: f32,
    x0: usize,
    y0: usize,
    stamp_mask: Option<&[bool]>,
) -> Option<MeasuredStar> {
    // Limit moments to star vicinity using pre-computed robust sigma
    let fitting_radius = 5.0_f64.max(4.0 * robust_sigma as f64);
    let fitting_radius_sq = fitting_radius * fitting_radius;

    let mut cx = init_cx as f64;
    let mut cy = init_cy as f64;

    let mut m_xx = 0.0_f64;
    let mut m_yy = 0.0_f64;
    let mut m_xy = 0.0_f64;

    // Iterate 2× to converge centroid
    for _ in 0..2 {
        m_xx = 0.0;
        m_yy = 0.0;
        m_xy = 0.0;
        let mut sum_w = 0.0_f64;
        let mut sum_wx = 0.0_f64;
        let mut sum_wy = 0.0_f64;

        for sy in 0..stamp_h {
            for sx in 0..stamp_w {
                // Skip non-green pixels in OSC
                if let Some(mask) = stamp_mask {
                    if !mask[sy * stamp_w + sx] {
                        continue;
                    }
                }
                let dx = sx as f64 - cx;
                let dy = sy as f64 - cy;
                if dx * dx + dy * dy > fitting_radius_sq {
                    continue;
                }
                let val = stamp[sy * stamp_w + sx].max(0.0) as f64;
                if val <= 0.0 {
                    continue;
                }

                // I-weighting for moments
                sum_w += val;
                sum_wx += val * sx as f64;
                sum_wy += val * sy as f64;
                m_xx += val * dx * dx;
                m_yy += val * dy * dy;
                m_xy += val * dx * dy;
            }
        }

        if sum_w < 1e-10 {
            return None;
        }

        // Update centroid
        cx = sum_wx / sum_w;
        cy = sum_wy / sum_w;

        // Normalize moments
        m_xx /= sum_w;
        m_yy /= sum_w;
        m_xy /= sum_w;
    }

    // Eigenvalue decomposition of 2×2 covariance matrix
    let trace = m_xx + m_yy;
    let det = m_xx * m_yy - m_xy * m_xy;
    let disc = (trace * trace - 4.0 * det).max(0.0);

    let lambda1 = (trace + disc.sqrt()) * 0.5; // Major axis
    let lambda2 = (trace - disc.sqrt()) * 0.5; // Minor axis

    if lambda1 <= 0.0 || lambda2 <= 0.0 {
        return None;
    }

    let fwhm_x = FWHM_FACTOR as f64 * lambda1.sqrt();
    let fwhm_y = FWHM_FACTOR as f64 * lambda2.sqrt();
    let fwhm = (fwhm_x * fwhm_y).sqrt(); // Geometric mean

    let eccentricity = (1.0 - lambda2 / lambda1).max(0.0).sqrt();

    // Position angle of major axis from eigenvector direction
    let theta = 0.5 * (2.0 * m_xy).atan2(m_xx - m_yy);

    Some(MeasuredStar {
        x: (cx + x0 as f64) as f32,
        y: (cy + y0 as f64) as f32,
        peak: star.peak,
        flux: star.flux,
        fwhm_x: fwhm_x as f32,
        fwhm_y: fwhm_y as f32,
        fwhm: fwhm as f32,
        eccentricity: eccentricity as f32,
        hfr,
        snr: 0.0,
        theta: theta as f32,
    })
}

/// Compute ellipse parameters from intensity-weighted second-order moments.
/// Returns `(theta, sigma_major, sigma_minor)` where theta is the position angle
/// of the major axis in radians, and sigmas are sqrt(eigenvalue) of the 2×2
/// covariance matrix.
fn moments_ellipse(stamp: &[f32], stamp_w: usize, cx: f32, cy: f32) -> (f64, f64, f64) {
    let stamp_h = stamp.len() / stamp_w;
    let mut m_xx = 0.0_f64;
    let mut m_yy = 0.0_f64;
    let mut m_xy = 0.0_f64;
    let mut sum_w = 0.0_f64;

    for sy in 0..stamp_h {
        for sx in 0..stamp_w {
            let val = stamp[sy * stamp_w + sx].max(0.0) as f64;
            if val <= 0.0 {
                continue;
            }
            let dx = sx as f64 - cx as f64;
            let dy = sy as f64 - cy as f64;
            m_xx += val * dx * dx;
            m_yy += val * dy * dy;
            m_xy += val * dx * dy;
            sum_w += val;
        }
    }

    if sum_w < 1e-10 {
        return (0.0, 1.0, 1.0);
    }

    m_xx /= sum_w;
    m_yy /= sum_w;
    m_xy /= sum_w;

    let theta = 0.5 * (2.0 * m_xy).atan2(m_xx - m_yy);

    let trace = m_xx + m_yy;
    let det = m_xx * m_yy - m_xy * m_xy;
    let disc = (trace * trace - 4.0 * det).max(0.0).sqrt();

    let lambda_major = ((trace + disc) * 0.5).max(0.01);
    let lambda_minor = ((trace - disc) * 0.5).max(0.01);

    (theta, lambda_major.sqrt(), lambda_minor.sqrt())
}

/// Measure using full 2D elliptical Gaussian fit.
fn measure_with_gaussian_fit(
    stamp: &[f32],
    stamp_w: usize,
    stamp_h: usize,
    init_cx: f32,
    init_cy: f32,
    init_sigma: f32,
    star: &DetectedStar,
    hfr: f32,
    x0: usize,
    y0: usize,
    stamp_mask: Option<&[bool]>,
) -> Option<MeasuredStar> {

    // Limit fitting radius to exclude distant nebula pixels
    let fitting_radius = 5.0_f64.max(4.0 * init_sigma as f64);
    let fitting_radius_sq = fitting_radius * fitting_radius;
    let cx64 = init_cx as f64;
    let cy64 = init_cy as f64;

    // Prepare pixel samples — only within fitting radius, no clamping.
    // For OSC, skip non-green pixels (interpolated) to avoid PSF broadening.
    let pixels: Vec<PixelSample> = (0..stamp_h)
        .flat_map(|sy| {
            (0..stamp_w).filter_map(move |sx| {
                // Skip non-green pixels in OSC
                if let Some(mask) = stamp_mask {
                    if !mask[sy * stamp_w + sx] {
                        return None;
                    }
                }
                let dx = sx as f64 - cx64;
                let dy = sy as f64 - cy64;
                if dx * dx + dy * dy <= fitting_radius_sq {
                    Some(PixelSample {
                        x: sx as f64,
                        y: sy as f64,
                        value: stamp[sy * stamp_w + sx] as f64,
                    })
                } else {
                    None
                }
            })
        })
        .collect();

    // Background estimate from annulus just outside the fitting radius.
    // Using stamp edges would underestimate in nebula regions where the local
    // background gradient is significant over the large stamp.
    let bg_inner_sq = fitting_radius_sq;
    let bg_outer = fitting_radius + 3.0;
    let bg_outer_sq = bg_outer * bg_outer;
    let mut annulus_bg_vals = Vec::new();
    for sy in 0..stamp_h {
        for sx in 0..stamp_w {
            let dx = sx as f64 - cx64;
            let dy = sy as f64 - cy64;
            let r_sq = dx * dx + dy * dy;
            if r_sq > bg_inner_sq && r_sq <= bg_outer_sq {
                annulus_bg_vals.push(stamp[sy * stamp_w + sx] as f64);
            }
        }
    }
    annulus_bg_vals.sort_by(|a, b| a.total_cmp(b));
    let init_bg = if annulus_bg_vals.is_empty() {
        // Fallback to edge pixels if annulus is empty (small stamp)
        let mut edge_vals = Vec::new();
        for sx in 0..stamp_w {
            edge_vals.push(stamp[sx] as f64);
            edge_vals.push(stamp[(stamp_h - 1) * stamp_w + sx] as f64);
        }
        for sy in 1..stamp_h - 1 {
            edge_vals.push(stamp[sy * stamp_w] as f64);
            edge_vals.push(stamp[sy * stamp_w + stamp_w - 1] as f64);
        }
        edge_vals.sort_by(|a, b| a.total_cmp(b));
        if edge_vals.is_empty() { 0.0 } else { edge_vals[edge_vals.len() / 2] }
    } else {
        annulus_bg_vals[annulus_bg_vals.len() / 2]
    };

    // Use moments to get initial theta and per-axis sigmas, then scale
    // to match the halfmax geometric mean (moments give good axis ratio +
    // angle, halfmax gives good absolute scale).
    let (mom_theta, mom_major, mom_minor) = moments_ellipse(stamp, stamp_w, init_cx, init_cy);
    let scale = init_sigma as f64 / (mom_major * mom_minor).sqrt().max(0.1);
    let init_sx = (mom_major * scale).max(0.5);
    let init_sy = (mom_minor * scale).max(0.5);

    let result = fit_gaussian_2d(
        &pixels,
        init_bg,
        star.peak as f64 - init_bg,
        init_cx as f64,
        init_cy as f64,
        init_sx,
        init_sy,
        mom_theta,
    )?;

    // Reject divergent fits: fitted sigma should not exceed 3x the initial estimate
    let max_sigma = (init_sigma as f64) * 3.0;
    if result.sigma_x > max_sigma || result.sigma_y > max_sigma
        || result.sigma_x <= 0.0 || result.sigma_y <= 0.0
    {
        return None;
    }

    let fwhm_x = FWHM_FACTOR as f64 * result.sigma_x;
    let fwhm_y = FWHM_FACTOR as f64 * result.sigma_y;
    let fwhm = (fwhm_x * fwhm_y).sqrt();

    let min_s = result.sigma_x.min(result.sigma_y);
    let max_s = result.sigma_x.max(result.sigma_y);
    let eccentricity = (1.0 - (min_s * min_s) / (max_s * max_s)).max(0.0).sqrt();

    Some(MeasuredStar {
        x: (result.x0 + x0 as f64) as f32,
        y: (result.y0 + y0 as f64) as f32,
        peak: star.peak,
        flux: star.flux,
        fwhm_x: fwhm_x as f32,
        fwhm_y: fwhm_y as f32,
        fwhm: fwhm as f32,
        eccentricity: eccentricity as f32,
        hfr,
        snr: 0.0,
        theta: result.theta as f32,
    })
}

/// Robust sigma estimate by walking from centroid to half-max in 8 directions.
///
/// Subtracts the local edge-median background to remove nebula pedestal, then
/// bilinearly interpolates pixel values along each ray from the centroid outward.
/// Finds where the value drops below `local_bg + 0.5 * (peak - local_bg)` and
/// averages those distances. Converts half-max radius to sigma via `r / sqrt(2 * ln(2))`.
///
/// Falls back to capped second-moment estimate if fewer than 4 directions cross half-max.
fn estimate_sigma_halfmax(stamp: &[f32], stamp_w: usize, cx: f32, cy: f32) -> f32 {
    let stamp_h = stamp.len() / stamp_w;

    // Local background from annulus at r=4..8 pixels from centroid.
    // This captures the nebula pedestal close to the star, unlike stamp edges
    // which can be 50px away where the mesh background is accurate.
    let inner_r_sq = 4.0_f32 * 4.0;
    let outer_r_sq = 8.0_f32 * 8.0;
    let mut annulus_vals = Vec::new();
    for sy in 0..stamp_h {
        for sx in 0..stamp_w {
            let dx = sx as f32 - cx;
            let dy = sy as f32 - cy;
            let r_sq = dx * dx + dy * dy;
            if r_sq >= inner_r_sq && r_sq <= outer_r_sq {
                annulus_vals.push(stamp[sy * stamp_w + sx]);
            }
        }
    }
    annulus_vals.sort_by(|a, b| a.total_cmp(b));
    let local_bg = if annulus_vals.is_empty() {
        0.0
    } else {
        annulus_vals[annulus_vals.len() / 2]
    };

    // Bilinear interpolation within the stamp
    let interp = |x: f32, y: f32| -> f32 {
        let ix = x.floor() as i32;
        let iy = y.floor() as i32;
        if ix < 0 || iy < 0 || ix + 1 >= stamp_w as i32 || iy + 1 >= stamp_h as i32 {
            return 0.0;
        }
        let fx = x - ix as f32;
        let fy = y - iy as f32;
        let (ix, iy) = (ix as usize, iy as usize);
        let v00 = stamp[iy * stamp_w + ix];
        let v10 = stamp[iy * stamp_w + ix + 1];
        let v01 = stamp[(iy + 1) * stamp_w + ix];
        let v11 = stamp[(iy + 1) * stamp_w + ix + 1];
        (1.0 - fy) * ((1.0 - fx) * v00 + fx * v10) + fy * ((1.0 - fx) * v01 + fx * v11)
    };

    let peak = interp(cx, cy);
    let amplitude = peak - local_bg;
    if amplitude <= 0.0 {
        return 2.0;
    }
    // Half-max threshold above local background: bg + 0.5 * amplitude
    let half_max_threshold = local_bg + amplitude * 0.5;

    // 8 directions: N, NE, E, SE, S, SW, W, NW
    const DIRS: [(f32, f32); 8] = [
        (0.0, -1.0),
        (std::f32::consts::FRAC_1_SQRT_2, -std::f32::consts::FRAC_1_SQRT_2),
        (1.0, 0.0),
        (std::f32::consts::FRAC_1_SQRT_2, std::f32::consts::FRAC_1_SQRT_2),
        (0.0, 1.0),
        (-std::f32::consts::FRAC_1_SQRT_2, std::f32::consts::FRAC_1_SQRT_2),
        (-1.0, 0.0),
        (-std::f32::consts::FRAC_1_SQRT_2, -std::f32::consts::FRAC_1_SQRT_2),
    ];

    let max_walk = (stamp_w.max(stamp_h) / 2) as f32;
    let step = 0.5_f32; // sub-pixel stepping

    let mut half_max_distances = Vec::with_capacity(8);

    for &(ddx, ddy) in &DIRS {
        let mut prev_val = peak;
        let mut t = step;
        while t <= max_walk {
            let x = cx + ddx * t;
            let y = cy + ddy * t;
            let val = interp(x, y);
            if val <= half_max_threshold {
                // Linear interpolation between previous and current step for sub-step accuracy
                if prev_val > half_max_threshold {
                    let frac = (prev_val - half_max_threshold) / (prev_val - val);
                    half_max_distances.push(t - step + step * frac);
                }
                break;
            }
            prev_val = val;
            t += step;
        }
    }

    if half_max_distances.len() >= 4 {
        let mean_r: f32 = half_max_distances.iter().sum::<f32>() / half_max_distances.len() as f32;
        // sigma = half_max_radius / sqrt(2 * ln(2))
        let sigma = mean_r / (2.0_f32 * 2.0_f32.ln()).sqrt();
        sigma.max(0.5)
    } else {
        // Fallback: second-moment estimate, but capped
        let stamp_radius = (stamp_w.min(stamp_h) / 2) as f32;
        let fallback = estimate_sigma_from_stamp(stamp, stamp_w, cx, cy);
        fallback.min(stamp_radius / 4.0).max(0.5)
    }
}

/// Quick sigma estimate from stamp using second moments around centroid.
fn estimate_sigma_from_stamp(stamp: &[f32], stamp_w: usize, cx: f32, cy: f32) -> f32 {
    let stamp_h = stamp.len() / stamp_w;
    let mut sum_wr2 = 0.0_f64;
    let mut sum_w = 0.0_f64;

    for sy in 0..stamp_h {
        for sx in 0..stamp_w {
            let val = stamp[sy * stamp_w + sx].max(0.0) as f64;
            if val > 0.0 {
                let dx = sx as f64 - cx as f64;
                let dy = sy as f64 - cy as f64;
                sum_wr2 += val * (dx * dx + dy * dy);
                sum_w += val;
            }
        }
    }

    if sum_w > 0.0 {
        (sum_wr2 / sum_w).sqrt() as f32
    } else {
        2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_gaussian_stamp(
        size: usize,
        cx: f32,
        cy: f32,
        amp: f32,
        sigma_x: f32,
        sigma_y: f32,
        theta: f32,
    ) -> Vec<f32> {
        let mut stamp = vec![0.0_f32; size * size];
        let (ct, st) = (theta.cos(), theta.sin());
        let inv_sx2 = 1.0 / (2.0 * sigma_x * sigma_x);
        let inv_sy2 = 1.0 / (2.0 * sigma_y * sigma_y);

        for y in 0..size {
            for x in 0..size {
                let dx = x as f32 - cx;
                let dy = y as f32 - cy;
                let u = dx * ct + dy * st;
                let v = -dx * st + dy * ct;
                stamp[y * size + x] = amp * (-u * u * inv_sx2 - v * v * inv_sy2).exp();
            }
        }
        stamp
    }

    #[test]
    fn test_moments_isotropic() {
        let size = 31;
        let sigma = 3.0;
        let stamp = make_gaussian_stamp(size, 15.0, 15.0, 5000.0, sigma, sigma, 0.0);
        let star = DetectedStar {
            x: 15.0,
            y: 15.0,
            peak: 5000.0,
            flux: 50000.0,
            area: 50,
            theta: 0.0,
            eccentricity: 0.0,
        };

        let result =
            measure_with_moments(&stamp, size, size, 15.0, 15.0, sigma, &star, 3.0, 0, 0, None)
                .unwrap();

        let expected_fwhm = FWHM_FACTOR * sigma;
        assert!(
            (result.fwhm - expected_fwhm).abs() < 0.3,
            "FWHM {} expected ~{}",
            result.fwhm,
            expected_fwhm
        );
        assert!(
            result.eccentricity < 0.05,
            "Eccentricity {} should be ~0 for isotropic",
            result.eccentricity
        );
    }

    #[test]
    fn test_moments_elongated() {
        let size = 41;
        let sx = 2.0_f32;
        let sy = 5.0_f32;
        let stamp = make_gaussian_stamp(size, 20.0, 20.0, 5000.0, sx, sy, 0.0);
        let star = DetectedStar {
            x: 20.0,
            y: 20.0,
            peak: 5000.0,
            flux: 50000.0,
            area: 80,
            theta: 0.0,
            eccentricity: 0.0,
        };

        let result =
            measure_with_moments(&stamp, size, size, 20.0, 20.0, 3.5, &star, 4.0, 0, 0, None).unwrap();

        // Expected eccentricity: sqrt(1 - 4/25) = sqrt(0.84) ≈ 0.917
        let expected_ecc = (1.0 - (sx * sx) / (sy * sy)).sqrt();
        assert!(
            (result.eccentricity - expected_ecc).abs() < 0.05,
            "Eccentricity {} expected ~{}",
            result.eccentricity,
            expected_ecc
        );
    }

    #[test]
    fn test_hfr_gaussian() {
        let size = 31;
        let sigma = 3.0_f32;
        let stamp = make_gaussian_stamp(size, 15.0, 15.0, 5000.0, sigma, sigma, 0.0);
        let hfr = compute_hfr(&stamp, size, 15.0, 15.0, 15.0, None);

        // For Gaussian: HFR ≈ σ * √(π/2) ≈ 3.76 for σ=3
        let expected = sigma * (std::f32::consts::PI / 2.0).sqrt();
        assert!(
            (hfr - expected).abs() < 0.3,
            "HFR {} expected ~{}",
            hfr,
            expected
        );
    }

    #[test]
    fn test_gaussian_fit_vs_moments() {
        let size = 31;
        let sigma = 3.0;
        let stamp = make_gaussian_stamp(size, 15.0, 15.0, 5000.0, sigma, sigma, 0.0);
        let star = DetectedStar {
            x: 15.0,
            y: 15.0,
            peak: 5000.0,
            flux: 50000.0,
            area: 50,
            theta: 0.0,
            eccentricity: 0.0,
        };

        let moments_result =
            measure_with_moments(&stamp, size, size, 15.0, 15.0, sigma, &star, 3.0, 0, 0, None)
                .unwrap();
        let fit_result =
            measure_with_gaussian_fit(&stamp, size, size, 15.0, 15.0, sigma, &star, 3.0, 0, 0, None).unwrap();

        // Should agree within 10%
        let diff_pct = (moments_result.fwhm - fit_result.fwhm).abs() / fit_result.fwhm * 100.0;
        assert!(
            diff_pct < 10.0,
            "Moments FWHM {} vs Fit FWHM {}: diff {}%",
            moments_result.fwhm,
            fit_result.fwhm,
            diff_pct
        );
    }

    #[test]
    fn test_gaussian_fit_45deg_eccentricity() {
        // End-to-end: elongated star at 45° through measure_with_gaussian_fit.
        // sigma_x=2, sigma_y=5, theta=π/4 — the exact failure case.
        let size = 41;
        let sx = 2.0_f32;
        let sy = 5.0_f32;
        let theta = std::f32::consts::FRAC_PI_4;
        let stamp = make_gaussian_stamp(size, 20.0, 20.0, 5000.0, sx, sy, theta);
        let star = DetectedStar {
            x: 20.0,
            y: 20.0,
            peak: 5000.0,
            flux: 50000.0,
            area: 80,
            theta: 0.0,
            eccentricity: 0.0,
        };

        let init_sigma = (sx * sy).sqrt(); // geometric mean ≈ 3.16
        let result =
            measure_with_gaussian_fit(&stamp, size, size, 20.0, 20.0, init_sigma, &star, 4.0, 0, 0, None)
                .unwrap();

        // Expected eccentricity: sqrt(1 - 4/25) ≈ 0.917
        let expected_ecc = (1.0 - (sx * sx) / (sy * sy)).sqrt();
        assert!(
            (result.eccentricity - expected_ecc).abs() < 0.1,
            "eccentricity {} expected ~{}",
            result.eccentricity,
            expected_ecc
        );

        // Theta should be near ±π/4 (equivalent representations with swapped axes)
        let theta_diff = (result.theta.abs() - theta.abs()).abs();
        assert!(
            theta_diff < 0.15,
            "theta {} expected ~±{}",
            result.theta,
            theta
        );
    }
}
