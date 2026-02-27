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
) -> Vec<MeasuredStar> {
    use rayon::prelude::*;

    stars
        .par_iter()
        .filter_map(|star| {
            measure_single_star(data, width, height, star, background, bg_map, use_gaussian_fit)
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

    // Compute HFR (always, regardless of fitting method)
    let hfr = compute_hfr(&stamp, stamp_w, rel_cx, rel_cy);

    if use_gaussian_fit {
        // Full 2D Gaussian fit
        measure_with_gaussian_fit(&stamp, stamp_w, stamp_h, rel_cx, rel_cy, star, hfr, x0, y0)
    } else {
        // Windowed moments
        measure_with_moments(&stamp, stamp_w, stamp_h, rel_cx, rel_cy, estimated_sigma, star, hfr, x0, y0)
    }
}

/// Compute Half-Flux Radius.
/// HFR = Σ(w_i * d_i) / Σ(w_i) where w_i = max(0, pixel - bg), d_i = distance from centroid.
fn compute_hfr(stamp: &[f32], stamp_w: usize, cx: f32, cy: f32) -> f32 {
    let stamp_h = stamp.len() / stamp_w;
    let mut sum_wd = 0.0_f64;
    let mut sum_w = 0.0_f64;

    for sy in 0..stamp_h {
        for sx in 0..stamp_w {
            let val = stamp[sy * stamp_w + sx].max(0.0) as f64;
            if val > 0.0 {
                let dx = sx as f64 - cx as f64;
                let dy = sy as f64 - cy as f64;
                let dist = (dx * dx + dy * dy).sqrt();
                sum_wd += val * dist;
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
    _sigma_w: f32,
    star: &DetectedStar,
    hfr: f32,
    x0: usize,
    y0: usize,
) -> Option<MeasuredStar> {
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
                let val = stamp[sy * stamp_w + sx].max(0.0) as f64;
                if val <= 0.0 {
                    continue;
                }
                let dx = sx as f64 - cx;
                let dy = sy as f64 - cy;

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
    })
}

/// Measure using full 2D elliptical Gaussian fit.
fn measure_with_gaussian_fit(
    stamp: &[f32],
    stamp_w: usize,
    stamp_h: usize,
    init_cx: f32,
    init_cy: f32,
    star: &DetectedStar,
    hfr: f32,
    x0: usize,
    y0: usize,
) -> Option<MeasuredStar> {
    // Prepare pixel samples
    let pixels: Vec<PixelSample> = (0..stamp_h)
        .flat_map(|sy| {
            (0..stamp_w).map(move |sx| PixelSample {
                x: sx as f64,
                y: sy as f64,
                value: stamp[sy * stamp_w + sx].max(0.0) as f64,
            })
        })
        .collect();

    // Estimate initial sigma from the stamp data (use windowed moments for robustness)
    let init_sigma = estimate_sigma_from_stamp(stamp, stamp_w, init_cx, init_cy);

    // Edge pixel median as background estimate
    let mut edge_vals = Vec::new();
    for sx in 0..stamp_w {
        edge_vals.push(stamp[sx].max(0.0) as f64);
        edge_vals.push(stamp[(stamp_h - 1) * stamp_w + sx].max(0.0) as f64);
    }
    for sy in 1..stamp_h - 1 {
        edge_vals.push(stamp[sy * stamp_w].max(0.0) as f64);
        edge_vals.push(stamp[sy * stamp_w + stamp_w - 1].max(0.0) as f64);
    }
    edge_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let init_bg = if edge_vals.is_empty() {
        0.0
    } else {
        edge_vals[edge_vals.len() / 2]
    };

    let result = fit_gaussian_2d(
        &pixels,
        init_bg,
        star.peak as f64 - init_bg,
        init_cx as f64,
        init_cy as f64,
        init_sigma as f64,
    )?;

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
    })
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
        };

        let result =
            measure_with_moments(&stamp, size, size, 15.0, 15.0, sigma, &star, 3.0, 0, 0)
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
        };

        let result =
            measure_with_moments(&stamp, size, size, 20.0, 20.0, 3.5, &star, 4.0, 0, 0).unwrap();

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
        let hfr = compute_hfr(&stamp, size, 15.0, 15.0);

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
        };

        let moments_result =
            measure_with_moments(&stamp, size, size, 15.0, 15.0, sigma, &star, 3.0, 0, 0)
                .unwrap();
        let fit_result =
            measure_with_gaussian_fit(&stamp, size, size, 15.0, 15.0, &star, 3.0, 0, 0).unwrap();

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
}
