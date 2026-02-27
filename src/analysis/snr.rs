/// Per-star aperture photometry SNR and image-wide metrics.

use super::metrics::MeasuredStar;
use crate::processing::stretch::find_median;

/// Compute per-star SNR using aperture photometry with local sky annulus.
///
/// `data`: single-channel f32 image.
/// `stars`: measured stars (with centroid positions).
/// `median_fwhm`: used to scale aperture/annulus radii.
/// `background`: global background level.
/// `bg_map`: optional per-pixel background map.
pub(crate) fn compute_star_snr(
    data: &[f32],
    width: usize,
    height: usize,
    stars: &mut [MeasuredStar],
    median_fwhm: f32,
) {
    use rayon::prelude::*;

    let r_ap = (1.5 * median_fwhm).max(3.0);
    let r_inner = 3.0 * median_fwhm;
    let r_outer = 5.0 * median_fwhm;

    stars.par_iter_mut().for_each(|star| {
        star.snr = compute_single_snr(data, width, height, star.x, star.y, r_ap, r_inner, r_outer);
    });
}

fn compute_single_snr(
    data: &[f32],
    width: usize,
    height: usize,
    cx: f32,
    cy: f32,
    r_ap: f32,
    r_inner: f32,
    r_outer: f32,
) -> f32 {
    let r_outer_ceil = r_outer.ceil() as i32;
    let ix = cx.round() as i32;
    let iy = cy.round() as i32;

    // Check bounds
    if ix - r_outer_ceil < 0
        || iy - r_outer_ceil < 0
        || ix + r_outer_ceil >= width as i32
        || iy + r_outer_ceil >= height as i32
    {
        return 0.0;
    }

    let r_ap_sq = r_ap * r_ap;
    let r_inner_sq = r_inner * r_inner;
    let r_outer_sq = r_outer * r_outer;

    // Single pass: collect aperture pixel values and annulus pixels simultaneously
    let mut aperture_vals = Vec::new();
    let mut annulus_vals = Vec::new();

    let frac_x = cx - ix as f32;
    let frac_y = cy - iy as f32;

    for dy in -r_outer_ceil..=r_outer_ceil {
        for dx in -r_outer_ceil..=r_outer_ceil {
            let ddx = dx as f32 - frac_x;
            let ddy = dy as f32 - frac_y;
            let dist_sq = ddx * ddx + ddy * ddy;

            if dist_sq <= r_ap_sq {
                let px = (ix + dx) as usize;
                let py = (iy + dy) as usize;
                aperture_vals.push(data[py * width + px]);
            } else if dist_sq >= r_inner_sq && dist_sq <= r_outer_sq {
                let px = (ix + dx) as usize;
                let py = (iy + dy) as usize;
                annulus_vals.push(data[py * width + px]);
            }
        }
    }

    let n_ap = aperture_vals.len();
    if n_ap == 0 || annulus_vals.len() < 5 {
        return 0.0;
    }

    // Sigma-clip annulus (2 rounds, 3σ)
    let mut clipped = annulus_vals;
    let mut devs = Vec::with_capacity(clipped.len());
    for _ in 0..2 {
        if clipped.len() < 5 {
            break;
        }
        let med = find_median(&mut clipped);
        devs.clear();
        devs.extend(clipped.iter().map(|&v| (v - med).abs()));
        let mad = find_median(&mut devs);
        let sigma = 1.4826 * mad;
        if sigma < 1e-6 {
            break;
        }
        let lo = med - 3.0 * sigma;
        let hi = med + 3.0 * sigma;
        clipped.retain(|&v| v >= lo && v <= hi);
    }

    if clipped.len() < 5 {
        return 0.0;
    }

    let bg_local = find_median(&mut clipped);
    devs.clear();
    devs.extend(clipped.iter().map(|&v| (v - bg_local).abs()));
    let mad = find_median(&mut devs);
    let sigma_local = (1.4826 * mad) as f64;

    // Compute star flux from cached aperture values
    let mut star_flux = 0.0_f64;
    for &val in &aperture_vals {
        star_flux += (val - bg_local).max(0.0) as f64;
    }

    if star_flux <= 0.0 {
        return 0.0;
    }

    // SNR = F_star / sqrt(F_star + n_pix * sigma_local²)
    let noise_sq = star_flux + n_ap as f64 * sigma_local * sigma_local;
    if noise_sq <= 0.0 {
        return 0.0;
    }

    (star_flux / noise_sq.sqrt()) as f32
}

/// Compute image-wide SNR Weight (PixInsight-style).
/// SNRWeight = (MeanDeviation / noise)²
pub(crate) fn compute_snr_weight(data: &[f32], background: f32, noise: f32) -> f32 {
    if noise <= 0.0 || data.is_empty() {
        return 0.0;
    }

    // Subsample for performance
    let stride = (data.len() / 500_000).max(1);
    let mut sum_dev = 0.0_f64;
    let mut count = 0u64;

    for (i, &val) in data.iter().enumerate() {
        if i % stride != 0 {
            continue;
        }
        sum_dev += (val - background).abs() as f64;
        count += 1;
    }

    if count == 0 {
        return 0.0;
    }

    let mean_dev = (sum_dev / count as f64) as f32;
    let ratio = mean_dev / noise;
    ratio * ratio
}

/// Compute whole-image SNR in decibels: 20 × log10(mean / noise).
/// Matches PixInsight SNRViews methodology.
pub(crate) fn compute_snr_db(data: &[f32], noise: f32) -> f32 {
    if noise <= 0.0 || data.is_empty() {
        return 0.0;
    }

    let stride = (data.len() / 500_000).max(1);
    let mut sum = 0.0_f64;
    let mut count = 0u64;

    for (i, &val) in data.iter().enumerate() {
        if i % stride != 0 {
            continue;
        }
        sum += val as f64;
        count += 1;
    }

    if count == 0 {
        return 0.0;
    }

    let mean = (sum / count as f64) as f32;
    if mean <= 0.0 {
        return 0.0;
    }

    20.0 * (mean / noise).log10()
}

/// Compute PSF Signal = median(star_peaks) / noise.
pub(crate) fn compute_psf_signal(stars: &[MeasuredStar], noise: f32) -> f32 {
    if stars.is_empty() || noise <= 0.0 {
        return 0.0;
    }

    let mut peaks: Vec<f32> = stars.iter().map(|s| s.peak).collect();
    let median_peak = find_median(&mut peaks);
    median_peak / noise
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_snr_known_values() {
        // Create simple image with a single star
        let width = 100;
        let height = 100;
        let background = 1000.0_f32;
        let _noise = 50.0_f32;
        let mut data = vec![background; width * height];

        // Add a Gaussian star at center
        let cx = 50.0_f32;
        let cy = 50.0_f32;
        let sigma = 3.0_f32;
        let amp = 5000.0_f32;
        for y in 0..height {
            for x in 0..width {
                let dx = x as f32 - cx;
                let dy = y as f32 - cy;
                let r2 = dx * dx + dy * dy;
                data[y * width + x] += amp * (-r2 / (2.0 * sigma * sigma)).exp();
            }
        }

        let fwhm = 2.3548 * sigma;
        let snr = compute_single_snr(
            &data, width, height, cx, cy, 1.5 * fwhm, 3.0 * fwhm, 5.0 * fwhm,
        );

        assert!(snr > 10.0, "SNR {} should be >> 10 for bright star", snr);
    }

    #[test]
    fn test_snr_weight() {
        let data = vec![1000.0_f32; 10000];
        let w = compute_snr_weight(&data, 1000.0, 50.0);
        // All pixels at background → MeanDev ≈ 0 → SNRWeight ≈ 0
        assert!(w < 0.01, "SNRWeight {} should be ~0 for flat image", w);
    }

    #[test]
    fn test_psf_signal() {
        let stars = vec![
            MeasuredStar {
                x: 0.0, y: 0.0, peak: 5000.0, flux: 50000.0,
                fwhm_x: 7.0, fwhm_y: 7.0, fwhm: 7.0, eccentricity: 0.0, hfr: 3.5, snr: 0.0,
            },
            MeasuredStar {
                x: 0.0, y: 0.0, peak: 3000.0, flux: 30000.0,
                fwhm_x: 7.0, fwhm_y: 7.0, fwhm: 7.0, eccentricity: 0.0, hfr: 3.5, snr: 0.0,
            },
            MeasuredStar {
                x: 0.0, y: 0.0, peak: 7000.0, flux: 70000.0,
                fwhm_x: 7.0, fwhm_y: 7.0, fwhm: 7.0, eccentricity: 0.0, hfr: 3.5, snr: 0.0,
            },
        ];

        let psf_sig = compute_psf_signal(&stars, 50.0);
        // Median of [5000, 3000, 7000] = 5000, PSF signal = 5000/50 = 100
        assert!(
            (psf_sig - 100.0).abs() < 1.0,
            "PSF signal {} expected ~100",
            psf_sig
        );
    }
}
