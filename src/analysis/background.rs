/// Background estimation: sigma-clipped statistics and optional mesh-grid spatial background.

use crate::processing::stretch::find_median;

/// Result of background estimation.
pub struct BackgroundResult {
    /// Global background level (ADU).
    pub background: f32,
    /// Background noise estimate (sigma, ADU).
    pub noise: f32,
    /// Per-pixel background map (only if mesh-grid mode was used).
    pub background_map: Option<Vec<f32>>,
    /// Per-pixel noise map (bilinear interpolation of per-cell sigma).
    /// Only present if mesh-grid mode was used.
    pub noise_map: Option<Vec<f32>>,
}

/// Compute auto-tuned mesh cell size from image dimensions.
/// Targets ~32 cells across the longest axis, minimum 16 pixels per cell.
pub fn auto_cell_size(width: usize, height: usize) -> usize {
    (width.max(height) / 32).max(16)
}

/// Estimate background and noise using sigma-clipped statistics.
/// Subsamples to ~500k pixels, runs 3 rounds of 3-sigma clipping.
/// Estimate background with SExtractor-style mesh grid for spatially varying backgrounds.
pub fn estimate_background_mesh(
    data: &[f32],
    width: usize,
    height: usize,
    cell_size: usize,
) -> BackgroundResult {
    let cell_size = cell_size.max(16);
    let nx = (width + cell_size - 1) / cell_size;
    let ny = (height + cell_size - 1) / cell_size;

    // Per-cell statistics (parallelized — each cell is independent)
    use rayon::prelude::*;
    let cells: Vec<(f32, f32, bool)> = (0..ny * nx)
        .into_par_iter()
        .map(|cell_idx| {
            let cy = cell_idx / nx;
            let cx = cell_idx % nx;
            let y0 = cy * cell_size;
            let y1 = (y0 + cell_size).min(height);
            let x0 = cx * cell_size;
            let x1 = (x0 + cell_size).min(width);

            let mut samples = Vec::with_capacity((y1 - y0) * (x1 - x0));
            for y in y0..y1 {
                for x in x0..x1 {
                    let val = data[y * width + x];
                    if val.is_finite() {
                        samples.push(val);
                    }
                }
            }

            if samples.len() < 10 {
                return (0.0, 0.0, false);
            }

            let original_len = samples.len();
            let (mode, sigma) = sigma_clipped_stats(&mut samples, 3, 3.0);

            // Reject star-contaminated cells (>30% clipped)
            let valid = samples.len() >= (original_len * 7 / 10);
            (mode, sigma, valid)
        })
        .collect();

    let mut cell_bg = vec![0.0_f32; nx * ny];
    let mut cell_sigma = vec![0.0_f32; nx * ny];
    let mut cell_valid = vec![true; nx * ny];
    for (i, &(bg, sigma, valid)) in cells.iter().enumerate() {
        cell_bg[i] = bg;
        cell_sigma[i] = sigma;
        cell_valid[i] = valid;
    }

    // Fill invalid cells with nearest valid neighbor
    for cy in 0..ny {
        for cx in 0..nx {
            let idx = cy * nx + cx;
            if cell_valid[idx] {
                continue;
            }
            // Search expanding ring
            let mut best_dist = usize::MAX;
            let mut best_val = 0.0_f32;
            let mut best_sig = 1.0_f32;
            let max_r = nx.max(ny);
            for r in 1..=max_r {
                let mut found = false;
                let y_lo = cy.saturating_sub(r);
                let y_hi = (cy + r + 1).min(ny);
                let x_lo = cx.saturating_sub(r);
                let x_hi = (cx + r + 1).min(nx);
                for sy in y_lo..y_hi {
                    for sx in x_lo..x_hi {
                        let sidx = sy * nx + sx;
                        if cell_valid[sidx] {
                            let d = cx.abs_diff(sx) + cy.abs_diff(sy);
                            if d < best_dist {
                                best_dist = d;
                                best_val = cell_bg[sidx];
                                best_sig = cell_sigma[sidx];
                                found = true;
                            }
                        }
                    }
                }
                if found {
                    break;
                }
            }
            cell_bg[idx] = best_val;
            cell_sigma[idx] = best_sig;
        }
    }

    // 3×3 median filter on grid to suppress remaining star contamination
    let mut filtered_bg = cell_bg.clone();
    let mut neighbors = [0.0_f32; 9];
    for cy in 0..ny {
        for cx in 0..nx {
            let mut count = 0;
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let sy = cy as i32 + dy;
                    let sx = cx as i32 + dx;
                    if sy >= 0 && sy < ny as i32 && sx >= 0 && sx < nx as i32 {
                        neighbors[count] = cell_bg[sy as usize * nx + sx as usize];
                        count += 1;
                    }
                }
            }
            filtered_bg[cy * nx + cx] = find_median(&mut neighbors[..count]);
        }
    }

    // Bicubic spline interpolation to full resolution (parallelized)
    let bg_map = interpolate_grid_to_map(&filtered_bg, nx, ny, width, height, cell_size);

    // Global noise = median of cell sigmas
    let mut valid_sigmas: Vec<f32> = cell_sigma
        .iter()
        .zip(cell_valid.iter())
        .filter(|(_, &v)| v)
        .map(|(&s, _)| s)
        .collect();

    let noise = if valid_sigmas.is_empty() {
        1.0
    } else {
        find_median(&mut valid_sigmas)
    };

    // Global background = median of filtered cell backgrounds
    let mut valid_bgs: Vec<f32> = (0..nx * ny).map(|i| filtered_bg[i]).collect();
    let background = find_median(&mut valid_bgs);

    // Bicubic interpolation of per-cell noise to full resolution
    let noise_map = interpolate_grid_to_map(&cell_sigma, nx, ny, width, height, cell_size);

    BackgroundResult {
        background,
        noise: noise.max(0.001),
        background_map: Some(bg_map),
        noise_map: Some(noise_map),
    }
}

/// Noise estimation via iterative Multiresolution Support (MRS).
///
/// Algorithm (matches PixInsight's approach):
/// 1. Compute à trous wavelet coefficients at layers 1..noise_layers
/// 2. At each layer, identify "significant" pixels (|w_j| > 3σ_j)
/// 3. Build a cumulative significance mask across all layers
/// 4. Re-estimate noise from layer-1 coefficients using only unmasked pixels
/// 5. Iterate until convergence (σ changes < 0.1%)
///
/// `noise_layers`: number of wavelet scales for significance detection (1-6).
/// Layer 1 = finest scale (pixel noise), higher layers = coarser structures.
/// More layers → better structure rejection → purer noise estimate.
/// Default: 4 (matches PixInsight). Minimum: 1.
pub fn estimate_noise_mrs(
    data: &[f32],
    width: usize,
    height: usize,
    noise_layers: usize,
) -> f32 {
    use super::convolution::{b3_spline_smooth, b3_spline_smooth_dilated};

    if data.is_empty() {
        return 0.001;
    }

    let total = width * height;
    let noise_layers = noise_layers.clamp(1, 6);

    // ── Layer 1: compute wavelet coefficients w1 = data - smooth(data) ──
    let mut smoothed = vec![0.0_f32; total];
    b3_spline_smooth(data, width, height, &mut smoothed);

    let mut w1 = vec![0.0_f32; total];
    for i in 0..total {
        w1[i] = data[i] - smoothed[i];
    }

    // If only 1 layer requested, use simple MAD on w1 (original behavior)
    if noise_layers == 1 {
        return simple_mad_noise(&w1, width, height);
    }

    // ── Build significance mask across all wavelet layers ──
    // true = significant (contains structure), excluded from noise estimate
    let mut sig_mask = vec![false; total];

    // Layer 1: mark significant pixels in w1
    let sigma1 = simple_mad_noise(&w1, width, height);
    let thresh1 = 3.0 * sigma1;
    for i in 0..total {
        if w1[i].abs() > thresh1 {
            sig_mask[i] = true;
        }
    }

    // Layers 2..N: compute dilated wavelet coefficients, mark significant pixels
    // Each layer operates on the smoothed output of the previous layer
    let mut prev_smooth = smoothed;

    for layer in 2..=noise_layers {
        let mut next_smooth = vec![0.0_f32; total];
        b3_spline_smooth_dilated(&prev_smooth, width, height, &mut next_smooth, layer);

        // Wavelet coefficients at this layer: w_j = c_{j-1} - c_j
        // Subsample for sigma estimation
        let target_samples = 500_000usize;
        let stride = ((total as f64 / target_samples as f64).sqrt() as usize).max(1);
        let border = (2_usize << (layer - 1)).min(width / 4);

        let mut wj_samples = Vec::with_capacity(target_samples);
        let mut y = border;
        while y < height.saturating_sub(border) {
            let mut x = border;
            while x < width.saturating_sub(border) {
                let idx = y * width + x;
                if !sig_mask[idx] {
                    let coeff = prev_smooth[idx] - next_smooth[idx];
                    if coeff.is_finite() {
                        wj_samples.push(coeff);
                    }
                }
                x += stride;
            }
            y += stride;
        }

        if wj_samples.len() >= 100 {
            let med = find_median(&mut wj_samples);
            let mut abs_devs: Vec<f32> = wj_samples.iter().map(|&v| (v - med).abs()).collect();
            let mad = find_median(&mut abs_devs);
            let sigma_j = 1.4826 * mad;
            let thresh_j = 3.0 * sigma_j;

            // Mark significant pixels at this scale
            for idx in 0..total {
                if !sig_mask[idx] {
                    let coeff = prev_smooth[idx] - next_smooth[idx];
                    if coeff.abs() > thresh_j {
                        sig_mask[idx] = true;
                    }
                }
            }
        }

        prev_smooth = next_smooth;
    }

    // ── Re-estimate noise from layer-1 coefficients, excluding masked pixels ──
    // Iterate until convergence
    let target_samples = 500_000usize;
    let stride = ((total as f64 / target_samples as f64).sqrt() as usize).max(1);
    let border = 2usize;

    let mut sigma = sigma1;
    for _iteration in 0..4 {
        let thresh = 3.0 * sigma;

        let mut samples = Vec::with_capacity(target_samples);
        let mut y = border;
        while y < height.saturating_sub(border) {
            let mut x = border;
            while x < width.saturating_sub(border) {
                let idx = y * width + x;
                if !sig_mask[idx] && w1[idx].abs() <= thresh {
                    samples.push(w1[idx]);
                }
                x += stride;
            }
            y += stride;
        }

        if samples.len() < 100 {
            break;
        }

        let med = find_median(&mut samples);
        let mut abs_devs: Vec<f32> = samples.iter().map(|&v| (v - med).abs()).collect();
        let mad = find_median(&mut abs_devs);
        let new_sigma = 1.4826 * mad;

        if new_sigma <= 0.0 || (new_sigma - sigma).abs() / sigma < 0.001 {
            sigma = new_sigma;
            break;
        }
        sigma = new_sigma;
    }

    sigma.max(0.001)
}

/// Simple MAD-based noise estimate from wavelet coefficients (subsampled).
fn simple_mad_noise(w1: &[f32], width: usize, height: usize) -> f32 {
    let total = width * height;
    let target_samples = 500_000usize;
    let stride = ((total as f64 / target_samples as f64).sqrt() as usize).max(1);
    let border = 2usize;

    let mut samples = Vec::with_capacity(target_samples);
    let mut y = border;
    while y < height.saturating_sub(border) {
        let mut x = border;
        while x < width.saturating_sub(border) {
            let idx = y * width + x;
            let coeff = w1[idx];
            if coeff.is_finite() {
                samples.push(coeff);
            }
            x += stride;
        }
        y += stride;
    }

    if samples.len() < 100 {
        return 0.001;
    }

    let median = find_median(&mut samples);
    let mut abs_devs: Vec<f32> = samples.iter().map(|&v| (v - median).abs()).collect();
    let mad = find_median(&mut abs_devs);
    (1.4826 * mad).max(0.001)
}

/// Iterative sigma-clipped statistics. Returns (mode, sigma_MAD).
/// Mode = 2.5*median - 1.5*mean (SExtractor formula).
/// Modifies the slice in-place (quickselect requirement).
fn sigma_clipped_stats(samples: &mut Vec<f32>, rounds: usize, kappa: f32) -> (f32, f32) {
    let mut median = 0.0_f32;
    let mut sigma = 0.0_f32;
    let mut abs_devs = Vec::with_capacity(samples.len());

    for _ in 0..rounds {
        if samples.len() < 10 {
            break;
        }

        median = find_median(samples);

        // MAD = median(|xi - median|)
        abs_devs.clear();
        abs_devs.extend(samples.iter().map(|&x| (x - median).abs()));
        let mad = find_median(&mut abs_devs);
        sigma = 1.4826 * mad;

        if sigma < 1e-6 {
            break;
        }

        // Clip
        let lo = median - kappa * sigma;
        let hi = median + kappa * sigma;
        samples.retain(|&x| x >= lo && x <= hi);
    }

    if samples.len() < 10 {
        return (median, sigma.max(0.001));
    }

    // Final statistics on clipped data
    median = find_median(samples);
    abs_devs.clear();
    abs_devs.extend(samples.iter().map(|&x| (x - median).abs()));
    let mad = find_median(&mut abs_devs);
    sigma = 1.4826 * mad;

    let mean: f32 = samples.iter().sum::<f32>() / samples.len() as f32;

    // Mode approximation (SExtractor formula) with asymmetry fallback.
    // When |mean - median| >= 0.3 * sigma, the distribution is too skewed
    // (bright nebulosity, galaxy core) for the mode formula to be reliable.
    let mode = if (mean - median).abs() < 0.3 * sigma {
        2.5 * median - 1.5 * mean
    } else {
        median
    };

    (mode, sigma.max(0.001))
}

/// Catmull-Rom cubic weights for parameter t ∈ [0, 1].
/// Returns weights for grid points at offsets [-1, 0, +1, +2] from the base cell.
#[inline]
fn cubic_weights(t: f32) -> [f32; 4] {
    let t2 = t * t;
    let t3 = t2 * t;
    [
        -0.5 * t3 + t2 - 0.5 * t,
        1.5 * t3 - 2.5 * t2 + 1.0,
        -1.5 * t3 + 2.0 * t2 + 0.5 * t,
        0.5 * t3 - 0.5 * t2,
    ]
}

/// Bicubic (Catmull-Rom) spline interpolation of a cell grid to full-resolution pixel map.
/// C¹-continuous, passes through grid values, no overshoot for monotone data.
fn interpolate_grid_to_map(
    grid: &[f32],
    nx: usize,
    ny: usize,
    width: usize,
    height: usize,
    cell_size: usize,
) -> Vec<f32> {
    let mut map = vec![0.0_f32; width * height];

    if nx < 2 || ny < 2 {
        map.fill(grid[0]);
        return map;
    }

    let half_cell = cell_size as f32 * 0.5;
    let inv_cell = 1.0 / cell_size as f32;
    let nx_i = nx as i32;
    let ny_i = ny as i32;

    use rayon::prelude::*;
    map.par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            let fy = (y as f32 - half_cell) * inv_cell;
            let iy = fy.floor() as i32;
            let ty = (fy - iy as f32).clamp(0.0, 1.0);
            let wy = cubic_weights(ty);

            for (x, dst) in row.iter_mut().enumerate() {
                let fx = (x as f32 - half_cell) * inv_cell;
                let ix = fx.floor() as i32;
                let tx = (fx - ix as f32).clamp(0.0, 1.0);
                let wx = cubic_weights(tx);

                let mut val = 0.0_f32;
                for j in 0..4_i32 {
                    let jy = (iy + j - 1).clamp(0, ny_i - 1) as usize;
                    for i in 0..4_i32 {
                        let jx = (ix + i - 1).clamp(0, nx_i - 1) as usize;
                        val += wy[j as usize] * wx[i as usize] * grid[jy * nx + jx];
                    }
                }
                *dst = val;
            }
        });

    map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_auto_cell_size() {
        assert_eq!(auto_cell_size(1024, 768), 32);
        assert_eq!(auto_cell_size(4096, 2160), 128);
        assert_eq!(auto_cell_size(256, 256), 16);
    }

    #[test]
    fn test_sigma_clipped_gaussian() {
        // Generate N(1000, 50) with 5% outliers at 5000
        let mut rng = 42u64;
        let mut samples = Vec::with_capacity(10000);
        for _ in 0..10000 {
            // Box-Muller
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let u1 = (rng >> 11) as f64 / (1u64 << 53) as f64;
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let u2 = (rng >> 11) as f64 / (1u64 << 53) as f64;
            let u1 = u1.max(1e-15);
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            let val = 1000.0 + 50.0 * z;

            // 5% outliers
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let outlier_chance = (rng >> 11) as f64 / (1u64 << 53) as f64;
            if outlier_chance < 0.05 {
                samples.push(5000.0_f32);
            } else {
                samples.push(val as f32);
            }
        }

        let (mode, sigma) = sigma_clipped_stats(&mut samples, 3, 3.0);
        assert!(
            (mode - 1000.0).abs() < 15.0,
            "mode {} should be ~1000",
            mode
        );
        assert!(
            (sigma - 50.0).abs() < 10.0,
            "sigma {} should be ~50",
            sigma
        );
    }

    #[test]
    fn test_background_estimation() {
        // Flat image with noise
        let width = 200;
        let height = 200;
        let bg_level = 1000.0_f32;
        let mut data = vec![bg_level; width * height];

        // Add some pseudo-noise
        let mut rng = 123u64;
        for val in data.iter_mut() {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let noise = ((rng >> 33) as f32 / (1u64 << 31) as f32 - 0.5) * 100.0;
            *val += noise;
        }

        let cell_size = auto_cell_size(width, height);
        let result = estimate_background_mesh(&data, width, height, cell_size);
        assert!(
            (result.background - bg_level).abs() < 30.0,
            "bg {} should be ~{}",
            result.background,
            bg_level
        );
        assert!(result.noise > 0.0, "noise should be > 0");
    }

    #[test]
    fn test_background_mesh() {
        // Linear gradient: 500 at left, 1500 at right
        let width = 256;
        let height = 256;
        let mut data = vec![0.0_f32; width * height];
        for y in 0..height {
            for x in 0..width {
                data[y * width + x] = 500.0 + 1000.0 * (x as f32 / (width - 1) as f32);
            }
        }

        let result = estimate_background_mesh(&data, width, height, 64);
        assert!(result.background_map.is_some());
        let bg_map = result.background_map.unwrap();

        // Check that the background map captures the gradient direction
        // (left side should be lower than right side)
        let left_bg = bg_map[height / 2 * width + width / 4];
        let right_bg = bg_map[height / 2 * width + 3 * width / 4];
        assert!(
            left_bg < right_bg,
            "left bg {} should be < right bg {}",
            left_bg,
            right_bg
        );
        assert!(
            left_bg < 900.0,
            "left bg {} should be < 900",
            left_bg
        );
        assert!(
            right_bg > 1100.0,
            "right bg {} should be > 1100",
            right_bg
        );
    }

    #[test]
    fn test_mrs_significance_masking() {
        // Image with bright extended structure on flat background.
        // Multi-layer MRS should give a LOWER noise estimate than single-layer
        // because it masks out the structure from the noise sample.
        let w = 200;
        let h = 200;
        let bg_level = 1000.0_f32;
        let mut data = vec![bg_level; w * h];

        // Add deterministic pseudo-noise
        let true_sigma = 30.0_f32;
        for i in 0..data.len() {
            let r = ((i as u64).wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407) >> 33) as f32;
            data[i] += true_sigma * (r / (1u64 << 31) as f32 - 0.5) * 2.0;
        }

        // Add a bright extended structure (simulated nebula patch)
        for y in 60..140 {
            for x in 60..140 {
                let dx = (x as f32 - 100.0) / 30.0;
                let dy = (y as f32 - 100.0) / 30.0;
                data[y * w + x] += 2000.0 * (-0.5 * (dx * dx + dy * dy)).exp();
            }
        }

        let noise_1layer = estimate_noise_mrs(&data, w, h, 1);
        let noise_4layer = estimate_noise_mrs(&data, w, h, 4);

        eprintln!("1-layer noise: {:.2}, 4-layer noise: {:.2}, true: {:.2}",
            noise_1layer, noise_4layer, true_sigma);

        // With bright structure, 1-layer estimate may be inflated.
        // 4-layer should be <= 1-layer (better structure rejection).
        assert!(
            noise_4layer <= noise_1layer,
            "4-layer ({:.2}) should be <= 1-layer ({:.2}) with bright structure",
            noise_4layer, noise_1layer,
        );
        // 4-layer should be reasonably close to truth
        assert!(
            (noise_4layer - true_sigma).abs() < true_sigma * 0.5,
            "4-layer noise {:.2} should be within 50% of true sigma {:.2}",
            noise_4layer, true_sigma,
        );
    }
}
