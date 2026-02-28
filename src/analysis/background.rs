/// Background estimation: sigma-clipped statistics and optional mesh-grid spatial background.

use crate::processing::stretch::find_median;

/// Result of background estimation.
pub(crate) struct BackgroundResult {
    /// Global background level (ADU).
    pub background: f32,
    /// Background noise estimate (sigma, ADU).
    pub noise: f32,
    /// Per-pixel background map (only if mesh-grid mode was used).
    pub background_map: Option<Vec<f32>>,
}

/// Estimate background and noise using sigma-clipped statistics.
/// Subsamples to ~500k pixels, runs 3 rounds of 3-sigma clipping.
pub(crate) fn estimate_background(data: &[f32], width: usize, height: usize) -> BackgroundResult {
    let total = width * height;
    let target_samples = 500_000usize;
    let stride = ((total as f64 / target_samples as f64).sqrt() as usize).max(1);

    // Subsample, skipping 2-pixel border
    let border = 2;
    let mut samples = Vec::with_capacity(target_samples);
    let mut y = border;
    while y < height.saturating_sub(border) {
        let mut x = border;
        while x < width.saturating_sub(border) {
            let val = data[y * width + x];
            if val.is_finite() {
                samples.push(val);
            }
            x += stride;
        }
        y += stride;
    }

    if samples.len() < 100 {
        return BackgroundResult {
            background: 0.0,
            noise: 1.0,
            background_map: None,
        };
    }

    // Iterative sigma-clipping: 3 rounds, 3σ threshold
    let (mode, sigma) = sigma_clipped_stats(&mut samples, 3, 3.0);

    BackgroundResult {
        background: mode,
        noise: sigma.max(0.001),
        background_map: None,
    }
}

/// Estimate background with SExtractor-style mesh grid for spatially varying backgrounds.
pub(crate) fn estimate_background_mesh(
    data: &[f32],
    width: usize,
    height: usize,
    cell_size: usize,
) -> BackgroundResult {
    let cell_size = cell_size.max(16);
    let nx = (width + cell_size - 1) / cell_size;
    let ny = (height + cell_size - 1) / cell_size;

    let mut cell_bg = vec![0.0_f32; nx * ny];
    let mut cell_sigma = vec![0.0_f32; nx * ny];
    let mut cell_valid = vec![true; nx * ny];

    // Per-cell statistics
    for cy in 0..ny {
        let y0 = cy * cell_size;
        let y1 = (y0 + cell_size).min(height);
        for cx in 0..nx {
            let x0 = cx * cell_size;
            let x1 = (x0 + cell_size).min(width);
            let cell_idx = cy * nx + cx;

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
                cell_valid[cell_idx] = false;
                continue;
            }

            let original_len = samples.len();
            let (mode, sigma) = sigma_clipped_stats(&mut samples, 3, 3.0);

            // Reject star-contaminated cells (>30% clipped)
            if samples.len() < (original_len * 7 / 10) {
                cell_valid[cell_idx] = false;
            }

            cell_bg[cell_idx] = mode;
            cell_sigma[cell_idx] = sigma;
        }
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

    // Bilinear interpolation to full resolution (parallelized)
    let mut bg_map = vec![0.0_f32; width * height];

    if nx < 2 || ny < 2 {
        // Too few cells for bilinear interpolation — fill with flat background
        let flat_val = filtered_bg[0];
        bg_map.fill(flat_val);
    } else {
        let half_cell = cell_size as f32 * 0.5;
        let inv_cell = 1.0 / cell_size as f32;
        let nx_clamp = nx as i32 - 2;

        use rayon::prelude::*;
        bg_map
            .par_chunks_mut(width)
            .enumerate()
            .for_each(|(y, row)| {
                let fy = (y as f32 - half_cell) * inv_cell;
                let iy = (fy.floor() as i32).clamp(0, ny as i32 - 2) as usize;
                let ty = (fy - iy as f32).clamp(0.0, 1.0);
                let row0 = &filtered_bg[iy * nx..(iy + 1) * nx];
                let row1 = &filtered_bg[(iy + 1) * nx..(iy + 2) * nx];
                let ty_inv = 1.0 - ty;

                for (x, dst) in row.iter_mut().enumerate() {
                    let fx = (x as f32 - half_cell) * inv_cell;
                    let ix = (fx.floor() as i32).clamp(0, nx_clamp) as usize;
                    let tx = (fx - ix as f32).clamp(0.0, 1.0);
                    let tx_inv = 1.0 - tx;

                    *dst = row0[ix] * tx_inv * ty_inv
                        + row0[ix + 1] * tx * ty_inv
                        + row1[ix] * tx_inv * ty
                        + row1[ix + 1] * tx * ty;
                }
            });
    }

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

    BackgroundResult {
        background,
        noise: noise.max(0.001),
        background_map: Some(bg_map),
    }
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

    // Mode approximation (SExtractor formula)
    let mode = 2.5 * median - 1.5 * mean;

    (mode, sigma.max(0.001))
}

#[cfg(test)]
mod tests {
    use super::*;

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

        let result = estimate_background(&data, width, height);
        assert!(
            (result.background - bg_level).abs() < 30.0,
            "bg {} should be ~{}",
            result.background,
            bg_level
        );
        assert!(result.noise > 0.0, "noise should be > 0");
        assert!(result.background_map.is_none());
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
}
