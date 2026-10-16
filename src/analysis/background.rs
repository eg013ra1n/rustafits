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
///
/// `want_noise_map`: the per-pixel noise map is a second bicubic
/// interpolation pass over the cell grid (the same cost as the background
/// map's own interpolation) — a caller that only reads `background`,
/// `noise` (the scalar) and `background_map` should pass `false` to skip
/// it. The interpolation is a pure leaf (`interpolate_grid_to_map`, no
/// shared mutable state) that reads `cell_sigma` and writes nothing back,
/// so skipping it changes no other field: `background_map`, `background`
/// and `noise` are computed earlier from `filtered_bg`/`cell_sigma`
/// directly, never from `noise_map` itself.
pub fn estimate_background_mesh(
    data: &[f32],
    width: usize,
    height: usize,
    cell_size: usize,
    want_noise_map: bool,
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

    // Bicubic interpolation of per-cell noise to full resolution — skipped
    // when the caller never reads it (see `want_noise_map` doc comment
    // above); this is a pure leaf call with no effect on `background`,
    // `noise` or `background_map`, all already finalized above.
    let noise_map = if want_noise_map {
        Some(interpolate_grid_to_map(
            &cell_sigma,
            nx,
            ny,
            width,
            height,
            cell_size,
        ))
    } else {
        None
    };

    BackgroundResult {
        background,
        noise: noise.max(0.001),
        background_map: Some(bg_map),
        noise_map,
    }
}

/// Noise estimation via iterative Multiresolution Support (MRS).
///
/// Algorithm (matches the reference implementation's approach):
/// 1. Compute à trous wavelet coefficients at layers 1..noise_layers
/// 2. At each layer, identify "significant" pixels (|w_j| > 3σ_j)
/// 3. Build a cumulative significance mask across all layers
/// 4. Re-estimate noise from layer-1 coefficients using only unmasked pixels
/// 5. Iterate until convergence (σ changes < 0.1%)
///
/// `noise_layers`: number of wavelet scales for significance detection (1-6).
/// Layer 1 = finest scale (pixel noise), higher layers = coarser structures.
/// More layers → better structure rejection → purer noise estimate.
/// Default: 4 (the usual default). Minimum: 1.
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
    for (m, &w) in sig_mask.iter_mut().zip(w1.iter()) {
        *m = w.abs() > thresh1;
    }

    // Layers 2..N: compute dilated wavelet coefficients, mark significant pixels
    // Each layer operates on the smoothed output of the previous layer.
    //
    // Tier A task-3b, item 1: `next_smooth` used to be a fresh `vec![0.0;
    // total]` allocated inside the loop on every layer (up to 5 more
    // full-plane allocations at `noise_layers = 6`), each one immediately
    // superseding the previous `prev_smooth` once `prev_smooth =
    // next_smooth` ran. Two named buffers now ping-pong via `mem::swap`
    // instead: `prev_smooth` starts as the layer-1 `smoothed` plane (moved
    // in, no copy), `next_smooth` is the ONE extra buffer this loop ever
    // allocates, and every iteration writes into whichever buffer currently
    // holds the STALE (already-consumed) values before swapping the names.
    // `b3_spline_smooth_dilated` writes every element of its `output`
    // parameter (`par_chunks_mut(width)` covers every row, every `row[x] =
    // sum`) — the same total-overwrite guarantee `vec![0.0; total]` gave,
    // so reusing a buffer that still holds a prior layer's numbers is safe:
    // nothing is ever read from it before this call fills it.
    let mut prev_smooth = smoothed;
    let mut next_smooth = vec![0.0_f32; total];

    for layer in 2..=noise_layers {
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
            // `mask |= |c| > t` is the old "unmasked and over the threshold -> mask" without the
            // branch (a masked pixel stays masked, NaN compares false), so it vectorizes.
            for ((m, &p), &n) in sig_mask.iter_mut().zip(prev_smooth.iter()).zip(next_smooth.iter()) {
                *m |= (p - n).abs() > thresh_j;
            }
        }

        std::mem::swap(&mut prev_smooth, &mut next_smooth);
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

/// The median of `|x - median|` over `samples`: the same element `find_median` returns on the
/// deviation vector (the `len/2`-th smallest; a deviation is a non-negative f32, whose bit pattern
/// orders like its value, and equal values are the same bits), found by an integer selection on
/// reusable scratch. `samples` is only read, so its permutation (which the caller's final mean
/// sum depends on) is untouched.
fn median_abs_dev(samples: &[f32], median: f32, scratch: &mut Vec<u32>) -> f32 {
    if samples.is_empty() {
        return 0.0;
    }
    scratch.clear();
    scratch.extend(samples.iter().map(|&x| (x - median).abs().to_bits()));
    let k = scratch.len() / 2;
    let (_, v, _) = scratch.select_nth_unstable(k);
    f32::from_bits(*v)
}

/// Iterative sigma-clipped statistics. Returns (mode, sigma_MAD).
/// Mode = 2.5*median - 1.5*mean (SExtractor formula).
/// Modifies the slice in-place (quickselect requirement).
fn sigma_clipped_stats(samples: &mut Vec<f32>, rounds: usize, kappa: f32) -> (f32, f32) {
    let mut median = 0.0_f32;
    let mut sigma = 0.0_f32;
    let mut abs_devs: Vec<u32> = Vec::with_capacity(samples.len());

    for _ in 0..rounds {
        if samples.len() < 10 {
            break;
        }

        median = find_median(samples);

        // MAD = median(|xi - median|)
        let mad = median_abs_dev(samples, median, &mut abs_devs);
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
    let mad = median_abs_dev(samples, median, &mut abs_devs);
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

    // Per-column tables: the x floor / clamp / weights / clamped grid columns depend on x only,
    // so they are computed once per column instead of once per pixel (the same expressions on
    // the same values, hence the same bits).
    let cols: Vec<([usize; 4], [f32; 4])> = (0..width)
        .map(|x| {
            let fx = (x as f32 - half_cell) * inv_cell;
            let ix = fx.floor() as i32;
            let tx = (fx - ix as f32).clamp(0.0, 1.0);
            let mut jx = [0usize; 4];
            for (i, slot) in jx.iter_mut().enumerate() {
                *slot = (ix + i as i32 - 1).clamp(0, nx_i - 1) as usize;
            }
            (jx, cubic_weights(tx))
        })
        .collect();

    use rayon::prelude::*;
    map.par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            let fy = (y as f32 - half_cell) * inv_cell;
            let iy = fy.floor() as i32;
            let ty = (fy - iy as f32).clamp(0.0, 1.0);
            let wy = cubic_weights(ty);
            let mut grows: [&[f32]; 4] = [&[]; 4];
            for (j, slot) in grows.iter_mut().enumerate() {
                let jy = (iy + j as i32 - 1).clamp(0, ny_i - 1) as usize;
                *slot = &grid[jy * nx..jy * nx + nx];
            }

            for (dst, (jx, wx)) in row.iter_mut().zip(cols.iter()) {
                let mut val = 0.0_f32;
                for j in 0..4 {
                    let g = grows[j];
                    for i in 0..4 {
                        val += wy[j] * wx[i] * g[jx[i]];
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
        let result = estimate_background_mesh(&data, width, height, cell_size, true);
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

        let result = estimate_background_mesh(&data, width, height, 64, true);
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

    /// D4: `want_noise_map = false` must not move `background`, `noise` or
    /// `background_map` by a single bit — the noise-map interpolation is a
    /// pure leaf with no side effect on the rest of `BackgroundResult`, and
    /// this pins that claim rather than just asserting it in a comment.
    #[test]
    fn background_map_identical_with_noise_map_flag_either_way() {
        // Star-field-like image: flat background plus a few bright blobs,
        // so the >30%-clipped / nearest-neighbor-fill path is exercised too.
        let width = 200;
        let height = 200;
        let bg_level = 800.0_f32;
        let mut data = vec![bg_level; width * height];
        let mut rng = 42u64;
        for val in data.iter_mut() {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let noise = ((rng >> 33) as f32 / (1u64 << 31) as f32 - 0.5) * 40.0;
            *val += noise;
        }
        // A few saturated "stars" to force some cells past the contamination
        // threshold and through the neighbor-fill branch.
        for &(cx, cy) in &[(30usize, 30usize), (150, 60), (90, 170)] {
            for dy in 0..6 {
                for dx in 0..6 {
                    data[(cy + dy) * width + (cx + dx)] = 50_000.0;
                }
            }
        }

        let cell_size = auto_cell_size(width, height);
        let with_noise = estimate_background_mesh(&data, width, height, cell_size, true);
        let without_noise = estimate_background_mesh(&data, width, height, cell_size, false);

        assert_eq!(with_noise.background, without_noise.background);
        assert_eq!(with_noise.noise, without_noise.noise);
        assert_eq!(with_noise.background_map, without_noise.background_map);
        assert!(with_noise.noise_map.is_some());
        assert!(without_noise.noise_map.is_none());
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

    /// Deterministic synthetic plane for the Tier A task-3b bit-identity
    /// pins below: a gradient background, xorshift pseudo-noise, several
    /// Gaussian "stars", and one saturated blob — the mix `estimate_noise_mrs`
    /// exists to be robust against (structure the significance mask must
    /// suppress, plus a non-finite-free saturated region to exercise the
    /// border/clamped index paths).
    fn synthetic_noise_fixture(w: usize, h: usize) -> Vec<f32> {
        let mut data = vec![0.0_f32; w * h];
        for y in 0..h {
            for x in 0..w {
                data[y * w + x] = 500.0 + 800.0 * (x as f32 / (w - 1) as f32);
            }
        }
        // Deterministic xorshift64 pseudo-noise.
        let mut rng: u64 = 88172645463325252;
        for v in data.iter_mut() {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            let u = (rng >> 40) as f32 / (1u64 << 24) as f32;
            *v += (u - 0.5) * 40.0;
        }
        // A handful of Gaussian "stars" at fixed positions/amplitudes.
        let stars: [(usize, usize, f32, f32); 6] = [
            (60, 60, 2.5, 3000.0),
            (200, 150, 2.0, 5000.0),
            (350, 100, 1.5, 2000.0),
            (100, 300, 3.0, 8000.0),
            (400, 400, 1.0, 1500.0),
            (250, 420, 1.8, 4000.0),
        ];
        for &(cx, cy, r, amp) in &stars {
            let radius = (r * 6.0).ceil() as i64;
            for dy in -radius..=radius {
                for dx in -radius..=radius {
                    let x = cx as i64 + dx;
                    let y = cy as i64 + dy;
                    if x < 0 || y < 0 || x as usize >= w || y as usize >= h {
                        continue;
                    }
                    let fx = dx as f32 / r;
                    let fy = dy as f32 / r;
                    data[y as usize * w + x as usize] += amp * (-0.5 * (fx * fx + fy * fy)).exp();
                }
            }
        }
        // One saturated blob.
        for y in 240..260.min(h) {
            for x in 240..260.min(w) {
                data[y * w + x] = 65000.0;
            }
        }
        data
    }

    /// Task 3b pin (controller ruling R-TA-2): `estimate_noise_mrs` returns
    /// the SAME `f32`, bit-exact, for every `noise_layers` this crate's
    /// caller exercises (1 — the `simple_mad_noise` early return; 2 — the
    /// shortest layer loop; 4 — production's `MRS_LAYERS`; 6 — the clamp
    /// ceiling), before and after the ping-pong/bit-packed-mask rewrite.
    /// Values recorded from the CURRENT (pre-task-3b) code.
    #[test]
    fn estimate_noise_mrs_bit_exact_pins() {
        let (w, h) = (512, 512);
        let data = synthetic_noise_fixture(w, h);

        let n1 = estimate_noise_mrs(&data, w, h, 1);
        let n2 = estimate_noise_mrs(&data, w, h, 2);
        let n4 = estimate_noise_mrs(&data, w, h, 4);
        let n6 = estimate_noise_mrs(&data, w, h, 6);

        assert_eq!(n1, 12.759915_f32, "noise_layers=1");
        assert_eq!(n2, 12.710779_f32, "noise_layers=2");
        assert_eq!(n4, 12.720099_f32, "noise_layers=4");
        assert_eq!(n6, 12.748694_f32, "noise_layers=6");
    }

    // ---- bit-for-bit oracles: the pre-optimization bodies, kept verbatim ----

    fn sigma_clipped_stats_ref(samples: &mut Vec<f32>, rounds: usize, kappa: f32) -> (f32, f32) {
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


    fn interpolate_ref(
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


    fn mesh_ref(
        data: &[f32],
        width: usize,
        height: usize,
        cell_size: usize,
        want_noise_map: bool,
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
                let (mode, sigma) = sigma_clipped_stats_ref(&mut samples, 3, 3.0);

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
        let bg_map = interpolate_ref(&filtered_bg, nx, ny, width, height, cell_size);

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

        // Bicubic interpolation of per-cell noise to full resolution — skipped
        // when the caller never reads it (see `want_noise_map` doc comment
        // above); this is a pure leaf call with no effect on `background`,
        // `noise` or `background_map`, all already finalized above.
        let noise_map = if want_noise_map {
            Some(interpolate_ref(
                &cell_sigma,
                nx,
                ny,
                width,
                height,
                cell_size,
            ))
        } else {
            None
        };

        BackgroundResult {
            background,
            noise: noise.max(0.001),
            background_map: Some(bg_map),
            noise_map,
        }
    }


    fn simple_mad_noise_ref(w1: &[f32], width: usize, height: usize) -> f32 {
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


    fn mrs_ref(
        data: &[f32],
        width: usize,
        height: usize,
        noise_layers: usize,
    ) -> f32 {
        use super::super::convolution::{b3_spline_smooth_ref as b3_spline_smooth, b3_spline_smooth_dilated_ref as b3_spline_smooth_dilated};

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
            return simple_mad_noise_ref(&w1, width, height);
        }

        // ── Build significance mask across all wavelet layers ──
        // true = significant (contains structure), excluded from noise estimate
        let mut sig_mask = vec![false; total];

        // Layer 1: mark significant pixels in w1
        let sigma1 = simple_mad_noise_ref(&w1, width, height);
        let thresh1 = 3.0 * sigma1;
        for i in 0..total {
            if w1[i].abs() > thresh1 {
                sig_mask[i] = true;
            }
        }

        // Layers 2..N: compute dilated wavelet coefficients, mark significant pixels
        // Each layer operates on the smoothed output of the previous layer.
        //
        // Tier A task-3b, item 1: `next_smooth` used to be a fresh `vec![0.0;
        // total]` allocated inside the loop on every layer (up to 5 more
        // full-plane allocations at `noise_layers = 6`), each one immediately
        // superseding the previous `prev_smooth` once `prev_smooth =
        // next_smooth` ran. Two named buffers now ping-pong via `mem::swap`
        // instead: `prev_smooth` starts as the layer-1 `smoothed` plane (moved
        // in, no copy), `next_smooth` is the ONE extra buffer this loop ever
        // allocates, and every iteration writes into whichever buffer currently
        // holds the STALE (already-consumed) values before swapping the names.
        // `b3_spline_smooth_dilated` writes every element of its `output`
        // parameter (`par_chunks_mut(width)` covers every row, every `row[x] =
        // sum`) — the same total-overwrite guarantee `vec![0.0; total]` gave,
        // so reusing a buffer that still holds a prior layer's numbers is safe:
        // nothing is ever read from it before this call fills it.
        let mut prev_smooth = smoothed;
        let mut next_smooth = vec![0.0_f32; total];

        for layer in 2..=noise_layers {
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

            std::mem::swap(&mut prev_smooth, &mut next_smooth);
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


    fn lcg(seed: &mut u64) -> f32 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 40) as f32) / 16_777_216.0
    }

    /// Star-and-gradient plane with NaN holes, ADU-scale values.
    fn scene(w: usize, h: usize, seed: u64) -> Vec<f32> {
        let mut s = seed;
        let mut d = vec![0.0_f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let n = (lcg(&mut s) + lcg(&mut s) + lcg(&mut s) - 1.5) * 40.0;
                d[y * w + x] = 500.0 + 0.02 * x as f32 + 0.01 * y as f32 + n;
            }
        }
        for _ in 0..(w * h / 900).max(3) {
            let cx = (lcg(&mut s) * w as f32) as isize;
            let cy = (lcg(&mut s) * h as f32) as isize;
            let amp = 200.0 + lcg(&mut s) * 30000.0;
            for dy in -6..=6isize {
                for dx in -6..=6isize {
                    let (x, y) = (cx + dx, cy + dy);
                    if x >= 0 && y >= 0 && (x as usize) < w && (y as usize) < h {
                        d[y as usize * w + x as usize] += amp * (-((dx * dx + dy * dy) as f32) / 6.0).exp();
                    }
                }
            }
        }
        for _ in 0..(w * h / 500) {
            let i = (lcg(&mut s) * (w * h) as f32) as usize % (w * h);
            d[i] = f32::NAN;
        }
        d
    }

    fn assert_bits_eq(a: &[f32], b: &[f32], what: &str) {
        assert_eq!(a.len(), b.len(), "{what} length");
        for i in 0..a.len() {
            assert_eq!(a[i].to_bits(), b[i].to_bits(), "{what} at {i}");
        }
    }

    fn check_plane(d: &[f32], w: usize, h: usize, cell: usize, tag: &str) {
        for &want_noise in &[false, true] {
            let a = estimate_background_mesh(d, w, h, cell, want_noise);
            let b = mesh_ref(d, w, h, cell, want_noise);
            assert_eq!(a.background.to_bits(), b.background.to_bits(), "{tag} background");
            assert_eq!(a.noise.to_bits(), b.noise.to_bits(), "{tag} noise");
            assert_bits_eq(a.background_map.as_ref().unwrap(), b.background_map.as_ref().unwrap(), &format!("{tag} bg map"));
            match (&a.noise_map, &b.noise_map) {
                (Some(x), Some(y)) => assert_bits_eq(x, y, &format!("{tag} noise map")),
                (None, None) => {}
                _ => panic!("{tag} noise map presence"),
            }
        }
        for layers in [1usize, 2, 4, 6] {
            let a = estimate_noise_mrs(d, w, h, layers);
            let b = mrs_ref(d, w, h, layers);
            assert_eq!(a.to_bits(), b.to_bits(), "{tag} mrs layers {layers}");
        }
    }

    #[test]
    fn mesh_and_mrs_are_bit_equal_to_the_pre_optimization_bodies() {
        for &(w, h, cell, seed) in &[
            (96usize, 80usize, 16usize, 1u64),
            (301, 211, 32, 2),
            (640, 480, 64, 3),
            (777, 333, 128, 4),
            (1024, 700, 128, 5),
            (50, 40, 16, 6),
        ] {
            let d = scene(w, h, seed);
            check_plane(&d, w, h, cell, &format!("{w}x{h} c{cell}"));
        }
        let d = scene(20, 20, 9);
        check_plane(&d, 20, 20, 16, "tiny");
        check_plane(&vec![f32::NAN; 64 * 64], 64, 64, 16, "all-nan");
    }

    #[test]
    fn sigma_clipped_stats_is_bit_equal_on_random_vectors() {
        let mut s = 77u64;
        for n in [0usize, 5, 10, 11, 64, 1000, 16384, 20000] {
            for shape in 0..4 {
                let mut v: Vec<f32> = (0..n)
                    .map(|i| match shape {
                        0 => 500.0 + (lcg(&mut s) - 0.5) * 100.0,
                        1 => if i % 17 == 0 { 40000.0 } else { 300.0 + lcg(&mut s) * 20.0 },
                        2 => 1.0,
                        _ => (lcg(&mut s) * 8.0).floor(),
                    })
                    .collect();
                let mut w = v.clone();
                let a = sigma_clipped_stats(&mut v, 3, 3.0);
                let b = sigma_clipped_stats_ref(&mut w, 3, 3.0);
                assert_eq!(a.0.to_bits(), b.0.to_bits(), "mode n={n} shape={shape}");
                assert_eq!(a.1.to_bits(), b.1.to_bits(), "sigma n={n} shape={shape}");
            }
        }
    }

    /// ATH_BG_REAL=<float32 FITS> pins both functions on a real calibrated plane (x65535).
    #[test]
    fn real_plane_is_bit_equal_when_given() {
        let Ok(path) = std::env::var("ATH_BG_REAL") else { return };
        let b = std::fs::read(&path).unwrap();
        let (mut w, mut h, mut off) = (0usize, 0usize, 0usize);
        'o: loop {
            for c in 0..36 {
                let card = std::str::from_utf8(&b[off + c * 80..off + c * 80 + 80]).unwrap();
                let num = |k: &str| card.starts_with(k).then(|| card[10..].split('/').next().unwrap().trim().parse::<usize>().unwrap());
                if let Some(x) = num("NAXIS1  ") { w = x; }
                if let Some(x) = num("NAXIS2  ") { h = x; }
                if card.starts_with("END ") { off += 2880; break 'o; }
            }
            off += 2880;
        }
        let d: Vec<f32> = (0..w * h)
            .map(|i| f32::from_be_bytes(b[off + i * 4..off + i * 4 + 4].try_into().unwrap()) * 65535.0)
            .collect();
        check_plane(&d, w, h, 128, "real");
    }
}
