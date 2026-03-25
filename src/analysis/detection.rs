/// Star detection: DAOFIND-inspired matched filter + proximity blend rejection.

/// A detected star candidate before metric computation.
#[derive(Clone)]
pub struct DetectedStar {
    /// Intensity-weighted centroid X (subpixel).
    pub x: f32,
    /// Intensity-weighted centroid Y (subpixel).
    pub y: f32,
    /// Background-subtracted peak value.
    pub peak: f32,
    /// Total background-subtracted flux.
    pub flux: f32,
    /// Number of connected pixels above threshold.
    pub area: usize,
    /// Position angle from stamp-based I-weighted second-order moments (radians).
    pub theta: f32,
    /// Eccentricity from stamp-based moments (0 = round, approaching 1 = elongated).
    pub eccentricity: f32,
}

/// Detection parameters.
pub struct DetectionParams {
    pub detection_sigma: f32,
    pub min_star_area: usize,
    pub max_star_area: usize,
    pub saturation_limit: f32,
}

impl Default for DetectionParams {
    fn default() -> Self {
        DetectionParams {
            detection_sigma: 5.0,
            min_star_area: 5,
            max_star_area: 2000,
            saturation_limit: 0.95 * 65535.0,
        }
    }
}

/// Detect stars in an image using DAOFIND-style matched filter + proximity blend rejection.
///
/// `data`: single-channel f32 image (raw ADU values, NOT background-subtracted).
/// `background`: global background level.
/// `noise`: background noise sigma.
/// `bg_map`: optional per-pixel background map (from mesh-grid estimation).
/// `noise_map`: optional per-pixel noise map for adaptive thresholds.
/// `fwhm`: estimated FWHM for matched filter kernel (pixels).
/// `field_fwhm`: calibrated field FWHM from Pass 1 (used by Pass 2 filters; None in Pass 1).
pub fn detect_stars(
    data: &[f32],
    width: usize,
    height: usize,
    background: f32,
    noise: f32,
    bg_map: Option<&[f32]>,
    noise_map: Option<&[f32]>,
    params: &DetectionParams,
    fwhm: f32,
    field_fwhm: Option<f32>,
) -> Vec<DetectedStar> {
    // Stage 1: Separable Gaussian convolution + peak detection
    let sigma = fwhm / 2.3548;
    let (kernel_1d, energy_1d) = super::convolution::generate_1d_kernel(sigma);
    let radius = kernel_1d.len() / 2;

    // Separable kernel energy: the 2D kernel K(x,y) = h(x)×v(y), so
    // Σ K² = (Σ h²) × (Σ v²) = energy_1d² (same kernel both axes).
    let kernel_energy_sq = energy_1d * energy_1d;

    // Convolution threshold: detection_sigma × noise × sqrt(Σ K²)
    let threshold = params.detection_sigma * noise * kernel_energy_sq.sqrt();

    // Separable convolution (SIMD-accelerated horizontal + vertical passes)
    let mut conv = vec![0.0_f32; width * height];
    super::convolution::separable_convolve(data, width, height, &kernel_1d, &mut conv);

    // Peak detection: conv > threshold AND conv > all 8 neighbors
    // When noise_map is available, use local adaptive threshold per pixel.
    let mut peaks: Vec<(usize, usize, f32)> = Vec::new();
    for y in (radius + 1)..(height - radius - 1) {
        for x in (radius + 1)..(width - radius - 1) {
            let c = conv[y * width + x];
            let local_threshold = if let Some(nm) = noise_map {
                params.detection_sigma * nm[y * width + x] * kernel_energy_sq.sqrt()
            } else {
                threshold
            };
            if c <= local_threshold {
                continue;
            }
            // 8-neighbor comparison + minimum neighbor count above threshold.
            // Requiring ≥3 neighbors above threshold rejects isolated noise peaks
            // while keeping real stars (which have extended PSF wings).  Matches
            // Siril's star_finder.c:329 candidate validation.
            let neighbors = [
                conv[(y - 1) * width + x - 1],
                conv[(y - 1) * width + x],
                conv[(y - 1) * width + x + 1],
                conv[y * width + x - 1],
                conv[y * width + x + 1],
                conv[(y + 1) * width + x - 1],
                conv[(y + 1) * width + x],
                conv[(y + 1) * width + x + 1],
            ];
            let is_max = neighbors.iter().all(|&v| c > v);
            if is_max {
                let above = neighbors.iter().filter(|&&v| v > local_threshold).count();
                if above >= 3 {
                    peaks.push((x, y, c));
                }
            }
        }
    }

    // Non-maximum suppression using spatial hashing.
    // Suppression radius = measurement stamp radius (max(8, 4σ+1)) so sidelobes
    // are removed before blend detection (which rejects both peaks unconditionally).
    peaks.sort_by(|a, b| b.2.total_cmp(&a.2));
    let nms_sigma = fwhm / 2.3548;
    let nms_blend_r = (4.0 * nms_sigma + 2.0).max(9.0).ceil() as usize;
    let sup_radius = radius.max(nms_blend_r);
    let sup_radius_sq = (sup_radius * sup_radius) as f32;
    let cell_size = sup_radius.max(1);

    // Grid dimensions
    let grid_w = (width + cell_size - 1) / cell_size;
    let grid_h = (height + cell_size - 1) / cell_size;
    let mut grid: Vec<Vec<usize>> = vec![Vec::new(); grid_w * grid_h];

    let mut peak_positions: Vec<(usize, usize, f32)> = Vec::new();

    for (i, &(px, py, conv_val)) in peaks.iter().enumerate() {
        let gx = px / cell_size;
        let gy = py / cell_size;

        // Check neighboring cells for suppression
        let gx_lo = gx.saturating_sub(1);
        let gy_lo = gy.saturating_sub(1);
        let gx_hi = (gx + 2).min(grid_w);
        let gy_hi = (gy + 2).min(grid_h);

        let mut is_suppressed = false;
        'outer: for ngy in gy_lo..gy_hi {
            for ngx in gx_lo..gx_hi {
                for &kept_idx in &grid[ngy * grid_w + ngx] {
                    let (kx, ky, _) = peaks[kept_idx];
                    let dx = px as f32 - kx as f32;
                    let dy = py as f32 - ky as f32;
                    if dx * dx + dy * dy <= sup_radius_sq {
                        is_suppressed = true;
                        break 'outer;
                    }
                }
            }
        }

        if !is_suppressed {
            grid[gy * grid_w + gx].push(i);
            peak_positions.push((px, py, conv_val));
        }
    }

    // ── Stage 2: Proximity-based blend rejection + stamp metrics ────────
    // Two peaks within blend_radius share overlapping PSFs — skip both.
    // blend_radius must cover the measurement stamp radius (max(8, 4σ) in metrics.rs)
    // so that no neighbor falls inside a star's fitting region.
    let sigma = fwhm / 2.3548;
    let meas_stamp_r = (4.0 * sigma + 2.0).max(9.0);
    let blend_radius = meas_stamp_r.max(2.0 * fwhm);
    let blend_radius_sq = blend_radius * blend_radius;
    let blend_cell = blend_radius.ceil().max(1.0) as usize;
    let blend_grid_w = (width + blend_cell - 1) / blend_cell;
    let blend_grid_h = (height + blend_cell - 1) / blend_cell;

    let mut blend_grid: Vec<Vec<usize>> = vec![Vec::new(); blend_grid_w * blend_grid_h];
    for (i, &(px, py, _)) in peak_positions.iter().enumerate() {
        let gx = (px / blend_cell).min(blend_grid_w - 1);
        let gy = (py / blend_cell).min(blend_grid_h - 1);
        blend_grid[gy * blend_grid_w + gx].push(i);
    }

    // Mark blended peaks: any two peaks within blend_radius reject both.
    // NMS already handles sidelobes — after NMS, surviving close pairs are
    // genuine separate stars whose overlapping PSFs corrupt measurements.
    let mut blended = vec![false; peak_positions.len()];
    for (i, &(px, py, _)) in peak_positions.iter().enumerate() {
        if blended[i] { continue; }
        let gx = px / blend_cell;
        let gy = py / blend_cell;
        let search_r = (blend_radius / blend_cell as f32).ceil() as usize + 1;
        let gx_lo = gx.saturating_sub(search_r);
        let gy_lo = gy.saturating_sub(search_r);
        let gx_hi = (gx + search_r + 1).min(blend_grid_w);
        let gy_hi = (gy + search_r + 1).min(blend_grid_h);

        for ngy in gy_lo..gy_hi {
            for ngx in gx_lo..gx_hi {
                for &j in &blend_grid[ngy * blend_grid_w + ngx] {
                    if j <= i { continue; }
                    let (jx, jy, _) = peak_positions[j];
                    let dx = px as f32 - jx as f32;
                    let dy = py as f32 - jy as f32;
                    if dx * dx + dy * dy <= blend_radius_sq {
                        blended[i] = true;
                        blended[j] = true;
                    }
                }
            }
        }
    }

    // Process non-blended peaks via stamp-based metrics.
    // Stamp radius = 1×FWHM (smaller than blend_radius to avoid neighbor contamination).
    let stamp_r = (fwhm.ceil() as i32).max(3);
    let mut stars = Vec::new();
    for (i, &(px, py, _conv_val)) in peak_positions.iter().enumerate() {
        if blended[i] { continue; }
        if let Some(star) = process_peak_stamp(
            px, py, stamp_r, data, width, height, background, bg_map, params,
        ) {
            // ── Pass 2 filters (when field_fwhm is provided) ──
            if let Some(ff) = field_fwhm {
                // 1. Sharpness
                if px >= 1 && py >= 1 && px < width - 1 && py < height - 1 {
                    let cp = conv[py * width + px];
                    if cp > 0.0 {
                        let neighbors_sum =
                            conv[(py - 1) * width + px - 1]
                            + conv[(py - 1) * width + px]
                            + conv[(py - 1) * width + px + 1]
                            + conv[py * width + px - 1]
                            + conv[py * width + px + 1]
                            + conv[(py + 1) * width + px - 1]
                            + conv[(py + 1) * width + px]
                            + conv[(py + 1) * width + px + 1];
                        let sharpness = (cp - neighbors_sum / 8.0) / cp;
                        if sharpness < 0.1 || sharpness > 0.9 {
                            continue;
                        }
                    }
                }

                // 2. Concentration index
                let field_sigma = ff / 2.3548;
                let ci_inner_r_sq = field_sigma * field_sigma;
                let ci_outer_r = 3.0 * field_sigma;
                let ci_outer_r_sq = ci_outer_r * ci_outer_r;
                let ci_radius = ci_outer_r.ceil() as i32;
                let scx = star.x.round() as i32;
                let scy = star.y.round() as i32;

                let mut flux_inner = 0.0_f64;
                let mut flux_outer = 0.0_f64;
                for dy in -ci_radius..=ci_radius {
                    let py2 = scy + dy;
                    if py2 < 0 || py2 >= height as i32 { continue; }
                    for dx in -ci_radius..=ci_radius {
                        let px2 = scx + dx;
                        if px2 < 0 || px2 >= width as i32 { continue; }
                        let r_sq = (dx * dx + dy * dy) as f32;
                        if r_sq > ci_outer_r_sq { continue; }
                        let bg = bg_map.map_or(background, |m| m[py2 as usize * width + px2 as usize]);
                        let val = (data[py2 as usize * width + px2 as usize] - bg).max(0.0) as f64;
                        flux_outer += val;
                        if r_sq <= ci_inner_r_sq {
                            flux_inner += val;
                        }
                    }
                }
                if flux_outer > 0.0 {
                    let ci = flux_inner / flux_outer;
                    if ci < 0.15 || ci > 0.75 {
                        continue;
                    }
                }

                // 3. Edge margin
                let margin = (2.0 * ff).max(8.0);
                if star.x < margin || star.y < margin
                    || star.x > (width as f32 - margin)
                    || star.y > (height as f32 - margin)
                {
                    continue;
                }
            }

            stars.push(star);
        }
    }

    // Post-detection centroid dedup: belt-and-suspenders safety net catches any
    // remaining duplicates from NMS edge cases.
    let dedup_radius = (fwhm * 0.5).max(1.5_f32);
    let dedup_radius_sq = dedup_radius * dedup_radius;
    stars.sort_by(|a, b| b.flux.total_cmp(&a.flux));
    let mut deduped: Vec<DetectedStar> = Vec::with_capacity(stars.len());
    for star in stars.drain(..) {
        let is_dup = deduped.iter().any(|k| {
            let dx = star.x - k.x;
            let dy = star.y - k.y;
            dx * dx + dy * dy <= dedup_radius_sq
        });
        if !is_dup {
            deduped.push(star);
        }
    }
    stars = deduped;

    stars
}

/// Compute star metrics from a stamp around a peak position.
fn process_peak_stamp(
    peak_x: usize,
    peak_y: usize,
    stamp_r: i32,
    data: &[f32],
    width: usize,
    height: usize,
    background: f32,
    bg_map: Option<&[f32]>,
    params: &DetectionParams,
) -> Option<DetectedStar> {
    let cx_i = peak_x as i32;
    let cy_i = peak_y as i32;

    if cx_i - stamp_r <= 0 || cy_i - stamp_r <= 0
        || cx_i + stamp_r >= width as i32 - 1
        || cy_i + stamp_r >= height as i32 - 1
    {
        return None;
    }

    let bg_at = |x: usize, y: usize| bg_map.map_or(background, |m| m[y * width + x]);

    // First pass: find peak value in stamp
    let mut peak = 0.0_f32;
    let mut raw_peak = 0.0_f32;
    for dy in -stamp_r..=stamp_r {
        let py = (cy_i + dy) as usize;
        for dx in -stamp_r..=stamp_r {
            let px = (cx_i + dx) as usize;
            let raw = data[py * width + px];
            let bg = bg_at(px, py);
            let val = raw - bg;
            if val > peak { peak = val; }
            if raw > raw_peak { raw_peak = raw; }
        }
    }

    if peak <= 0.0 { return None; }

    // Second pass: compute metrics using 5% of peak threshold.
    // Only includes pixels that clearly belong to THIS star, not neighbors.
    let threshold = peak * 0.05;
    let mut sum_w = 0.0_f64;
    let mut sum_wx = 0.0_f64;
    let mut sum_wy = 0.0_f64;
    let mut flux = 0.0_f64;
    let mut area = 0_usize;

    for dy in -stamp_r..=stamp_r {
        let py = (cy_i + dy) as usize;
        for dx in -stamp_r..=stamp_r {
            let px = (cx_i + dx) as usize;
            let bg = bg_at(px, py);
            let val = data[py * width + px] - bg;
            if val <= threshold { continue; }

            area += 1;
            let w = (val as f64).powi(2);
            sum_w += w;
            sum_wx += w * px as f64;
            sum_wy += w * py as f64;
            flux += val as f64;
        }
    }

    if area < params.min_star_area || area > params.max_star_area {
        return None;
    }
    if raw_peak > params.saturation_limit {
        return None;
    }
    if sum_w < 1e-10 {
        return None;
    }

    let cx = (sum_wx / sum_w) as f32;
    let cy = (sum_wy / sum_w) as f32;

    // Stamp-based I-weighted second moments for theta and eccentricity
    let mut sf = 0.0_f64;
    let mut six = 0.0_f64;
    let mut siy = 0.0_f64;
    let mut sixx = 0.0_f64;
    let mut siyy = 0.0_f64;
    let mut sixy = 0.0_f64;
    for dy in -stamp_r..=stamp_r {
        let py = (cy_i + dy) as usize;
        for dx in -stamp_r..=stamp_r {
            let px = (cx_i + dx) as usize;
            let bg = bg_at(px, py);
            let v = data[py * width + px] - bg;
            if v <= threshold { continue; }
            let v = v as f64;
            sf += v;
            six += v * px as f64;
            siy += v * py as f64;
            sixx += v * (px as f64) * (px as f64);
            siyy += v * (py as f64) * (py as f64);
            sixy += v * (px as f64) * (py as f64);
        }
    }
    let (theta, ecc) = if sf > 1e-10 {
        let icx = six / sf;
        let icy = siy / sf;
        let mxx = sixx / sf - icx * icx;
        let myy = siyy / sf - icy * icy;
        let mxy = sixy / sf - icx * icy;
        let t = (0.5 * (2.0 * mxy).atan2(mxx - myy)) as f32;
        let trace = mxx + myy;
        let det = mxx * myy - mxy * mxy;
        let disc = (trace * trace - 4.0 * det).max(0.0);
        let l1 = (trace + disc.sqrt()) * 0.5;
        let l2 = (trace - disc.sqrt()) * 0.5;
        let e = if l1 > 0.0 { (1.0 - l2 / l1).max(0.0).sqrt() as f32 } else { 0.0 };
        (t, e)
    } else {
        (0.0, 0.0)
    };

    Some(DetectedStar {
        x: cx,
        y: cy,
        peak,
        flux: flux as f32,
        area,
        theta,
        eccentricity: ecc,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Generate a synthetic star field for testing.
    fn make_star_field(
        width: usize,
        height: usize,
        stars: &[(f32, f32, f32, f32)], // (x, y, amplitude, sigma)
        background: f32,
        noise_sigma: f32,
    ) -> Vec<f32> {
        let mut data = vec![background; width * height];

        // Add stars (Gaussian profiles)
        for &(sx, sy, amp, sigma) in stars {
            let r = (4.0 * sigma).ceil() as i32;
            let inv_2s2 = 1.0 / (2.0 * sigma * sigma);
            for dy in -r..=r {
                for dx in -r..=r {
                    let px = sx as i32 + dx;
                    let py = sy as i32 + dy;
                    if px >= 0 && px < width as i32 && py >= 0 && py < height as i32 {
                        let ddx = px as f32 - sx;
                        let ddy = py as f32 - sy;
                        data[py as usize * width + px as usize] +=
                            amp * (-inv_2s2 * (ddx * ddx + ddy * ddy)).exp();
                    }
                }
            }
        }

        // Add pseudo-noise
        if noise_sigma > 0.0 {
            let mut rng = 99u64;
            for val in data.iter_mut() {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                let u1 = ((rng >> 11) as f64 / (1u64 << 53) as f64).max(1e-15);
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                let u2 = (rng >> 11) as f64 / (1u64 << 53) as f64;
                let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
                *val += noise_sigma * z as f32;
            }
        }

        data
    }

    #[test]
    fn test_detect_synthetic_stars() {
        let width = 200;
        let height = 200;
        let star_defs = vec![
            (50.0, 50.0, 5000.0, 3.0),
            (100.0, 80.0, 3000.0, 3.0),
            (150.0, 120.0, 7000.0, 3.0),
            (30.0, 160.0, 4000.0, 3.0),
            (170.0, 40.0, 6000.0, 3.0),
        ];
        let background = 1000.0;
        let noise = 50.0;

        let data = make_star_field(width, height, &star_defs, background, noise);

        let params = DetectionParams {
            detection_sigma: 5.0,
            min_star_area: 5,
            max_star_area: 2000,
            saturation_limit: 0.95 * 65535.0,
        };

        let fwhm = 3.0 * 2.3548; // sigma=3.0 → fwhm≈7.06
        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, fwhm, None);

        // Should detect all 5 stars
        assert!(
            stars.len() >= 4,
            "Expected at least 4 stars, got {}",
            stars.len()
        );
        assert!(
            stars.len() <= 10,
            "Expected at most 10 detections, got {}",
            stars.len()
        );

        // Check centroids are near truth (within 1 pixel)
        for &(sx, sy, _, _) in &star_defs {
            let closest = stars
                .iter()
                .map(|s| {
                    let dx = s.x - sx;
                    let dy = s.y - sy;
                    (dx * dx + dy * dy).sqrt()
                })
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap();
            assert!(
                closest < 1.5,
                "Star at ({}, {}) not found within 1.5px (closest={})",
                sx,
                sy,
                closest
            );
        }
    }

    #[test]
    fn test_reject_hot_pixels() {
        let width = 100;
        let height = 100;
        let background = 1000.0;
        let noise = 50.0;

        let mut data = vec![background; width * height];
        // Single hot pixel (should be rejected by min_star_area)
        data[50 * width + 50] = 10000.0;

        // Real star
        let sigma = 3.0_f32;
        let inv_2s2 = 1.0 / (2.0 * sigma * sigma);
        for dy in -10i32..=10 {
            for dx in -10i32..=10 {
                let px = (30 + dx) as usize;
                let py = (30 + dy) as usize;
                if px < width && py < height {
                    data[py * width + px] +=
                        5000.0 * (-inv_2s2 * (dx as f32 * dx as f32 + dy as f32 * dy as f32)).exp();
                }
            }
        }

        let params = DetectionParams {
            min_star_area: 5,
            ..DetectionParams::default()
        };

        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, 3.0, None);

        // Should detect the real star but not the hot pixel
        assert!(
            stars.len() >= 1,
            "Should detect at least 1 star"
        );

        // No detection should be at (50, 50)
        for s in &stars {
            let dist = ((s.x - 50.0).powi(2) + (s.y - 50.0).powi(2)).sqrt();
            assert!(
                dist > 3.0,
                "Hot pixel at (50,50) should have been rejected, got star at ({}, {})",
                s.x,
                s.y
            );
        }
    }

    #[test]
    fn test_skip_blended_close_stars() {
        // Two stars 8px apart with sigma=1.5 — their wings overlap at 1.5σ threshold,
        // merging into one CCL component with two peaks. With the PI approach, the
        // multi-peak component is skipped entirely (blended profiles are unreliable).
        let width = 100;
        let height = 100;
        let background = 1000.0;
        let noise = 30.0;
        let star_defs = vec![
            (40.0, 50.0, 8000.0, 1.5),
            (48.0, 50.0, 8000.0, 1.5),
        ];

        let data = make_star_field(width, height, &star_defs, background, noise);

        let params = DetectionParams {
            detection_sigma: 5.0,
            min_star_area: 3,
            max_star_area: 2000,
            saturation_limit: 0.95 * 65535.0,
        };

        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, 3.0, None);

        // With multi-peak skip, blended pair may produce 0 detections (if CCL merges
        // them) or 2 (if CCL keeps them separate). Either way: no duplicates.
        let mut dup_count = 0;
        for (i, a) in stars.iter().enumerate() {
            for b in stars.iter().skip(i + 1) {
                let dx = a.x - b.x;
                let dy = a.y - b.y;
                if (dx * dx + dy * dy).sqrt() < 3.0 {
                    dup_count += 1;
                }
            }
        }
        assert_eq!(
            dup_count, 0,
            "No duplicate detections expected, got {} close pairs (centroids: {:?})",
            dup_count,
            stars.iter().map(|s| (s.x, s.y)).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_no_duplicate_detections() {
        // Single bright star — must produce exactly 1 detection, never a duplicate.
        let width = 100;
        let height = 100;
        let background = 1000.0;
        let noise = 30.0;
        let star_defs = vec![
            (50.0, 50.0, 10000.0, 2.5),
        ];

        let data = make_star_field(width, height, &star_defs, background, noise);

        let params = DetectionParams {
            detection_sigma: 5.0,
            min_star_area: 3,
            max_star_area: 2000,
            saturation_limit: 0.95 * 65535.0,
        };

        let fwhm = 2.5 * 2.3548; // sigma=2.5
        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, fwhm, None);

        // Count detections within one FWHM of the true position
        let near: Vec<_> = stars
            .iter()
            .filter(|s| {
                let dx = s.x - 50.0;
                let dy = s.y - 50.0;
                (dx * dx + dy * dy).sqrt() < fwhm
            })
            .collect();

        assert_eq!(
            near.len(), 1,
            "Expected exactly 1 detection near (50,50), got {} (centroids: {:?})",
            near.len(),
            stars.iter().map(|s| (s.x, s.y)).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_no_deblend_extended_object() {
        // Large extended Gaussian "coma" at center + 5 real point-source stars
        // in the dark sky. The coma should NOT be deblended into false stars.
        let width = 300;
        let height = 300;
        let background = 1000.0;
        let noise = 30.0;

        // Point-source stars well away from center
        let star_defs = vec![
            (30.0, 30.0, 5000.0, 1.5),
            (270.0, 30.0, 6000.0, 1.5),
            (30.0, 270.0, 4000.0, 1.5),
            (270.0, 270.0, 7000.0, 1.5),
            (150.0, 30.0, 5500.0, 1.5),
        ];

        let mut data = make_star_field(width, height, &star_defs, background, noise);

        // Add a large extended coma (sigma=20) at the center — well above max_star_area
        let coma_x = 150.0_f32;
        let coma_y = 150.0_f32;
        let coma_amp = 3000.0_f32;
        let coma_sigma = 20.0_f32;
        let coma_r = (4.0 * coma_sigma).ceil() as i32;
        let inv_2s2 = 1.0 / (2.0 * coma_sigma * coma_sigma);
        for dy in -coma_r..=coma_r {
            for dx in -coma_r..=coma_r {
                let px = coma_x as i32 + dx;
                let py = coma_y as i32 + dy;
                if px >= 0 && px < width as i32 && py >= 0 && py < height as i32 {
                    let ddx = px as f32 - coma_x;
                    let ddy = py as f32 - coma_y;
                    data[py as usize * width + px as usize] +=
                        coma_amp * (-inv_2s2 * (ddx * ddx + ddy * ddy)).exp();
                }
            }
        }

        let params = DetectionParams {
            detection_sigma: 5.0,
            min_star_area: 5,
            max_star_area: 2000,
            saturation_limit: 0.95 * 65535.0,
        };

        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, 3.0, None);

        // All 5 real stars should be detected
        for &(sx, sy, _, _) in &star_defs {
            let closest = stars
                .iter()
                .map(|s| {
                    let dx = s.x - sx;
                    let dy = s.y - sy;
                    (dx * dx + dy * dy).sqrt()
                })
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap();
            assert!(
                closest < 2.0,
                "Star at ({}, {}) not found within 2px (closest={})",
                sx, sy, closest
            );
        }

        // No false detections in the coma region (within 40px of center)
        let coma_detections: Vec<_> = stars
            .iter()
            .filter(|s| {
                let dx = s.x - coma_x;
                let dy = s.y - coma_y;
                (dx * dx + dy * dy).sqrt() < 40.0
            })
            .collect();
        assert!(
            coma_detections.is_empty(),
            "Expected no detections in coma region, got {} (positions: {:?})",
            coma_detections.len(),
            coma_detections.iter().map(|s| (s.x, s.y)).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_sharpness_rejects_diffuse_blob() {
        let width = 200;
        let height = 200;
        let background = 1000.0;
        let noise = 30.0;
        let star_defs = vec![
            (50.0, 50.0, 5000.0, 2.0),   // real star
            (130.0, 130.0, 3000.0, 15.0), // diffuse blob
        ];
        let data = make_star_field(width, height, &star_defs, background, noise);
        let params = DetectionParams {
            detection_sigma: 3.0,
            min_star_area: 5,
            max_star_area: 2000,
            saturation_limit: 0.95 * 65535.0,
        };
        let field_fwhm = 2.0 * 2.3548;
        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, 3.0, Some(field_fwhm));
        let has_real = stars.iter().any(|s| {
            let dx = s.x - 50.0;
            let dy = s.y - 50.0;
            (dx * dx + dy * dy).sqrt() < 3.0
        });
        assert!(has_real, "Real star at (50,50) should be detected");
        let has_blob = stars.iter().any(|s| {
            let dx = s.x - 130.0;
            let dy = s.y - 130.0;
            (dx * dx + dy * dy).sqrt() < 20.0
        });
        assert!(!has_blob, "Diffuse blob should be rejected by sharpness/CI filters");
    }

    #[test]
    fn test_edge_margin_rejects_near_border_stars() {
        let width = 200;
        let height = 200;
        let background = 1000.0;
        let noise = 30.0;
        // Use sigma=1.5 stars (sharp point sources) so sharpness/CI filters pass.
        // field_fwhm = 1.5 * 2.3548 ≈ 3.53; edge margin = max(2 * 3.53, 8) = 8.0.
        // Edge star at x=5 is within margin=8, center star at x=100 is safely inside.
        let field_fwhm = 1.5 * 2.3548;
        let star_defs = vec![
            (5.0, 100.0, 8000.0, 1.5),  // near left edge (x < margin=8)
            (100.0, 100.0, 8000.0, 1.5), // safely inside
        ];
        let data = make_star_field(width, height, &star_defs, background, noise);
        let params = DetectionParams {
            detection_sigma: 3.0,
            ..DetectionParams::default()
        };
        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, field_fwhm, Some(field_fwhm));
        let has_center = stars.iter().any(|s| {
            let dx = s.x - 100.0;
            let dy = s.y - 100.0;
            (dx * dx + dy * dy).sqrt() < 3.0
        });
        assert!(has_center, "Center star at (100,100) should be detected");
        let has_edge = stars.iter().any(|s| {
            let dx = s.x - 5.0;
            let dy = s.y - 100.0;
            (dx * dx + dy * dy).sqrt() < 3.0
        });
        assert!(!has_edge, "Edge star at (5,100) should be rejected by edge margin");
    }

    #[test]
    fn test_pass1_no_filtering() {
        // With field_fwhm=None (Pass 1), no sharpness/CI/edge filters apply
        let width = 200;
        let height = 200;
        let background = 1000.0;
        let noise = 30.0;
        let star_defs = vec![
            (50.0, 50.0, 5000.0, 2.0),
            (10.0, 50.0, 5000.0, 2.0),  // near edge — should NOT be rejected in Pass 1
        ];
        let data = make_star_field(width, height, &star_defs, background, noise);
        let params = DetectionParams { detection_sigma: 3.0, ..DetectionParams::default() };
        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, 4.7, None);
        // Both stars should be detected in Pass 1
        assert!(stars.len() >= 2, "Pass 1 (None) should not apply edge/sharpness filters, got {} stars", stars.len());
    }
}
