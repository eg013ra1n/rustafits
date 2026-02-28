/// Star detection: DAOFIND-inspired matched filter + connected component labeling.

/// A detected star candidate before metric computation.
pub(crate) struct DetectedStar {
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
pub(crate) struct DetectionParams {
    pub detection_sigma: f32,
    pub min_star_area: usize,
    pub max_star_area: usize,
    pub saturation_limit: f32,
    pub max_stars: usize,
}

impl Default for DetectionParams {
    fn default() -> Self {
        DetectionParams {
            detection_sigma: 5.0,
            min_star_area: 5,
            max_star_area: 2000,
            saturation_limit: 0.95 * 65535.0,
            max_stars: 200,
        }
    }
}

/// Detect stars in an image using DAOFIND-style matched filter + CCL.
///
/// `data`: single-channel f32 image (raw ADU values, NOT background-subtracted).
/// `background`: global background level.
/// `noise`: background noise sigma.
/// `bg_map`: optional per-pixel background map (from mesh-grid estimation).
pub(crate) fn detect_stars(
    data: &[f32],
    width: usize,
    height: usize,
    background: f32,
    noise: f32,
    bg_map: Option<&[f32]>,
    params: &DetectionParams,
) -> Vec<DetectedStar> {
    // Stage 1: Gaussian convolution + peak detection
    let estimated_fwhm = 3.0_f32; // Default initial estimate
    let sigma = estimated_fwhm / 2.3548;
    let radius = (2.0 * sigma).ceil() as usize;
    let ksize = 2 * radius + 1;

    // Generate zero-sum Gaussian kernel (DAOFIND key insight)
    let mut kernel = vec![0.0_f32; ksize * ksize];
    let inv_2s2 = 1.0 / (2.0 * sigma * sigma);
    let mut ksum = 0.0_f32;
    for ky in 0..ksize {
        for kx in 0..ksize {
            let dy = ky as f32 - radius as f32;
            let dx = kx as f32 - radius as f32;
            let g = (-inv_2s2 * (dx * dx + dy * dy)).exp();
            kernel[ky * ksize + kx] = g;
            ksum += g;
        }
    }
    // Subtract mean to make zero-sum
    let kmean = ksum / (ksize * ksize) as f32;
    let mut kernel_energy_sq = 0.0_f32;
    for v in kernel.iter_mut() {
        *v -= kmean;
        kernel_energy_sq += *v * *v;
    }

    // Convolution threshold: detection_sigma × noise × sqrt(Σ K²)
    let threshold = params.detection_sigma * noise * kernel_energy_sq.sqrt();

    // Convolve image with kernel and find local maxima
    let mut conv = vec![0.0_f32; width * height];
    for y in radius..(height - radius) {
        for x in radius..(width - radius) {
            let mut sum = 0.0_f32;
            for ky in 0..ksize {
                let iy = y + ky - radius;
                let row_off = iy * width;
                let k_row_off = ky * ksize;
                for kx in 0..ksize {
                    let ix = x + kx - radius;
                    sum += data[row_off + ix] * kernel[k_row_off + kx];
                }
            }
            conv[y * width + x] = sum;
        }
    }

    // Peak detection: conv > threshold AND conv > all 8 neighbors
    let mut peaks: Vec<(usize, usize, f32)> = Vec::new();
    for y in (radius + 1)..(height - radius - 1) {
        for x in (radius + 1)..(width - radius - 1) {
            let c = conv[y * width + x];
            if c <= threshold {
                continue;
            }
            // 8-neighbor comparison
            let is_max = c > conv[(y - 1) * width + x - 1]
                && c > conv[(y - 1) * width + x]
                && c > conv[(y - 1) * width + x + 1]
                && c > conv[y * width + x - 1]
                && c > conv[y * width + x + 1]
                && c > conv[(y + 1) * width + x - 1]
                && c > conv[(y + 1) * width + x]
                && c > conv[(y + 1) * width + x + 1];
            if is_max {
                peaks.push((x, y, c));
            }
        }
    }

    // Non-maximum suppression using spatial hashing
    peaks.sort_by(|a, b| b.2.total_cmp(&a.2));
    let sup_radius = radius;
    let sup_radius_sq = (sup_radius * sup_radius) as f32;
    let cell_size = sup_radius.max(1);

    // Grid dimensions
    let grid_w = (width + cell_size - 1) / cell_size;
    let grid_h = (height + cell_size - 1) / cell_size;
    let mut grid: Vec<Vec<usize>> = vec![Vec::new(); grid_w * grid_h];

    let mut peak_positions: Vec<(usize, usize)> = Vec::new();

    for (i, &(px, py, _)) in peaks.iter().enumerate() {
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
                    if dx * dx + dy * dy < sup_radius_sq {
                        is_suppressed = true;
                        break 'outer;
                    }
                }
            }
        }

        if !is_suppressed {
            grid[gy * grid_w + gx].push(i);
            peak_positions.push((px, py));
        }
    }

    // Stage 2: CCL with Union-Find on thresholded original image
    let low_threshold = 1.5 * noise; // Lower threshold to capture star wings

    // Build binary mask and label via two-pass CCL
    let mut labels = vec![0u32; width * height];
    let mut parent = vec![0u32; 1]; // Index 0 unused (label 0 = background)
    let mut next_label = 1u32;

    // First pass
    for y in 0..height {
        for x in 0..width {
            let bg = bg_map.map_or(background, |m| m[y * width + x]);
            let val = data[y * width + x] - bg;
            if !val.is_finite() || val <= low_threshold {
                continue;
            }

            let idx = y * width + x;
            let left_label = if x > 0 { labels[idx - 1] } else { 0 };
            let top_label = if y > 0 { labels[idx - width] } else { 0 };

            match (left_label, top_label) {
                (0, 0) => {
                    labels[idx] = next_label;
                    parent.push(next_label);
                    next_label += 1;
                }
                (l, 0) | (0, l) => {
                    labels[idx] = find(&mut parent, l);
                }
                (l, t) => {
                    let rl = find(&mut parent, l);
                    let rt = find(&mut parent, t);
                    labels[idx] = rl.min(rt);
                    if rl != rt {
                        union(&mut parent, rl, rt);
                    }
                }
            }
        }
    }

    // Second pass: resolve labels
    for l in labels.iter_mut() {
        if *l > 0 {
            *l = find(&mut parent, *l);
        }
    }

    // Collect components: only those near a detected peak
    // Build set of peak labels
    let mut peak_labels = std::collections::HashSet::new();
    for &(px, py) in &peak_positions {
        let l = labels[py * width + px];
        if l > 0 {
            peak_labels.insert(l);
        } else {
            // Peak might be just off — check neighbors
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nx = px as i32 + dx;
                    let ny = py as i32 + dy;
                    if nx >= 0 && nx < width as i32 && ny >= 0 && ny < height as i32 {
                        let l2 = labels[ny as usize * width + nx as usize];
                        if l2 > 0 {
                            peak_labels.insert(l2);
                        }
                    }
                }
            }
        }
    }

    // Gather pixel lists per component
    use std::collections::HashMap;
    let mut components: HashMap<u32, Vec<(usize, usize)>> = HashMap::new();
    for y in 0..height {
        for x in 0..width {
            let l = labels[y * width + x];
            if l > 0 && peak_labels.contains(&l) {
                components.entry(l).or_default().push((x, y));
            }
        }
    }

    // Validate components and compute centroids
    let mut stars = Vec::new();
    for (_label, pixels) in &components {
        let area = pixels.len();

        // Area filter
        if area < params.min_star_area || area > params.max_star_area {
            continue;
        }

        // Border rejection
        let touches_border = pixels.iter().any(|&(x, y)| x == 0 || y == 0 || x == width - 1 || y == height - 1);
        if touches_border {
            continue;
        }

        // Compute peak and check saturation
        let mut peak = 0.0_f32;
        let mut sum_w = 0.0_f64;
        let mut sum_wx = 0.0_f64;
        let mut sum_wy = 0.0_f64;
        let mut flux = 0.0_f64;

        // Bounding box for aspect ratio check
        let mut min_x = usize::MAX;
        let mut max_x = 0usize;
        let mut min_y = usize::MAX;
        let mut max_y = 0usize;

        for &(x, y) in pixels {
            let bg = bg_map.map_or(background, |m| m[y * width + x]);
            let val = data[y * width + x] - bg;
            if val > peak {
                peak = val;
            }

            // I² weighting for centroid
            let w = (val.max(0.0) as f64).powi(2);
            sum_w += w;
            sum_wx += w * x as f64;
            sum_wy += w * y as f64;
            flux += val.max(0.0) as f64;

            min_x = min_x.min(x);
            max_x = max_x.max(x);
            min_y = min_y.min(y);
            max_y = max_y.max(y);
        }

        // Saturation check
        if peak > params.saturation_limit {
            continue;
        }

        // Aspect ratio check (reject cosmic rays/satellite trails)
        let bbox_w = (max_x - min_x + 1) as f32;
        let bbox_h = (max_y - min_y + 1) as f32;
        let aspect = bbox_w.max(bbox_h) / bbox_w.min(bbox_h);
        if aspect > 3.0 {
            continue;
        }

        if sum_w < 1e-10 {
            continue;
        }

        let cx = (sum_wx / sum_w) as f32;
        let cy = (sum_wy / sum_w) as f32;

        // Theta from stamp-based I-weighted second moments (not CCL pixels).
        // CCL blobs have too few pixels → grid-induced theta coherence.
        // A stamp over the continuous image gives reliable orientation.
        let stamp_r = ((2.0 * (area as f64 / std::f64::consts::PI).sqrt()) as i32).max(5);
        let cx_i = cx.round() as i32;
        let cy_i = cy.round() as i32;
        let (theta, ecc) = {
            let mut sf = 0.0_f64;
            let mut six = 0.0_f64;
            let mut siy = 0.0_f64;
            let mut sixx = 0.0_f64;
            let mut siyy = 0.0_f64;
            let mut sixy = 0.0_f64;
            for dy in -stamp_r..=stamp_r {
                let py = cy_i + dy;
                if py < 0 || py >= height as i32 { continue; }
                for dx in -stamp_r..=stamp_r {
                    let px = cx_i + dx;
                    if px < 0 || px >= width as i32 { continue; }
                    let bg = bg_map.map_or(background, |m| m[py as usize * width + px as usize]);
                    let v = (data[py as usize * width + px as usize] - bg).max(0.0) as f64;
                    sf += v;
                    six += v * px as f64;
                    siy += v * py as f64;
                    sixx += v * (px as f64) * (px as f64);
                    siyy += v * (py as f64) * (py as f64);
                    sixy += v * (px as f64) * (py as f64);
                }
            }
            if sf > 1e-10 {
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
            }
        };

        stars.push(DetectedStar {
            x: cx,
            y: cy,
            peak,
            flux: flux as f32,
            area,
            theta,
            eccentricity: ecc,
        });
    }

    // Sort by flux descending, keep top max_stars
    stars.sort_by(|a, b| b.flux.total_cmp(&a.flux));
    stars.truncate(params.max_stars);

    stars
}

// ── Union-Find ──────────────────────────────────────────────────────────────

fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        parent[x as usize] = parent[parent[x as usize] as usize]; // path halving
        x = parent[x as usize];
    }
    x
}

fn union(parent: &mut [u32], a: u32, b: u32) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra != rb {
        parent[ra.max(rb) as usize] = ra.min(rb);
    }
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
            max_stars: 200,
        };

        let stars = detect_stars(&data, width, height, background, noise, None, &params);

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

        let stars = detect_stars(&data, width, height, background, noise, None, &params);

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
}
