/// Star detection: DAOFIND-inspired matched filter + connected component labeling.

/// A detected star candidate before metric computation.
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

/// Detect stars in an image using DAOFIND-style matched filter + CCL.
///
/// `data`: single-channel f32 image (raw ADU values, NOT background-subtracted).
/// `background`: global background level.
/// `noise`: background noise sigma.
/// `bg_map`: optional per-pixel background map (from mesh-grid estimation).
/// `noise_map`: optional per-pixel noise map for adaptive thresholds.
/// `fwhm`: estimated FWHM for matched filter kernel (pixels).
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
    let low_threshold_global = 1.5 * noise; // Lower threshold to capture star wings

    // Build binary mask and label via two-pass CCL
    let mut labels = vec![0u32; width * height];
    let mut parent = vec![0u32; 1]; // Index 0 unused (label 0 = background)
    let mut next_label = 1u32;

    // First pass
    for y in 0..height {
        for x in 0..width {
            let bg = bg_map.map_or(background, |m| m[y * width + x]);
            let val = data[y * width + x] - bg;
            let low_threshold = if let Some(nm) = noise_map {
                1.5 * nm[y * width + x]
            } else {
                low_threshold_global
            };
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
    // Build map of label → peak positions (for deblending multi-peak components)
    use std::collections::HashMap;
    let mut peak_map: HashMap<u32, Vec<(usize, usize)>> = HashMap::new();
    for &(px, py) in &peak_positions {
        let l = labels[py * width + px];
        if l > 0 {
            peak_map.entry(l).or_default().push((px, py));
        } else {
            // Peak might be just off — check neighbors, add to first found label
            'search: for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nx = px as i32 + dx;
                    let ny = py as i32 + dy;
                    if nx >= 0 && nx < width as i32 && ny >= 0 && ny < height as i32 {
                        let l2 = labels[ny as usize * width + nx as usize];
                        if l2 > 0 {
                            peak_map.entry(l2).or_default().push((px, py));
                            break 'search;
                        }
                    }
                }
            }
        }
    }

    // Gather pixel lists per component
    let mut components: HashMap<u32, Vec<(usize, usize)>> = HashMap::new();
    for y in 0..height {
        for x in 0..width {
            let l = labels[y * width + x];
            if l > 0 && peak_map.contains_key(&l) {
                components.entry(l).or_default().push((x, y));
            }
        }
    }

    // Validate components, deblend multi-peak blobs, and compute centroids
    let mut stars = Vec::new();
    for (label, pixels) in &components {
        let peaks = &peak_map[label];
        if peaks.len() <= 1 || pixels.len() > params.max_star_area {
            // Single peak, or component too large to be merged stars: process as whole.
            // Extended objects (comet comae, nebulae) create huge CCL blobs with many
            // convolution peaks — deblending would produce false detections.
            if let Some(star) = process_component(pixels, data, width, height, background, bg_map, params) {
                stars.push(star);
            }
        } else {
            // Multi-peak with reasonable per-peak size: split by nearest peak (Voronoi tessellation)
            let mut subs: Vec<Vec<(usize, usize)>> = vec![Vec::new(); peaks.len()];
            for &(x, y) in pixels {
                let nearest = peaks
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, &(px, py))| {
                        let dx = x as i32 - px as i32;
                        let dy = y as i32 - py as i32;
                        (dx * dx + dy * dy) as u32
                    })
                    .unwrap()
                    .0;
                subs[nearest].push((x, y));
            }
            for sub in &subs {
                if let Some(star) = process_component(sub, data, width, height, background, bg_map, params) {
                    stars.push(star);
                }
            }
        }
    }

    // Sort by flux descending (useful for two-pass FWHM estimation on top-20 brightest)
    stars.sort_by(|a, b| b.flux.total_cmp(&a.flux));

    stars
}

// ── Per-component processing ────────────────────────────────────────────────

/// Validate a component (pixel list) and compute centroid / shape metrics.
/// Returns `None` if the component is rejected (area, border, saturation, aspect ratio).
fn process_component(
    pixels: &[(usize, usize)],
    data: &[f32],
    width: usize,
    height: usize,
    background: f32,
    bg_map: Option<&[f32]>,
    params: &DetectionParams,
) -> Option<DetectedStar> {
    let area = pixels.len();

    // Area filter
    if area < params.min_star_area || area > params.max_star_area {
        return None;
    }

    // Border rejection
    let touches_border = pixels.iter().any(|&(x, y)| x == 0 || y == 0 || x == width - 1 || y == height - 1);
    if touches_border {
        return None;
    }

    // Compute peak and check saturation
    let mut peak = 0.0_f32;
    let mut raw_peak = 0.0_f32;
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
        let raw = data[y * width + x];
        let bg = bg_map.map_or(background, |m| m[y * width + x]);
        let val = raw - bg;
        if val > peak {
            peak = val;
        }
        if raw > raw_peak {
            raw_peak = raw;
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

    // Saturation check — use raw (pre-background-subtracted) peak
    if raw_peak > params.saturation_limit {
        return None;
    }

    // Aspect ratio check (reject cosmic rays/satellite trails)
    let bbox_w = (max_x - min_x + 1) as f32;
    let bbox_h = (max_y - min_y + 1) as f32;
    let aspect = bbox_w.max(bbox_h) / bbox_w.min(bbox_h);
    if aspect > 4.0 {
        return None;
    }

    if sum_w < 1e-10 {
        return None;
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
        };

        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, 3.0);

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

        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, 3.0);

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
    fn test_deblend_close_stars() {
        // Two stars 8px apart with sigma=1.5 — their wings overlap at 1.5σ threshold,
        // merging into one CCL component. Deblending should split them.
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

        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, 3.0);

        assert!(
            stars.len() >= 2,
            "Expected at least 2 deblended stars, got {} (centroids: {:?})",
            stars.len(),
            stars.iter().map(|s| (s.x, s.y)).collect::<Vec<_>>()
        );

        // Each true star should have a detection within 2px
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

        let stars = detect_stars(&data, width, height, background, noise, None, None, &params, 3.0);

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
}
