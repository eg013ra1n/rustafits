//! Adaptive multi-level star detection — the detector used by plate solving.
//!
//! A single global threshold either drowns in the diffuse light of a bright
//! nebula/galaxy or thresholds the diffuse light itself as one giant blob, so
//! it under-detects exactly the long-focal-length / extended-object frames
//! the solver most needs to handle. This detector instead deepens
//! adaptively:
//!
//!  1. **Falling-threshold ladder with an occupancy mask.** Bright stars are
//!     found first at a high, star-count-derived level and their `3·HFD`
//!     disks marked claimed; progressively lower thresholds then add only
//!     fainter, not-yet-claimed stars. The deepest pass recomputes
//!     background/noise on a 12×N mesh, so faint stars sitting *on top of*
//!     nebulosity clear a local threshold instead of a global one.
//!  2. **Saturation-aware level clip** `level = max(3.5σ, level − bg − 1)` so
//!     a flat-topped saturated core still presents supra-threshold pixels
//!     and is detected/centroided rather than vanishing.
//!  3. **MAD-annulus local background + shrink-to-symmetry flux-weighted
//!     centroid**, with a ≥2/4 cross-neighbour hot-pixel gate and
//!     `snr>10`, `hfd∈(hfd_min,30]` accept tests, for clean centroids on a
//!     contaminated/gradient background.

use super::detection::DetectedStar;

/// One detected star (superset of what the solver consumes).
struct Star {
    x: f32,
    y: f32,
    peak: f32,
    flux: f32,
    hfd: f32,
    snr: f32,
}

/// Histogram-mode background + sigma-clipped global noise σ. The histogram
/// mode is robust to the bright extended light of a nebula/galaxy (which a
/// mean would chase); the noise is a 3σ-clipped RMS about that background.
pub(crate) fn background_and_noise(
    lum: &[f32],
    width: usize,
    height: usize,
) -> (f32, f32) {
    // Integer histogram, bins 1..=64999 (ignore 0 and ≥65000).
    let mut hist = vec![0u32; 65000];
    let mut sum = 0.0f64;
    let mut n = 0u64;
    for &v in lum {
        if v <= 0.0 || v >= 65000.0 {
            continue;
        }
        let b = v.round() as usize;
        if b >= 1 && b < 65000 {
            hist[b] += 1;
            sum += v as f64;
            n += 1;
        }
    }
    if n == 0 {
        return (0.0, 1.0);
    }
    let mut mode = 1usize;
    let mut best = 0u32;
    for (b, &c) in hist.iter().enumerate().take(65000).skip(1) {
        if c > best {
            best = c;
            mode = b;
        }
    }
    let mean = (sum / n as f64) as f32;
    // Strange low-value spike → use the mean instead.
    let backgr = if mean > 1.5 * mode as f32 { mean } else { mode as f32 };

    // Sigma-clipped RMS about the background, sparse global sample.
    let mut step = (height as f32 / 71.0).round() as usize;
    if step < 1 {
        step = 1;
    }
    if step % 2 == 0 {
        step += 1;
    }
    let mut sd = 0.0f32;
    for iter in 0..7 {
        let mut acc = 0.0f64;
        let mut cnt = 0u64;
        let mut y = 0;
        while y < height {
            let mut x = 0;
            while x < width {
                let v = lum[y * width + x];
                if v != 0.0 && v < 2.0 * backgr {
                    let d = v - backgr;
                    if iter == 0 || d.abs() <= 3.0 * sd {
                        acc += (d as f64) * (d as f64);
                        cnt += 1;
                    }
                }
                x += step;
            }
            y += step;
        }
        let sd_new = if cnt > 0 {
            (acc / cnt as f64).sqrt() as f32
        } else {
            sd.max(1.0)
        };
        let converged = iter > 0 && (sd - sd_new).abs() < 0.05 * sd_new.max(1e-6);
        sd = sd_new.max(1e-6);
        if converged {
            break;
        }
    }
    (backgr, sd.max(1e-6))
}

/// Count-driven `star_level` / `star_level2` with the saturation-aware clip.
/// Levels are picked so the brightest ~`6·max_stars` / `24·max_stars` pixels
/// fall above them, then pulled below the saturation ceiling (`− bg − 1`) so
/// flat-topped cores still register.
fn star_levels(
    lum: &[f32],
    backgr: f32,
    sd: f32,
    max_stars: usize,
) -> (f32, f32) {
    let mut hist = vec![0u32; 65536];
    for &v in lum {
        if v > 0.0 && v < 65535.0 {
            hist[v.round() as usize] += 1;
        }
    }
    let factor = (6 * max_stars) as u64;
    let factor2 = (24 * max_stars) as u64;
    let mut cum = 0u64;
    let mut level1 = backgr;
    let mut level2 = backgr;
    let mut got1 = false;
    for b in (1..65535).rev() {
        cum += hist[b] as u64;
        if !got1 && cum >= factor {
            level1 = b as f32;
            got1 = true;
        }
        if cum >= factor2 {
            level2 = b as f32;
            break;
        }
    }
    let floor = (3.5 * sd).max(1.0);
    let clip = |lvl: f32| (lvl - backgr - 1.0).max(floor);
    (clip(level1), clip(level2))
}

/// Sigma-clipped (background, noise) over an explicit sample slice.
fn local_bg_noise(samples: &[f32]) -> (f32, f32) {
    if samples.is_empty() {
        return (0.0, 1.0);
    }
    let mut mean = samples.iter().sum::<f32>() / samples.len() as f32;
    let mut sd = 1.0f32;
    for it in 0..5 {
        let mut acc = 0.0f64;
        let mut cnt = 0u64;
        for &v in samples {
            if it == 0 || (v - mean).abs() <= 3.0 * sd {
                acc += (v - mean) as f64 * (v - mean) as f64;
                cnt += 1;
            }
        }
        if cnt == 0 {
            break;
        }
        let sd_new = (acc / cnt as f64).sqrt() as f32;
        // recompute clipped mean
        let mut msum = 0.0f64;
        let mut mcnt = 0u64;
        for &v in samples {
            if (v - mean).abs() <= 3.0 * sd_new.max(1e-6) {
                msum += v as f64;
                mcnt += 1;
            }
        }
        if mcnt > 0 {
            mean = (msum / mcnt as f64) as f32;
        }
        let converged = it > 0 && (sd - sd_new).abs() < 0.1 * sd_new.max(1e-6);
        sd = sd_new.max(1e-6);
        if converged {
            break;
        }
    }
    (mean, sd.max(1e-6))
}

fn median(v: &mut [f32]) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[v.len() / 2]
}

/// Half-flux-diameter + flux-weighted centroid at a seed pixel. Uses a
/// MAD-estimated local annulus background and a shrink-to-symmetry box so
/// the centroid stays clean on a contaminated/gradient background and
/// blended/extended detections are rejected.
/// Returns `(cx, cy, hfd, flux, peak, snr)`.
#[allow(clippy::too_many_arguments)]
fn hfd_at(
    lum: &[f32],
    width: usize,
    height: usize,
    sx: usize,
    sy: usize,
) -> Option<(f32, f32, f32, f32, f32, f32)> {
    let mut rs: i32 = 14;
    let annulus_w: i32 = 3;

    // Local background = median of the rs..rs+annulus_w annulus; σ via MAD.
    let mut ann: Vec<f32> = Vec::new();
    for dy in -(rs + annulus_w)..=(rs + annulus_w) {
        for dx in -(rs + annulus_w)..=(rs + annulus_w) {
            let rr = ((dx * dx + dy * dy) as f32).sqrt();
            if rr > rs as f32 && rr <= (rs + annulus_w) as f32 {
                let x = sx as i32 + dx;
                let y = sy as i32 + dy;
                if x >= 0 && y >= 0 && (x as usize) < width && (y as usize) < height {
                    ann.push(lum[y as usize * width + x as usize]);
                }
            }
        }
    }
    if ann.len() < 8 {
        return None;
    }
    let star_bg = median(&mut ann);
    let mut devs: Vec<f32> = ann.iter().map(|&v| (v - star_bg).abs()).collect();
    let sd_bg = (1.4826 * median(&mut devs)).max(1.0);

    // Flux-weighted centroid + shrink-to-symmetry box.
    let (cx, cy);
    let (mut sum_val, mut sum_x, mut sum_y);
    loop {
        sum_val = 0.0f64;
        sum_x = 0.0f64;
        sum_y = 0.0f64;
        let mut signal = 0u32;
        for dy in -rs..=rs {
            for dx in -rs..=rs {
                let x = sx as i32 + dx;
                let y = sy as i32 + dy;
                if x < 0 || y < 0 || x as usize >= width || y as usize >= height {
                    continue;
                }
                let val = lum[y as usize * width + x as usize] - star_bg;
                if val > 3.0 * sd_bg {
                    sum_val += val as f64;
                    sum_x += val as f64 * x as f64;
                    sum_y += val as f64 * y as f64;
                    signal += 1;
                }
            }
        }
        if sum_val <= 12.0 * sd_bg as f64 || signal <= 1 {
            return None;
        }
        let box_px = ((2 * rs + 1) * (2 * rs + 1)) as f32;
        let boxed = signal as f32 >= (2.0 / 9.0) * box_px;
        if boxed || rs <= 4 {
            cx = (sum_x / sum_val) as f32;
            cy = (sum_y / sum_val) as f32;
            break;
        }
        rs -= if rs <= 4 { 1 } else { 2 };
    }

    // Aperture: grow until the radial signal falls to ≤10% of peak.
    let mut peak = 0.0f32;
    for dy in -rs..=rs {
        for dx in -rs..=rs {
            let x = sx as i32 + dx;
            let y = sy as i32 + dy;
            if x < 0 || y < 0 || x as usize >= width || y as usize >= height {
                continue;
            }
            let val = lum[y as usize * width + x as usize] - star_bg;
            if val > peak {
                peak = val;
            }
        }
    }
    if peak <= 0.0 {
        return None;
    }
    let mut r_ap = 1i32;
    while r_ap < rs {
        // Mean signal in the (r_ap, r_ap+1] annulus.
        let mut s = 0.0f32;
        let mut c = 0u32;
        for dy in -(r_ap + 1)..=(r_ap + 1) {
            for dx in -(r_ap + 1)..=(r_ap + 1) {
                let rr = ((dx * dx + dy * dy) as f32).sqrt();
                if rr > r_ap as f32 && rr <= (r_ap + 1) as f32 {
                    let x = sx as i32 + dx;
                    let y = sy as i32 + dy;
                    if x >= 0
                        && y >= 0
                        && (x as usize) < width
                        && (y as usize) < height
                    {
                        s += lum[y as usize * width + x as usize] - star_bg;
                        c += 1;
                    }
                }
            }
        }
        let m = if c > 0 { s / c as f32 } else { 0.0 };
        if m <= 0.1 * peak {
            break;
        }
        r_ap += 1;
    }

    // HFD (Miyashita) over the r_ap box; flux & snr.
    let mut sum_v = 0.0f64;
    let mut sum_vr = 0.0f64;
    let mut flux = 0.0f64;
    for dy in -r_ap..=r_ap {
        for dx in -r_ap..=r_ap {
            let x = sx as i32 + dx;
            let y = sy as i32 + dy;
            if x < 0 || y < 0 || x as usize >= width || y as usize >= height {
                continue;
            }
            let val = (lum[y as usize * width + x as usize] - star_bg).max(0.0);
            let r = (((x as f32 - cx).powi(2)) + ((y as f32 - cy).powi(2))).sqrt();
            sum_v += val as f64;
            sum_vr += val as f64 * r as f64;
            flux += val as f64;
        }
    }
    if sum_v <= 0.0 {
        return None;
    }
    let hfd = ((2.0 * sum_vr / sum_v) as f32).max(0.8);
    let flux = flux as f32;
    let snr = if flux >= 1.0 {
        flux / (flux + std::f32::consts::PI * (r_ap as f32).powi(2) * sd_bg * sd_bg).sqrt()
    } else {
        0.0
    };
    Some((cx, cy, hfd, flux, peak, snr))
}

/// Stamp a filled disk of radius `r` into the occupancy mask.
fn stamp(mask: &mut [u8], width: usize, height: usize, cx: f32, cy: f32, r: f32) {
    let r = r.max(1.0);
    let ri = r.ceil() as i32;
    let cxi = cx.round() as i32;
    let cyi = cy.round() as i32;
    for dy in -ri..=ri {
        for dx in -ri..=ri {
            if ((dx * dx + dy * dy) as f32) <= r * r {
                let x = cxi + dx;
                let y = cyi + dy;
                if x >= 0 && y >= 0 && (x as usize) < width && (y as usize) < height {
                    mask[y as usize * width + x as usize] = 1;
                }
            }
        }
    }
}

/// Scan a rectangular region at a fixed detection level, appending accepted
/// stars and masking their disks. `bg`/`noise` are the reference background
/// and global noise for the hot-pixel cross test.
#[allow(clippy::too_many_arguments)]
fn scan_region(
    lum: &[f32],
    width: usize,
    height: usize,
    x0: usize,
    y0: usize,
    x1: usize,
    y1: usize,
    bg: f32,
    detection_level: f32,
    noise: f32,
    hfd_min: f32,
    mask: &mut [u8],
    out: &mut Vec<Star>,
    budget: usize,
) {
    let cross = 4.0 * noise;
    for y in y0.max(1)..y1.min(height.saturating_sub(1)) {
        for x in x0.max(1)..x1.min(width.saturating_sub(1)) {
            let i = y * width + x;
            if mask[i] != 0 {
                continue;
            }
            let v = lum[i] - bg;
            if v <= detection_level {
                continue;
            }
            // ≥2 of 4 cross-neighbours above bg+4σ (hot-pixel rejection).
            let n = ((lum[i - 1] - bg > cross) as u32)
                + ((lum[i + 1] - bg > cross) as u32)
                + ((lum[i - width] - bg > cross) as u32)
                + ((lum[i + width] - bg > cross) as u32);
            if n < 2 {
                continue;
            }
            if let Some((cx, cy, hfd, flux, peak, snr)) =
                hfd_at(lum, width, height, x, y)
            {
                if hfd > hfd_min && hfd <= 30.0 && snr > 10.0 {
                    let mi = (cy.round() as i32).clamp(0, height as i32 - 1) as usize
                        * width
                        + (cx.round() as i32).clamp(0, width as i32 - 1) as usize;
                    if mask[mi] == 0 {
                        stamp(mask, width, height, cx, cy, (3.0 * hfd).round());
                        out.push(Star {
                            x: cx,
                            y: cy,
                            peak,
                            flux,
                            hfd,
                            snr,
                        });
                        // Safety bound only — the star cap is applied AFTER a
                        // completed scan (flux sort + truncate). Returning at
                        // max_stars here truncated dense frames in row-major
                        // scan order: cocoon.fits came back as a horizontal
                        // stripe (all 500 stars at y < 1470 of 4176).
                        if out.len() >= budget {
                            return;
                        }
                    }
                }
            }
        }
    }
}

/// Adaptive multi-level detection. Returns stars sorted brightest-first
/// (by flux), trimmed to `max_stars`.
pub fn detect_stars_adaptive(
    lum: &[f32],
    width: usize,
    height: usize,
    max_stars: usize,
    hfd_min: f32,
) -> Vec<(DetectedStar, f32)> {
    if width < 8 || height < 8 {
        return Vec::new();
    }
    let max_stars = max_stars.max(8);
    // Per-scan safety bound: well above any real frame's star count at the
    // ladder levels, only there to bound memory on pathological input.
    let scan_budget = (max_stars * 64).max(8192);
    let (bg, noise) = background_and_noise(lum, width, height);
    let (star_level, star_level2) = star_levels(lum, bg, noise, max_stars);

    let mut mask = vec![0u8; width * height];
    let mut stars: Vec<Star> = Vec::with_capacity(max_stars);

    // retries 4 → 1, stop once we have enough or the ladder is exhausted.
    let mut retries = 4i32;
    while retries >= 1 && stars.len() < max_stars {
        match retries {
            4 => {
                if star_level > 30.0 * noise {
                    scan_region(
                        lum, width, height, 0, 0, width, height, bg, star_level,
                        noise, hfd_min, &mut mask, &mut stars, scan_budget,
                    );
                }
            }
            3 => {
                if star_level2 > 30.0 * noise {
                    scan_region(
                        lum, width, height, 0, 0, width, height, bg,
                        star_level2, noise, hfd_min, &mut mask, &mut stars,
                        scan_budget,
                    );
                }
            }
            2 => {
                scan_region(
                    lum, width, height, 0, 0, width, height, bg,
                    30.0 * noise, noise, hfd_min, &mut mask, &mut stars,
                    scan_budget,
                );
            }
            _ => {
                // retry 1: per-tile adaptive (12 columns on the long axis).
                let raster = 12usize;
                let (nx, ny) = if width >= height {
                    (raster, (raster * height / width).max(1))
                } else {
                    ((raster * width / height).max(1), raster)
                };
                let tw = width.div_ceil(nx);
                let th = height.div_ceil(ny);
                for ty in 0..ny {
                    for tx in 0..nx {
                        let x0 = tx * tw;
                        let y0 = ty * th;
                        let x1 = (x0 + tw).min(width);
                        let y1 = (y0 + th).min(height);
                        if x1 <= x0 || y1 <= y0 {
                            continue;
                        }
                        // Local bg/noise from a strided tile sample.
                        let mut samp: Vec<f32> = Vec::new();
                        let sstep = (((x1 - x0).max(y1 - y0)) / 64).max(1);
                        let mut yy = y0;
                        while yy < y1 {
                            let mut xx = x0;
                            while xx < x1 {
                                samp.push(lum[yy * width + xx]);
                                xx += sstep;
                            }
                            yy += sstep;
                        }
                        let (lbg, lnoise) = local_bg_noise(&samp);
                        scan_region(
                            lum, width, height, x0, y0, x1, y1, lbg,
                            7.0 * lnoise, lnoise.max(noise), hfd_min, &mut mask,
                            &mut stars, scan_budget,
                        );
                    }
                }
            }
        }
        retries -= 1;
    }

    // Brightest-first; trim to max_stars.
    stars.sort_by(|a, b| {
        b.flux.partial_cmp(&a.flux).unwrap_or(std::cmp::Ordering::Equal)
    });
    stars.truncate(max_stars);

    stars
        .into_iter()
        .map(|s| {
            (
                DetectedStar {
                    x: s.x,
                    y: s.y,
                    peak: s.peak,
                    flux: s.flux,
                    area: (s.hfd * s.hfd).max(1.0) as usize,
                    theta: 0.0,
                    eccentricity: 0.0,
                },
                s.snr,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn add_star(data: &mut [f32], width: usize, x: f32, y: f32, amp: f32, sigma: f32) {
        let height = data.len() / width;
        let r = (4.0 * sigma).ceil() as i32;
        let inv = 1.0 / (2.0 * sigma * sigma);
        for dy in -r..=r {
            for dx in -r..=r {
                let px = x as i32 + dx;
                let py = y as i32 + dy;
                if px >= 1 && (px as usize) < width - 1 && py >= 1 && (py as usize) < height - 1 {
                    let ddx = px as f32 - x;
                    let ddy = py as f32 - y;
                    data[py as usize * width + px as usize] +=
                        amp * (-(ddx * ddx + ddy * ddy) * inv).exp();
                }
            }
        }
    }

    /// The max_stars cap must select the BRIGHTEST stars frame-wide, not the
    /// first ones in scan order. Regression: on dense frames the row-major
    /// scan filled the cap from the top of the image and returned — the fast
    /// detector only ever saw a horizontal stripe (cocoon.fits: all 500
    /// stars at y < 1470 of 4176, 3% overlap with the slow path's top-100).
    #[test]
    fn cap_keeps_brightest_not_first_in_scan_order() {
        // Dense-frame regression (cocoon.fits): one bright wide halo eats the
        // brightest-pixel budget the ladder levels derive from, the working
        // level drops below hundreds of ordinary stars, and the row-major
        // scan fills the max_stars cap from the TOP of the frame and returns
        // — the fast detector only ever saw a horizontal stripe (all 500
        // cocoon stars at y < 1470 of 4176; 3% overlap with the slow path's
        // top-100). The cap must keep the brightest stars frame-wide.
        let width = 512;
        let height = 512;
        let mut data = vec![100.0_f32; width * height];
        // Wide bright halo: keeps levels 1-2 high (they only ever see it).
        add_star(&mut data, width, 256.0, 40.0, 30000.0, 12.0);
        // 10×18 grid of EQUAL ordinary stars over the whole frame — far more
        // than max_stars, all admitted by the same ladder level.
        for gy in 0..18 {
            for gx in 0..10 {
                add_star(
                    &mut data, width,
                    30.0 + gx as f32 * 48.0,
                    85.0 + gy as f32 * 21.0,
                    6000.0,
                    1.6,
                );
            }
        }
        // A row of clearly BRIGHTER stars at the very bottom — scanned last,
        // below the halo-driven levels 1-2, above the ordinary grid.
        for gx in 0..8 {
            add_star(&mut data, width, 50.0 + gx as f32 * 56.0, 486.0, 12000.0, 1.6);
        }
        // Deterministic mild noise so levels are sane.
        let mut rng = 12345u64;
        for v in data.iter_mut() {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            *v += ((rng >> 40) as f32 / (1u64 << 24) as f32 - 0.5) * 6.0;
        }

        let max_stars = 24;
        let stars = detect_stars_adaptive(&data, width, height, max_stars, 1.0);
        assert!(stars.len() >= 16, "expected detections, got {}", stars.len());
        let bright_bottom = stars.iter().filter(|(s, _)| s.y > 470.0).count();
        assert!(
            bright_bottom >= 6,
            "scan-order bias: only {}/8 of the brightest (bottom-row) stars made the cap; got {} stars, y range [{:.0},{:.0}]",
            bright_bottom,
            stars.len(),
            stars.iter().map(|(s, _)| s.y).fold(f32::MAX, f32::min),
            stars.iter().map(|(s, _)| s.y).fold(f32::MIN, f32::max),
        );
    }
}

#[cfg(all(test, feature = "debug-pipeline"))]
mod diag {
    use super::*;

    #[test]
    #[ignore]
    fn cocoon_top_structure() {
        let path = "tests/cocoon.fits";
        if !std::path::Path::new(path).exists() { return; }
        let (meta, pixels) = crate::formats::read_image(std::path::Path::new(path)).unwrap();
        let (lum, w, h, _, _) = crate::analysis::prepare_luminance(&meta, &pixels, true);
        let stars = detect_stars_adaptive(&lum, w, h, 8000, 0.8);
        eprintln!("total {}", stars.len());
        eprintln!("top 25 by flux (x y flux area~hfd2 snr):");
        for (s, snr) in stars.iter().take(25) {
            eprintln!("  {:7.1} {:7.1} {:10.0} {:5} {:7.1}", s.x, s.y, s.flux, s.area, snr);
        }
        // hfd distribution among top 500 vs rest
        let hfd = |a: usize| (a as f32).sqrt();
        let mut top: Vec<f32> = stars.iter().take(500).map(|(s, _)| hfd(s.area)).collect();
        let mut rest: Vec<f32> = stars.iter().skip(500).map(|(s, _)| hfd(s.area)).collect();
        top.sort_by(|a, b| a.total_cmp(b)); rest.sort_by(|a, b| a.total_cmp(b));
        eprintln!("hfd p10/p50/p90 top500: {:.1}/{:.1}/{:.1}  rest: {:.1}/{:.1}/{:.1}",
            top[top.len()/10], top[top.len()/2], top[top.len()*9/10],
            rest[rest.len()/10], rest[rest.len()/2], rest[rest.len()*9/10]);
        // surface-brightness ranking: flux/hfd^2 = flux/area
        let mut by_sb: Vec<&(DetectedStar, f32)> = stars.iter().collect();
        by_sb.sort_by(|a, b| (b.0.flux / b.0.area as f32).total_cmp(&(a.0.flux / a.0.area as f32)));
        eprintln!("top 25 by flux/area (x y flux area snr):");
        for e in by_sb.iter().take(25) {
            eprintln!("  {:7.1} {:7.1} {:10.0} {:5} {:7.1}", e.0.x, e.0.y, e.0.flux, e.0.area, e.1);
        }
    }
}
