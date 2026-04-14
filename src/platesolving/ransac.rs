use crate::platesolving::types::{ImageStar, ProjectedStar, RansacConfig, StarMatch};

/// RANSAC outlier rejection for Phase 2 individual star matching.
///
/// Given matched star pairs (image pixel ↔ catalog tangent plane),
/// finds the largest consistent set using random sampling + affine fitting.
pub struct RansacFilter;

impl RansacFilter {
    pub fn filter(
        matches: &[StarMatch],
        image_stars: &[ImageStar],
        catalog_stars: &[ProjectedStar],
        config: &RansacConfig,
    ) -> Vec<StarMatch> {
        if matches.len() < config.min_inliers {
            return Vec::new();
        }

        let n = matches.len();
        let mut best_inliers: Vec<usize> = Vec::new();
        let mut rng_state: u64 = 42;

        for _ in 0..config.max_iterations {
            let sample = random_sample(n, 3, &mut rng_state);

            // Fit affine from 3 sample pairs
            let Some(transform) = fit_affine_3pt(
                &sample,
                matches,
                image_stars,
                catalog_stars,
            ) else {
                continue;
            };

            let inliers = find_inliers(
                &transform,
                matches,
                image_stars,
                catalog_stars,
                config.threshold_px,
            );

            if inliers.len() > best_inliers.len() {
                best_inliers = inliers;
            }

            if best_inliers.len() as f64 > n as f64 * 0.8 {
                break;
            }
        }

        if best_inliers.len() < config.min_inliers {
            return Vec::new();
        }

        best_inliers
            .iter()
            .map(|&idx| matches[idx].clone())
            .collect()
    }
}

fn fit_affine_3pt(
    sample: &[usize],
    matches: &[StarMatch],
    image_stars: &[ImageStar],
    catalog_stars: &[ProjectedStar],
) -> Option<[f64; 6]> {
    if sample.len() < 3 {
        return None;
    }

    // Collect 3 point pairs
    let mut ix = [0.0; 3];
    let mut iy = [0.0; 3];
    let mut cx = [0.0; 3];
    let mut cy = [0.0; 3];

    for (k, &idx) in sample.iter().take(3).enumerate() {
        let m = &matches[idx];
        let is = &image_stars[m.image_idx];
        let cs = &catalog_stars[m.catalog_idx];
        ix[k] = is.x;
        iy[k] = is.y;
        cx[k] = cs.xi;
        cy[k] = cs.eta;
    }

    // Solve: [a, b, c] such that cx = a*ix + b*iy + c (and same for y)
    let det = ix[0] * (iy[1] - iy[2]) - iy[0] * (ix[1] - ix[2]) + (ix[1] * iy[2] - ix[2] * iy[1]);
    if det.abs() < 1e-20 {
        return None;
    }
    let inv = 1.0 / det;

    let a1 = (cx[0] * (iy[1] - iy[2]) + cx[1] * (iy[2] - iy[0]) + cx[2] * (iy[0] - iy[1])) * inv;
    let b1 = (cx[0] * (ix[2] - ix[1]) + cx[1] * (ix[0] - ix[2]) + cx[2] * (ix[1] - ix[0])) * inv;
    let c1 = (cx[0] * (ix[1] * iy[2] - ix[2] * iy[1]) + cx[1] * (ix[2] * iy[0] - ix[0] * iy[2]) + cx[2] * (ix[0] * iy[1] - ix[1] * iy[0])) * inv;

    let a2 = (cy[0] * (iy[1] - iy[2]) + cy[1] * (iy[2] - iy[0]) + cy[2] * (iy[0] - iy[1])) * inv;
    let b2 = (cy[0] * (ix[2] - ix[1]) + cy[1] * (ix[0] - ix[2]) + cy[2] * (ix[1] - ix[0])) * inv;
    let c2 = (cy[0] * (ix[1] * iy[2] - ix[2] * iy[1]) + cy[1] * (ix[2] * iy[0] - ix[0] * iy[2]) + cy[2] * (ix[0] * iy[1] - ix[1] * iy[0])) * inv;

    Some([a1, b1, c1, a2, b2, c2])
}

fn find_inliers(
    transform: &[f64; 6],
    matches: &[StarMatch],
    image_stars: &[ImageStar],
    catalog_stars: &[ProjectedStar],
    threshold_px: f64,
) -> Vec<usize> {
    let [a1, b1, c1, a2, b2, c2] = *transform;

    // Estimate scale to convert tangent-plane residuals to pixel-like units
    let scale = (a1 * a1 + a2 * a2).sqrt();
    let threshold = if scale > 1e-15 {
        threshold_px * scale
    } else {
        threshold_px * 1e-5
    };

    matches
        .iter()
        .enumerate()
        .filter_map(|(idx, m)| {
            let is = &image_stars[m.image_idx];
            let cs = &catalog_stars[m.catalog_idx];
            let pred_xi = a1 * is.x + b1 * is.y + c1;
            let pred_eta = a2 * is.x + b2 * is.y + c2;
            let dist = ((pred_xi - cs.xi).powi(2) + (pred_eta - cs.eta).powi(2)).sqrt();
            if dist < threshold {
                Some(idx)
            } else {
                None
            }
        })
        .collect()
}

fn xorshift64(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    x
}

fn random_sample(n: usize, k: usize, rng: &mut u64) -> Vec<usize> {
    let mut indices = Vec::with_capacity(k);
    while indices.len() < k {
        let idx = (xorshift64(rng) as usize) % n;
        if !indices.contains(&idx) {
            indices.push(idx);
        }
    }
    indices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ransac_filters_outliers() {
        let scale = 1e-4;
        let mut image_stars = Vec::new();
        let mut catalog_stars = Vec::new();
        let mut matches = Vec::new();

        for i in 0..15 {
            let x = 100.0 + (i as f64) * 50.0;
            let y = 200.0 + (i as f64) * 30.0;
            image_stars.push(ImageStar { x, y, flux: 1000.0 });
            catalog_stars.push(ProjectedStar {
                xi: x * scale,
                eta: y * scale,
                mag: 8.0,
                ra: 180.0,
                dec: 45.0,
            });
            matches.push(StarMatch {
                image_idx: i,
                catalog_idx: i,
                residual_px: 0.0,
            });
        }

        // Add 5 outliers
        for i in 0..5 {
            let idx = 15 + i;
            image_stars.push(ImageStar { x: 500.0, y: 500.0, flux: 500.0 });
            catalog_stars.push(ProjectedStar {
                xi: 0.1, eta: 0.1, mag: 10.0, ra: 180.0, dec: 45.0,
            });
            matches.push(StarMatch {
                image_idx: idx,
                catalog_idx: idx,
                residual_px: 0.0,
            });
        }

        let config = RansacConfig::default();
        let filtered = RansacFilter::filter(&matches, &image_stars, &catalog_stars, &config);
        assert!(
            filtered.len() >= 10,
            "Should keep most inliers, got {}",
            filtered.len()
        );
    }
}
