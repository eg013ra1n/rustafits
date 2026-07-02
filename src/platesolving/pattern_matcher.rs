/// Quad-based plate solving: nearest-neighbour quad matching with center correspondences.
///
/// Algorithm:
/// 1. For each star, find 3 nearest neighbors → 1 quad per star
/// 2. Compute 6 pairwise distances, sort, normalize by longest → scale-invariant ratios
/// 3. Brute-force compare ratios with tolerance (early exit on first mismatch)
/// 4. Each match gives a center correspondence (no star ordering ambiguity)
/// 5. Fit affine from 3+ center correspondences

use std::collections::HashSet;

use nalgebra::{DMatrix, DVector};

pub const QUAD_SIZE: usize = 4;
pub const NUM_EDGES: usize = 6; // C(4,2)

/// A quad: one star + its 3 nearest neighbors.
#[derive(Clone, Debug)]
pub struct Quad {
    pub star_indices: [usize; QUAD_SIZE],
    pub ratios: [f64; NUM_EDGES - 1], // first 5 ratios (last is always 1.0)
    pub center: (f64, f64),
    pub longest_dist: f64,
}

/// A matched quad pair with center correspondence.
#[derive(Clone, Debug)]
pub struct QuadMatch {
    pub image_center: (f64, f64),
    pub catalog_center: (f64, f64),
    pub scale_ratio: f64,
}

/// 6-parameter affine transform.
#[derive(Clone, Debug)]
pub struct AffineTransform {
    pub a1: f64,
    pub b1: f64,
    pub c1: f64,
    pub a2: f64,
    pub b2: f64,
    pub c2: f64,
}

impl AffineTransform {
    pub fn pixel_scale_x(&self) -> f64 {
        (self.a1 * self.a1 + self.b1 * self.b1).sqrt()
    }

    pub fn pixel_scale_y(&self) -> f64 {
        (self.a2 * self.a2 + self.b2 * self.b2).sqrt()
    }

    pub fn rotation_rad(&self) -> f64 {
        (-self.b1).atan2(self.a1)
    }

    pub fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        (
            self.a1 * x + self.b1 * y + self.c1,
            self.a2 * x + self.b2 * y + self.c2,
        )
    }
}

fn distance(a: (f64, f64), b: (f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
}

/// Build quads using nearest-neighbor selection (1 quad per star: the star
/// plus its 3 nearest neighbors). Thin wrapper over [`build_quads_multi`]
/// with `group_size = QUAD_SIZE`; behaviour is byte-for-byte identical to the
/// original implementation so dense-field solving is unaffected.
pub fn build_quads(stars: &[(f64, f64)], max_stars: usize) -> Vec<Quad> {
    build_quads_multi(stars, max_stars, QUAD_SIZE)
}

/// Append the 4-edge ratios + center for the 4 stars `idx` to `quads`,
/// deduplicating on the sorted star-index set and (legacy) on the center.
#[allow(clippy::too_many_arguments)]
fn push_quad_from_indices(
    stars: &[(f64, f64)],
    idx: [usize; QUAD_SIZE],
    quads: &mut Vec<Quad>,
    seen: &mut HashSet<[usize; QUAD_SIZE]>,
    dedup_epsilon: f64,
) {
    let mut key = idx;
    key.sort_unstable();
    if !seen.insert(key) {
        return;
    }

    let cx =
        (stars[idx[0]].0 + stars[idx[1]].0 + stars[idx[2]].0 + stars[idx[3]].0) / 4.0;
    let cy =
        (stars[idx[0]].1 + stars[idx[1]].1 + stars[idx[2]].1 + stars[idx[3]].1) / 4.0;

    // Legacy center dedup (preserves original build_quads output exactly).
    if quads.iter().any(|q: &Quad| {
        (q.center.0 - cx).abs() < dedup_epsilon && (q.center.1 - cy).abs() < dedup_epsilon
    }) {
        return;
    }

    let mut dists = [0.0f64; NUM_EDGES];
    let mut k = 0;
    for a in 0..QUAD_SIZE {
        for b in (a + 1)..QUAD_SIZE {
            dists[k] = distance(stars[idx[a]], stars[idx[b]]);
            k += 1;
        }
    }
    dists.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let longest = dists[NUM_EDGES - 1];
    if longest < 1e-15 {
        return;
    }
    let ratios = [
        dists[0] / longest,
        dists[1] / longest,
        dists[2] / longest,
        dists[3] / longest,
        dists[4] / longest,
    ];

    quads.push(Quad {
        star_indices: idx,
        ratios,
        center: (cx, cy),
        longest_dist: longest,
    });
}

/// Quad construction with a tunable per-star group size.
///
/// For each of the `max_stars` brightest stars, take its `group_size - 1`
/// nearest neighbors to form a local group of `group_size` stars, then emit
/// every C(group_size, 4) sub-quad. `group_size == 4` reproduces the classic
/// "1 quad per star" behaviour exactly; larger groups (6 → 15 quads/star,
/// 5 → 5 quads/star) massively densify the quad pool for sparse,
/// long-focal-length fields where only a handful of catalog-depth stars are
/// present — without which a correct quad almost never forms. Quads are
/// deduplicated globally by their sorted star-index set (and, for
/// `group_size == 4`, additionally by center to stay identical to the
/// original).
pub fn build_quads_multi(
    stars: &[(f64, f64)],
    max_stars: usize,
    group_size: usize,
) -> Vec<Quad> {
    let n = stars.len().min(max_stars);
    let g = group_size.max(QUAD_SIZE);
    if n < QUAD_SIZE {
        return Vec::new();
    }

    let mut quads = Vec::with_capacity(n);
    let mut seen: HashSet<[usize; QUAD_SIZE]> = HashSet::new();
    let dedup_epsilon = 1e-6;
    let n_neighbors = g - 1;

    // Scratch buffer of (index, dist) reused per star.
    let mut nbrs: Vec<(usize, f64)> = Vec::with_capacity(n.saturating_sub(1));

    for i in 0..n {
        nbrs.clear();
        for j in 0..n {
            if i == j {
                continue;
            }
            nbrs.push((j, distance(stars[i], stars[j])));
        }
        if nbrs.len() < n_neighbors {
            continue;
        }
        // Partial selection: the `n_neighbors` smallest distances, then sort
        // those so group membership is deterministic. For group_size == 4
        // this yields exactly the original [i, n0, n1, n2] (the 3 nearest).
        nbrs.select_nth_unstable_by(n_neighbors - 1, |a, b| {
            a.1.partial_cmp(&b.1).unwrap()
        });
        nbrs.truncate(n_neighbors);
        nbrs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // Local group: center star + its nearest neighbors.
        let mut group = [0usize; 6]; // g ≤ 6 in practice; clamp below
        let gsz = g.min(group.len());
        group[0] = i;
        for (slot, &(idx, _)) in nbrs.iter().take(gsz - 1).enumerate() {
            group[slot + 1] = idx;
        }

        // Emit every 4-subset of the group.
        for a in 0..gsz {
            for b in (a + 1)..gsz {
                for c in (b + 1)..gsz {
                    for d in (c + 1)..gsz {
                        push_quad_from_indices(
                            stars,
                            [group[a], group[b], group[c], group[d]],
                            &mut quads,
                            &mut seen,
                            dedup_epsilon,
                        );
                    }
                }
            }
        }
    }

    quads
}

/// Brute-force match image quads vs catalog quads.
///
/// For each ratio pair: `min(a,b)/max(a,b) > (1.0 - tolerance)`.
/// Early exit on first mismatch — most pairs fail on the first ratio.
pub fn match_quads(
    image_quads: &[Quad],
    catalog_quads: &[Quad],
    tolerance: f64,
) -> Vec<QuadMatch> {
    let threshold = 1.0 - tolerance;
    let mut matches = Vec::new();

    for iq in image_quads {
        for cq in catalog_quads {
            let mut similar = true;
            for k in 0..5 {
                let a = iq.ratios[k];
                let b = cq.ratios[k];
                let ratio = if a > b { b / a } else { a / b };
                if ratio < threshold {
                    similar = false;
                    break;
                }
            }
            if similar {
                let scale_ratio = if iq.longest_dist > 1e-15 {
                    cq.longest_dist / iq.longest_dist
                } else {
                    0.0
                };
                matches.push(QuadMatch {
                    image_center: iq.center,
                    catalog_center: cq.center,
                    scale_ratio,
                });
            }
        }
    }

    // Outlier removal: filter by scale_ratio consistency
    if matches.len() >= 3 {
        let mut scale_ratios: Vec<f64> = matches.iter().map(|m| m.scale_ratio).collect();
        scale_ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = scale_ratios[scale_ratios.len() / 2];

        if median > 1e-15 {
            matches.retain(|m| {
                let dev = (m.scale_ratio - median).abs() / median;
                dev < 0.2 // 20% deviation from median
            });
        }
    }

    matches
}

/// Fit 6-parameter affine from 3+ matched quad centers.
///
/// Xref = a1*x + b1*y + c1
/// Yref = a2*x + b2*y + c2
///
/// Sanity check: (a1²+b1²)/(a2²+b2²) must be in [0.9, 1.1].
pub fn fit_affine_from_centers(matches: &[QuadMatch]) -> Option<AffineTransform> {
    let n = matches.len();
    if n < 3 {
        return None;
    }

    // Build design matrix A (n×3) and targets (n-vectors) for the pair
    // [x, y, 1] · [a, b, c]ᵀ = target. Solve via SVD pseudo-inverse instead
    // of normal-equation LU: the condition number of AᵀA is the SQUARE of
    // A's, so LU on AᵀA is numerically unstable whenever points are
    // near-collinear. SVD degrades gracefully and returns the minimum-norm
    // solution even in rank-deficient cases.
    let mut a_mat = DMatrix::zeros(n, 3);
    let mut b_vec_x = DVector::zeros(n);
    let mut b_vec_y = DVector::zeros(n);

    for (i, m) in matches.iter().enumerate() {
        a_mat[(i, 0)] = m.image_center.0;
        a_mat[(i, 1)] = m.image_center.1;
        a_mat[(i, 2)] = 1.0;
        b_vec_x[i] = m.catalog_center.0;
        b_vec_y[i] = m.catalog_center.1;
    }

    let svd = a_mat.clone().svd(true, true);
    let eps = 1e-12;
    let params_x = svd.solve(&b_vec_x, eps).ok()?;
    let params_y = svd.solve(&b_vec_y, eps).ok()?;

    let transform = AffineTransform {
        a1: params_x[0],
        b1: params_x[1],
        c1: params_x[2],
        a2: params_y[0],
        b2: params_y[1],
        c2: params_y[2],
    };

    // Sanity: X/Y axis scales should match to within ~5%. Real pinhole
    // optics have essentially no aspect skew; the old 0.8–1.2 band let
    // through badly-conditioned fits.
    let scale_x_sq = transform.a1 * transform.a1 + transform.b1 * transform.b1;
    let scale_y_sq = transform.a2 * transform.a2 + transform.b2 * transform.b2;
    if scale_y_sq < 1e-30 {
        return None;
    }
    let ratio = scale_x_sq / scale_y_sq;
    if !(0.9..=1.1).contains(&ratio) {
        return None;
    }

    Some(transform)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_quads_basic() {
        // 5 stars in a cross pattern
        let stars = vec![
            (100.0, 100.0),
            (200.0, 100.0),
            (150.0, 50.0),
            (150.0, 150.0),
            (300.0, 200.0),
        ];
        let quads = build_quads(&stars, 5);
        assert!(!quads.is_empty(), "Should produce at least 1 quad");
        for q in &quads {
            assert!(q.longest_dist > 0.0);
            for r in &q.ratios {
                assert!(*r >= 0.0 && *r <= 1.0, "Ratio should be in [0,1]: {}", r);
            }
        }
    }

    #[test]
    fn matching_identical_patterns() {
        // Same geometric pattern at different scales and positions
        let image_stars: Vec<(f64, f64)> = vec![
            (100.0, 100.0),
            (200.0, 100.0),
            (150.0, 50.0),
            (150.0, 150.0),
        ];

        // Same pattern scaled by 0.001 and shifted
        let scale = 0.001;
        let catalog_stars: Vec<(f64, f64)> = image_stars
            .iter()
            .map(|(x, y)| (x * scale + 5.0, y * scale + 3.0))
            .collect();

        let img_quads = build_quads(&image_stars, 4);
        let cat_quads = build_quads(&catalog_stars, 4);

        assert!(!img_quads.is_empty());
        assert!(!cat_quads.is_empty());

        let matches = match_quads(&img_quads, &cat_quads, 0.01);
        assert!(
            !matches.is_empty(),
            "Identical patterns at different scales should match"
        );
    }

    #[test]
    fn matching_different_patterns_no_match() {
        let image_stars: Vec<(f64, f64)> = vec![
            (100.0, 100.0),
            (200.0, 100.0),
            (150.0, 50.0),
            (150.0, 150.0),
        ];

        // Completely different pattern
        let catalog_stars: Vec<(f64, f64)> = vec![
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 5.0),
            (3.0, 3.0),
        ];

        let img_quads = build_quads(&image_stars, 4);
        let cat_quads = build_quads(&catalog_stars, 4);
        let matches = match_quads(&img_quads, &cat_quads, 0.01);

        // Might get 0 or very few false matches
        assert!(
            matches.len() <= 1,
            "Different patterns should get few/no matches, got {}",
            matches.len()
        );
    }

    #[test]
    fn affine_fit_from_centers() {
        // Known transform: xref = 0.001*x + 5.0, yref = 0.001*y + 3.0
        let matches = vec![
            QuadMatch {
                image_center: (100.0, 100.0),
                catalog_center: (5.1, 3.1),
                scale_ratio: 0.001,
            },
            QuadMatch {
                image_center: (200.0, 100.0),
                catalog_center: (5.2, 3.1),
                scale_ratio: 0.001,
            },
            QuadMatch {
                image_center: (150.0, 200.0),
                catalog_center: (5.15, 3.2),
                scale_ratio: 0.001,
            },
            QuadMatch {
                image_center: (300.0, 50.0),
                catalog_center: (5.3, 3.05),
                scale_ratio: 0.001,
            },
        ];

        let transform = fit_affine_from_centers(&matches).expect("Fit should succeed");

        // Check: applying transform to image center should give catalog center
        let (xr, yr) = transform.apply(250.0, 150.0);
        let expected_x = 0.001 * 250.0 + 5.0;
        let expected_y = 0.001 * 150.0 + 3.0;
        assert!(
            (xr - expected_x).abs() < 0.01,
            "X: expected {expected_x}, got {xr}"
        );
        assert!(
            (yr - expected_y).abs() < 0.01,
            "Y: expected {expected_y}, got {yr}"
        );
    }
}
