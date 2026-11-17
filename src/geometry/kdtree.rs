//! A static 2-D KD-tree for nearest-neighbour lookups within a radius —
//! replaces the O(N·M) scan in frame-to-frame correspondence. Built once
//! per reference star list (≤ a few thousand points); queries are O(log n).

use rayon::prelude::*;

/// A range this long builds its two halves on two threads.
const PARALLEL_BUILD_MIN: usize = 2048;

/// A range this long is sorted on the pool (the levels above it, where one sort is most of the build's critical path).
const PARALLEL_SORT_MIN: usize = 8192;

/// One point while the tree is built.
struct Item {
    x: f64,
    y: f64,
    idx: usize,
    /// Its position in the range being sorted, set before each level's sort (the stable sort's tie-break).
    pos: usize,
}

#[derive(Debug)]
pub struct KdTree2 {
    /// Points in tree order; `idx[i]` is the original index of `pts[i]`.
    pts: Vec<(f64, f64)>,
    idx: Vec<usize>,
}

impl KdTree2 {
    pub fn build(points: &[(f64, f64)]) -> KdTree2 {
        // NaN keys make `partial_cmp` an inconsistent order, so the layout would depend on the sort's internals: that
        // input keeps the original recursion. Every caller hands finite points.
        if points.iter().any(|p| p.0.is_nan() || p.1.is_nan()) {
            let mut idx: Vec<usize> = (0..points.len()).collect();
            let mut pts: Vec<(f64, f64)> = points.to_vec();
            Self::build_rec(&mut pts, &mut idx, 0);
            return KdTree2 { pts, idx };
        }
        let mut items: Vec<Item> = points
            .iter()
            .enumerate()
            .map(|(i, &(x, y))| Item { x, y, idx: i, pos: 0 })
            .collect();
        Self::build_items(&mut items, 0);
        KdTree2 {
            pts: items.iter().map(|it| (it.x, it.y)).collect(),
            idx: items.iter().map(|it| it.idx).collect(),
        }
    }

    /// [`Self::build_rec`]'s tree, arranged bit for bit alike, on one array of records: a stable sort by the axis key
    /// is an unstable sort by (key, the record's position in the range as the level found it), so the sort moves one
    /// 32-byte record per swap and allocates nothing; the two halves are independent, so a large range builds them on
    /// two threads (the layout does not depend on which). Keys must not be NaN.
    fn build_items(items: &mut [Item], depth: usize) {
        let n = items.len();
        if n <= 1 {
            return;
        }
        let axis = depth % 2;
        let mid = n / 2;
        for (p, it) in items.iter_mut().enumerate() {
            it.pos = p;
        }
        // The order is total (no NaN, the position breaks every tie), so the sorted sequence is unique and a parallel sort
        // of a long range gives the sequential one's.
        let by_x = |a: &Item, b: &Item| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal).then(a.pos.cmp(&b.pos));
        let by_y = |a: &Item, b: &Item| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal).then(a.pos.cmp(&b.pos));
        match (axis == 0, n >= PARALLEL_SORT_MIN) {
            (true, true) => items.par_sort_unstable_by(by_x),
            (true, false) => items.sort_unstable_by(by_x),
            (false, true) => items.par_sort_unstable_by(by_y),
            (false, false) => items.sort_unstable_by(by_y),
        }
        let (left, right) = items.split_at_mut(mid);
        let right = &mut right[1..];
        if n >= PARALLEL_BUILD_MIN {
            rayon::join(|| Self::build_items(left, depth + 1), || Self::build_items(right, depth + 1));
        } else {
            Self::build_items(left, depth + 1);
            Self::build_items(right, depth + 1);
        }
    }

    pub fn len(&self) -> usize {
        self.pts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }

    /// Recursively arrange `[lo, hi)` so the median on the split axis sits
    /// in the middle, alternating axes by depth. Implicit balanced tree:
    /// node = middle index of its range.
    fn build_rec(pts: &mut [(f64, f64)], idx: &mut [usize], depth: usize) {
        let n = pts.len();
        if n <= 1 {
            return;
        }
        let axis = depth % 2;
        let mid = n / 2;
        // Selection by axis: sort the pair arrays together (n log n per
        // level is fine at these sizes and keeps the code obvious).
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            let ka = if axis == 0 { pts[a].0 } else { pts[a].1 };
            let kb = if axis == 0 { pts[b].0 } else { pts[b].1 };
            ka.partial_cmp(&kb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let sorted_pts: Vec<_> = order.iter().map(|&o| pts[o]).collect();
        let sorted_idx: Vec<_> = order.iter().map(|&o| idx[o]).collect();
        pts.copy_from_slice(&sorted_pts);
        idx.copy_from_slice(&sorted_idx);
        let (left_p, right_p) = pts.split_at_mut(mid);
        let (left_i, right_i) = idx.split_at_mut(mid);
        Self::build_rec(left_p, left_i, depth + 1);
        Self::build_rec(&mut right_p[1..], &mut right_i[1..], depth + 1);
    }

    /// Nearest point within `radius` (inclusive), as `(original index,
    /// distance)`.
    pub fn nearest_within(&self, x: f64, y: f64, radius: f64) -> Option<(usize, f64)> {
        if self.pts.is_empty() {
            return None;
        }
        let mut best: Option<(usize, f64)> = None;
        let mut best_d2 = radius * radius;
        self.search(0, self.pts.len(), 0, x, y, &mut best, &mut best_d2);
        best.map(|(i, d2)| (i, d2.sqrt()))
    }

    /// Every point within `radius` (inclusive) of `(x, y)`, as `(original index, distance)`, nearest first (ties by
    /// original index).
    pub fn within(&self, x: f64, y: f64, radius: f64) -> Vec<(usize, f64)> {
        let mut found: Vec<(usize, f64)> = Vec::new();
        self.collect(0, self.pts.len(), 0, x, y, radius * radius, &mut found);
        for f in &mut found {
            f.1 = f.1.sqrt();
        }
        found.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        found
    }

    #[allow(clippy::too_many_arguments)]
    fn collect(
        &self,
        lo: usize,
        hi: usize,
        depth: usize,
        x: f64,
        y: f64,
        r2: f64,
        out: &mut Vec<(usize, f64)>,
    ) {
        if lo >= hi {
            return;
        }
        let mid = lo + (hi - lo) / 2;
        let (px, py) = self.pts[mid];
        let d2 = (px - x).powi(2) + (py - y).powi(2);
        if d2 <= r2 {
            out.push((self.idx[mid], d2));
        }
        let diff = if depth % 2 == 0 { x - px } else { y - py };
        let (near, far) = if diff < 0.0 {
            ((lo, mid), (mid + 1, hi))
        } else {
            ((mid + 1, hi), (lo, mid))
        };
        self.collect(near.0, near.1, depth + 1, x, y, r2, out);
        if diff * diff <= r2 {
            self.collect(far.0, far.1, depth + 1, x, y, r2, out);
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn search(
        &self,
        lo: usize,
        hi: usize,
        depth: usize,
        x: f64,
        y: f64,
        best: &mut Option<(usize, f64)>,
        best_d2: &mut f64,
    ) {
        if lo >= hi {
            return;
        }
        let mid = lo + (hi - lo) / 2;
        let (px, py) = self.pts[mid];
        let d2 = (px - x).powi(2) + (py - y).powi(2);
        if d2 <= *best_d2 && best.map(|b| d2 < b.1).unwrap_or(true) {
            *best = Some((self.idx[mid], d2));
            *best_d2 = d2;
        }
        let axis = depth % 2;
        let diff = if axis == 0 { x - px } else { y - py };
        let (near, far) = if diff < 0.0 {
            ((lo, mid), (mid + 1, hi))
        } else {
            ((mid + 1, hi), (lo, mid))
        };
        self.search(near.0, near.1, depth + 1, x, y, best, best_d2);
        if diff * diff <= *best_d2 {
            self.search(far.0, far.1, depth + 1, x, y, best, best_d2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn points(n: usize, seed: u64) -> Vec<(f64, f64)> {
        let mut s = seed;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s % 1_000_000) as f64 / 1_000_000.0
        };
        (0..n).map(|_| (next() * 6000.0, next() * 4000.0)).collect()
    }

    fn brute(pts: &[(f64, f64)], x: f64, y: f64, r: f64) -> Option<(usize, f64)> {
        let mut best: Option<(usize, f64)> = None;
        for (i, (px, py)) in pts.iter().enumerate() {
            let d = ((px - x).powi(2) + (py - y).powi(2)).sqrt();
            if d <= r && best.map(|b| d < b.1).unwrap_or(true) {
                best = Some((i, d));
            }
        }
        best
    }

    #[test]
    fn within_matches_brute_force() {
        let pts = points(300, 7);
        let tree = KdTree2::build(&pts);
        let got = tree.within(3000.0, 2000.0, 900.0);
        let mut want: Vec<(usize, f64)> = pts
            .iter()
            .enumerate()
            .map(|(i, &(px, py))| (i, ((px - 3000.0).powi(2) + (py - 2000.0).powi(2)).sqrt()))
            .filter(|&(_, d)| d <= 900.0)
            .collect();
        want.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        assert!(!want.is_empty());
        assert_eq!(got.len(), want.len());
        for (g, w) in got.iter().zip(&want) {
            assert_eq!(g.0, w.0);
            assert!((g.1 - w.1).abs() < 1e-9);
        }
    }

    #[test]
    fn nearest_within_matches_brute_force() {
        let pts = points(2000, 77);
        let tree = KdTree2::build(&pts);
        assert_eq!(tree.len(), 2000);
        for (qx, qy) in points(500, 99) {
            let a = tree.nearest_within(qx, qy, 40.0);
            let b = brute(&pts, qx, qy, 40.0);
            match (a, b) {
                (None, None) => {}
                (Some((ia, da)), Some((ib, db))) => {
                    assert!((da - db).abs() < 1e-9, "distance mismatch");
                    // Ties at the same distance are legal either way.
                    if ia != ib {
                        assert!((da - db).abs() < 1e-9);
                    }
                }
                other => panic!("mismatch at ({qx},{qy}): {other:?}"),
            }
        }
    }

    /// The faster build lays the tree out exactly as the original recursion (stable sorts by the level's axis) does,
    /// ties on a coordinate included (a lattice of centroids has many).
    #[test]
    fn build_layout_equals_the_original_recursion() {
        let mut cases: Vec<Vec<(f64, f64)>> = vec![points(1, 3), points(2, 3), points(7, 5), points(300, 7), points(5000, 11)];
        // Ties: half-integer lattice, repeated points, signed zeros.
        let lattice: Vec<(f64, f64)> = (0..3000).map(|i| (((i * 7) % 40) as f64 * 0.5, ((i * 13) % 25) as f64 * 0.5)).collect();
        cases.push(lattice);
        cases.push(vec![(0.0, 0.0), (-0.0, 0.0), (0.0, -0.0), (1.0, 1.0), (1.0, 1.0), (0.0, 1.0)]);
        for pts in cases {
            let fast = KdTree2::build(&pts);
            let mut idx: Vec<usize> = (0..pts.len()).collect();
            let mut slow_pts = pts.clone();
            KdTree2::build_rec(&mut slow_pts, &mut idx, 0);
            assert_eq!(fast.idx, idx, "{} points", pts.len());
            assert_eq!(
                fast.pts.iter().map(|p| (p.0.to_bits(), p.1.to_bits())).collect::<Vec<_>>(),
                slow_pts.iter().map(|p| (p.0.to_bits(), p.1.to_bits())).collect::<Vec<_>>(),
                "{} points",
                pts.len()
            );
        }
    }

    #[test]
    fn empty_tree_returns_none() {
        let tree = KdTree2::build(&[]);
        assert!(tree.nearest_within(1.0, 1.0, 10.0).is_none());
    }

    #[test]
    fn radius_is_inclusive_and_exact() {
        let tree = KdTree2::build(&[(0.0, 0.0), (10.0, 0.0)]);
        assert_eq!(tree.nearest_within(4.0, 0.0, 4.0).unwrap().0, 0);
        assert!(tree.nearest_within(5.0, 0.0, 4.9).is_none());
        assert_eq!(tree.nearest_within(7.0, 0.0, 4.0).unwrap().0, 1);
    }
}
