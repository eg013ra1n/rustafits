//! Green-site helpers for measuring a raw CFA mosaic without interpolation.
//!
//! Star measurement on a one-shot-colour frame can work on the green sites
//! alone: half the pixels, all of one colour, and no interpolated values. These
//! helpers address those sites for a given [`BayerPattern`], in file row order
//! (the order the reader returns the samples in — the same order the pattern
//! addresses).
//!
//! A 2×2 cell holds two green sites: **G1** on the cell's first row, **G2** on
//! its second (RGGB: G1 = (1, 0), G2 = (0, 1); GRBG: G1 = (0, 0), G2 = (1, 1)).
//! The half-resolution outputs use the super-pixel grid of
//! [`super::debayer::super_pixel_debayer_f32`]: `width / 2 × height / 2` cells,
//! an odd last row or column left out.
//!
//! Every helper needs a real pattern; [`BayerPattern::None`] gives an empty
//! result (and `0 × 0` dimensions).

use rayon::prelude::*;

use crate::types::BayerPattern;

/// `(x, y)` of G1 and G2 inside the pattern's 2×2 cell.
fn green_sites(pattern: BayerPattern) -> Option<[(usize, usize); 2]> {
    match pattern {
        BayerPattern::Rggb | BayerPattern::Bggr => Some([(1, 0), (0, 1)]),
        BayerPattern::Grbg | BayerPattern::Gbrg => Some([(0, 0), (1, 1)]),
        BayerPattern::None => None,
    }
}

/// `true` on the green CFA sites of a `width × height` mosaic, row-major.
/// Empty for [`BayerPattern::None`].
pub fn green_mask(width: usize, height: usize, pattern: BayerPattern) -> Vec<bool> {
    let Some([(g1x, _), _]) = green_sites(pattern) else {
        return Vec::new();
    };
    // The green sites form a checkerboard: green where (x + y) has the parity
    // of G1's column (G1 sits on row 0).
    let mut mask = vec![false; width * height];
    mask.par_chunks_mut(width.max(1)).enumerate().for_each(|(y, row)| {
        for (x, m) in row.iter_mut().enumerate() {
            *m = (x + y) & 1 == g1x;
        }
    });
    mask
}

/// `G1 + G2` of every 2×2 cell: a half-resolution green plane of real samples,
/// `(plane, width / 2, height / 2)`. Equal to twice the green plane of
/// [`super::debayer::super_pixel_debayer_f32`]. Empty for [`BayerPattern::None`].
pub fn green_sum(mosaic: &[f32], width: usize, height: usize, pattern: BayerPattern) -> (Vec<f32>, usize, usize) {
    let Some([(g1x, g1y), (g2x, g2y)]) = green_sites(pattern) else {
        return (Vec::new(), 0, 0);
    };
    assert_eq!(mosaic.len(), width * height, "mosaic must hold width × height samples");
    let (out_w, out_h) = (width / 2, height / 2);
    let mut out = vec![0f32; out_w * out_h];
    out.par_chunks_mut(out_w.max(1)).enumerate().for_each(|(y, row)| {
        let g1_row = &mosaic[(2 * y + g1y) * width..][..2 * out_w];
        let g2_row = &mosaic[(2 * y + g2y) * width..][..2 * out_w];
        for ((v, a), b) in row.iter_mut().zip(g1_row.chunks_exact(2)).zip(g2_row.chunks_exact(2)) {
            *v = a[g1x] + b[g2x];
        }
    });
    (out, out_w, out_h)
}

/// One green phase: the G1 sample of every 2×2 cell, `(plane, width / 2,
/// height / 2)` — a rectangular plane of real samples, none interpolated and
/// none red or blue, for sky and noise statistics. Empty for
/// [`BayerPattern::None`].
pub fn green_phase(mosaic: &[f32], width: usize, height: usize, pattern: BayerPattern) -> (Vec<f32>, usize, usize) {
    let Some([(g1x, g1y), _]) = green_sites(pattern) else {
        return (Vec::new(), 0, 0);
    };
    assert_eq!(mosaic.len(), width * height, "mosaic must hold width × height samples");
    let (out_w, out_h) = (width / 2, height / 2);
    let mut out = vec![0f32; out_w * out_h];
    out.par_chunks_mut(out_w.max(1)).enumerate().for_each(|(y, row)| {
        let g1_row = &mosaic[(2 * y + g1y) * width..][..2 * out_w];
        for (v, a) in row.iter_mut().zip(g1_row.chunks_exact(2)) {
            *v = a[g1x];
        }
    });
    (out, out_w, out_h)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::processing::debayer::super_pixel_debayer_f32;

    const PATTERNS: [(BayerPattern, &str); 4] = [
        (BayerPattern::Rggb, "RGGB"),
        (BayerPattern::Grbg, "GRBG"),
        (BayerPattern::Gbrg, "GBRG"),
        (BayerPattern::Bggr, "BGGR"),
    ];

    /// Even, odd and degenerate sizes.
    const SIZES: [(usize, usize); 8] = [(8, 6), (7, 5), (6, 7), (9, 4), (1, 1), (2, 2), (1, 6), (0, 0)];

    /// Whether `(x, y)` is green, read off the pattern's name.
    fn is_green(name: &str, x: usize, y: usize) -> bool {
        name.as_bytes()[(y & 1) * 2 + (x & 1)] == b'G'
    }

    /// Every sample distinct, so an output value names the site it came from.
    fn indexed_mosaic(w: usize, h: usize) -> Vec<f32> {
        (0..w * h).map(|i| i as f32).collect()
    }

    #[test]
    fn green_mask_marks_the_green_sites() {
        for (pattern, name) in PATTERNS {
            for (w, h) in SIZES {
                let mask = green_mask(w, h, pattern);
                assert_eq!(mask.len(), w * h);
                for (i, &m) in mask.iter().enumerate() {
                    assert_eq!(m, is_green(name, i % w, i / w), "{name} {w}x{h} at ({}, {})", i % w, i / w);
                }
            }
        }
    }

    #[test]
    fn green_sum_adds_the_two_greens_of_each_cell() {
        for (pattern, name) in PATTERNS {
            for (w, h) in SIZES {
                let mosaic = indexed_mosaic(w, h);
                let (sum, sw, sh) = green_sum(&mosaic, w, h, pattern);
                assert_eq!((sw, sh, sum.len()), (w / 2, h / 2, (w / 2) * (h / 2)), "{name} {w}x{h}");
                for cy in 0..sh {
                    for cx in 0..sw {
                        let want: f32 = [(0, 0), (1, 0), (0, 1), (1, 1)]
                            .iter()
                            .map(|&(dx, dy)| (2 * cx + dx, 2 * cy + dy))
                            .filter(|&(x, y)| is_green(name, x, y))
                            .map(|(x, y)| mosaic[y * w + x])
                            .sum();
                        assert_eq!(sum[cy * sw + cx], want, "{name} {w}x{h} cell ({cx}, {cy})");
                    }
                }
            }
        }
    }

    #[test]
    fn green_sum_is_twice_the_super_pixel_green() {
        for (pattern, name) in PATTERNS {
            for (w, h) in [(64, 48), (65, 47), (130, 3)] {
                let mosaic: Vec<f32> = (0..w * h).map(|i| ((i * 7919) % 4093) as f32 + 0.25).collect();
                let (sum, sw, sh) = green_sum(&mosaic, w, h, pattern);
                let (rgb, rw, rh) = super_pixel_debayer_f32(&mosaic, w, h, pattern);
                assert_eq!((sw, sh), (rw, rh));
                let green = &rgb[rw * rh..2 * rw * rh];
                for (i, (&s, &g)) in sum.iter().zip(green).enumerate() {
                    assert_eq!(s, 2.0 * g, "{name} {w}x{h} cell {i}");
                }
            }
        }
    }

    #[test]
    fn green_phase_takes_the_first_row_green_of_each_cell() {
        for (pattern, name) in PATTERNS {
            for (w, h) in SIZES {
                let mosaic = indexed_mosaic(w, h);
                let (phase, pw, ph) = green_phase(&mosaic, w, h, pattern);
                assert_eq!((pw, ph, phase.len()), (w / 2, h / 2, (w / 2) * (h / 2)), "{name} {w}x{h}");
                for cy in 0..ph {
                    for cx in 0..pw {
                        let i = phase[cy * pw + cx] as usize;
                        let (x, y) = (i % w, i / w);
                        assert!(is_green(name, x, y), "{name} {w}x{h} cell ({cx}, {cy}) took ({x}, {y})");
                        assert_eq!((x / 2, y), (cx, 2 * cy), "{name} {w}x{h}: G1 of cell ({cx}, {cy})");
                    }
                }
            }
        }
    }

    #[test]
    fn no_pattern_gives_empty_results() {
        let mosaic = indexed_mosaic(4, 4);
        assert!(green_mask(4, 4, BayerPattern::None).is_empty());
        assert_eq!(green_sum(&mosaic, 4, 4, BayerPattern::None), (Vec::new(), 0, 0));
        assert_eq!(green_phase(&mosaic, 4, 4, BayerPattern::None), (Vec::new(), 0, 0));
    }
}
