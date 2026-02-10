use rayon::prelude::*;

use crate::types::BayerPattern;

/// Super-pixel 2x2 debayer for u16 input. Output is planar f32 RGB at half resolution.
pub fn super_pixel_debayer_u16(
    input: &[u16],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> (Vec<f32>, usize, usize) {
    let out_w = width / 2;
    let out_h = height / 2;
    let plane_size = out_w * out_h;
    let mut output = vec![0f32; plane_size * 3];

    let (r_plane, rest) = output.split_at_mut(plane_size);
    let (g_plane, b_plane) = rest.split_at_mut(plane_size);

    r_plane
        .par_chunks_mut(out_w)
        .zip(g_plane.par_chunks_mut(out_w))
        .zip(b_plane.par_chunks_mut(out_w))
        .enumerate()
        .for_each(|(y, ((r_row, g_row), b_row))| {
            debayer_u16_row(input, width, y, out_w, r_row, g_row, b_row, pattern);
        });

    (output, out_w, out_h)
}

fn debayer_u16_row(
    input: &[u16],
    width: usize,
    y: usize,
    out_w: usize,
    r_row: &mut [f32],
    g_row: &mut [f32],
    b_row: &mut [f32],
    pattern: BayerPattern,
) {
    let in_y = y * 2;
    for x in 0..out_w {
        let in_x = x * 2;

        let p00 = input[in_y * width + in_x] as f32;
        let p01 = input[in_y * width + in_x + 1] as f32;
        let p10 = input[(in_y + 1) * width + in_x] as f32;
        let p11 = input[(in_y + 1) * width + in_x + 1] as f32;

        let (r, g, b) = match pattern {
            BayerPattern::Rggb => (p00, (p01 + p10) * 0.5, p11),
            BayerPattern::Bggr => (p11, (p01 + p10) * 0.5, p00),
            BayerPattern::Gbrg => (p10, (p00 + p11) * 0.5, p01),
            BayerPattern::Grbg => (p01, (p00 + p11) * 0.5, p10),
            BayerPattern::None => (0.0, 0.0, 0.0),
        };

        r_row[x] = r;
        g_row[x] = g;
        b_row[x] = b;
    }
}

/// Super-pixel 2x2 debayer for f32 input. Output is planar f32 RGB at half resolution.
pub fn super_pixel_debayer_f32(
    input: &[f32],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> (Vec<f32>, usize, usize) {
    if pattern == BayerPattern::None {
        let out_w = width / 2;
        let out_h = height / 2;
        return (vec![0f32; out_w * out_h * 3], out_w, out_h);
    }

    let out_w = width / 2;
    let out_h = height / 2;
    let plane_size = out_w * out_h;
    let mut output = vec![0f32; plane_size * 3];

    let (r_plane, rest) = output.split_at_mut(plane_size);
    let (g_plane, b_plane) = rest.split_at_mut(plane_size);

    r_plane
        .par_chunks_mut(out_w)
        .zip(g_plane.par_chunks_mut(out_w))
        .zip(b_plane.par_chunks_mut(out_w))
        .enumerate()
        .for_each(|(y, ((r_row, g_row), b_row))| {
            debayer_f32_row(input, width, y, out_w, r_row, g_row, b_row, pattern);
        });

    (output, out_w, out_h)
}

/// Map Bayer pattern to indices: for each 2x2 block [p00, p01, p10, p11],
/// returns (r_idx, ga_idx, gb_idx, b_idx) where ga and gb are the two green positions.
/// Green = (ga + gb) * 0.5
fn pattern_indices(pattern: BayerPattern) -> (usize, usize, usize, usize) {
    match pattern {
        // (r, ga, gb, b) positions in [p00=0, p01=1, p10=2, p11=3]
        BayerPattern::Rggb => (0, 1, 2, 3), // R=p00, G=(p01+p10)/2, B=p11
        BayerPattern::Bggr => (3, 1, 2, 0), // R=p11, G=(p01+p10)/2, B=p00
        BayerPattern::Gbrg => (2, 0, 3, 1), // R=p10, G=(p00+p11)/2, B=p01
        BayerPattern::Grbg => (1, 0, 3, 2), // R=p01, G=(p00+p11)/2, B=p10
        BayerPattern::None => (0, 0, 0, 0),
    }
}

/// Per-row debayer dispatch for f32 input.
fn debayer_f32_row(
    input: &[f32],
    width: usize,
    y: usize,
    out_w: usize,
    r_row: &mut [f32],
    g_row: &mut [f32],
    b_row: &mut [f32],
    pattern: BayerPattern,
) {
    #[cfg(target_arch = "aarch64")]
    {
        debayer_f32_neon_row(input, width, y, out_w, r_row, g_row, b_row, pattern);
        return;
    }

    #[cfg(target_arch = "x86_64")]
    {
        debayer_f32_sse2_row(input, width, y, out_w, r_row, g_row, b_row, pattern);
        return;
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        debayer_f32_scalar_row(input, width, y, out_w, r_row, g_row, b_row, pattern);
    }
}

#[allow(dead_code)]
fn debayer_f32_scalar_row(
    input: &[f32],
    width: usize,
    y: usize,
    out_w: usize,
    r_row: &mut [f32],
    g_row: &mut [f32],
    b_row: &mut [f32],
    pattern: BayerPattern,
) {
    let in_y = y * 2;
    for x in 0..out_w {
        let in_x = x * 2;

        let p00 = input[in_y * width + in_x];
        let p01 = input[in_y * width + in_x + 1];
        let p10 = input[(in_y + 1) * width + in_x];
        let p11 = input[(in_y + 1) * width + in_x + 1];

        let (r, g, b) = match pattern {
            BayerPattern::Rggb => (p00, (p01 + p10) * 0.5, p11),
            BayerPattern::Bggr => (p11, (p01 + p10) * 0.5, p00),
            BayerPattern::Gbrg => (p10, (p00 + p11) * 0.5, p01),
            BayerPattern::Grbg => (p01, (p00 + p11) * 0.5, p10),
            BayerPattern::None => unreachable!(),
        };

        r_row[x] = r;
        g_row[x] = g;
        b_row[x] = b;
    }
}

#[cfg(target_arch = "aarch64")]
fn debayer_f32_neon_row(
    input: &[f32],
    width: usize,
    y: usize,
    out_w: usize,
    r_row: &mut [f32],
    g_row: &mut [f32],
    b_row: &mut [f32],
    pattern: BayerPattern,
) {
    use std::arch::aarch64::*;

    let (ri, gai, gbi, bi) = pattern_indices(pattern);
    let half = unsafe { vdupq_n_f32(0.5) };

    let in_y = y * 2;
    let row0_base = in_y * width;
    let row1_base = (in_y + 1) * width;
    let mut x = 0;

    unsafe {
        while x + 2 <= out_w {
            let in_x = x * 2;

            let r0 = vld1q_f32(input.as_ptr().add(row0_base + in_x));
            let r1 = vld1q_f32(input.as_ptr().add(row1_base + in_x));

            let r0_deint = vuzp_f32(vget_low_f32(r0), vget_high_f32(r0));
            let r1_deint = vuzp_f32(vget_low_f32(r1), vget_high_f32(r1));

            let vals: [float64x1_t; 4] = [
                vreinterpret_f64_f32(r0_deint.0),
                vreinterpret_f64_f32(r0_deint.1),
                vreinterpret_f64_f32(r1_deint.0),
                vreinterpret_f64_f32(r1_deint.1),
            ];

            let r_vals = vreinterpret_f32_f64(vals[ri]);
            let ga_vals = vreinterpret_f32_f64(vals[gai]);
            let gb_vals = vreinterpret_f32_f64(vals[gbi]);
            let b_vals = vreinterpret_f32_f64(vals[bi]);

            let g_vals = vmul_f32(vadd_f32(ga_vals, gb_vals), vget_low_f32(half));

            vst1_f32(r_row.as_mut_ptr().add(x), r_vals);
            vst1_f32(g_row.as_mut_ptr().add(x), g_vals);
            vst1_f32(b_row.as_mut_ptr().add(x), b_vals);

            x += 2;
        }
    }

    // Scalar remainder
    let in_y = y * 2;
    for rx in x..out_w {
        let in_x = rx * 2;
        let p00 = input[in_y * width + in_x];
        let p01 = input[in_y * width + in_x + 1];
        let p10 = input[(in_y + 1) * width + in_x];
        let p11 = input[(in_y + 1) * width + in_x + 1];

        let (r, g, b) = match pattern {
            BayerPattern::Rggb => (p00, (p01 + p10) * 0.5, p11),
            BayerPattern::Bggr => (p11, (p01 + p10) * 0.5, p00),
            BayerPattern::Gbrg => (p10, (p00 + p11) * 0.5, p01),
            BayerPattern::Grbg => (p01, (p00 + p11) * 0.5, p10),
            BayerPattern::None => unreachable!(),
        };

        r_row[rx] = r;
        g_row[rx] = g;
        b_row[rx] = b;
    }
}

#[cfg(target_arch = "x86_64")]
fn debayer_f32_sse2_row(
    input: &[f32],
    width: usize,
    y: usize,
    out_w: usize,
    r_row: &mut [f32],
    g_row: &mut [f32],
    b_row: &mut [f32],
    pattern: BayerPattern,
) {
    use std::arch::x86_64::*;

    let (ri, gai, gbi, bi) = pattern_indices(pattern);

    let in_y = y * 2;
    let row0_base = in_y * width;
    let row1_base = (in_y + 1) * width;
    let mut x = 0;

    unsafe {
        let v_half = _mm_set1_ps(0.5);

        while x + 2 <= out_w {
            let in_x = x * 2;

            let r0 = _mm_loadu_ps(input.as_ptr().add(row0_base + in_x));
            let r1 = _mm_loadu_ps(input.as_ptr().add(row1_base + in_x));

            let p00 = _mm_shuffle_ps(r0, r0, 0b00_00_10_00);
            let p01 = _mm_shuffle_ps(r0, r0, 0b00_00_11_01);
            let p10 = _mm_shuffle_ps(r1, r1, 0b00_00_10_00);
            let p11 = _mm_shuffle_ps(r1, r1, 0b00_00_11_01);

            let vals = [p00, p01, p10, p11];

            let r_vals = vals[ri];
            let ga_vals = vals[gai];
            let gb_vals = vals[gbi];
            let b_vals = vals[bi];

            let g_vals = _mm_mul_ps(_mm_add_ps(ga_vals, gb_vals), v_half);

            _mm_store_sd(
                r_row.as_mut_ptr().add(x) as *mut f64,
                _mm_castps_pd(r_vals),
            );
            _mm_store_sd(
                g_row.as_mut_ptr().add(x) as *mut f64,
                _mm_castps_pd(g_vals),
            );
            _mm_store_sd(
                b_row.as_mut_ptr().add(x) as *mut f64,
                _mm_castps_pd(b_vals),
            );

            x += 2;
        }
    }

    // Scalar remainder
    let in_y = y * 2;
    for rx in x..out_w {
        let in_x = rx * 2;
        let p00 = input[in_y * width + in_x];
        let p01 = input[in_y * width + in_x + 1];
        let p10 = input[(in_y + 1) * width + in_x];
        let p11 = input[(in_y + 1) * width + in_x + 1];

        let (r, g, b) = match pattern {
            BayerPattern::Rggb => (p00, (p01 + p10) * 0.5, p11),
            BayerPattern::Bggr => (p11, (p01 + p10) * 0.5, p00),
            BayerPattern::Gbrg => (p10, (p00 + p11) * 0.5, p01),
            BayerPattern::Grbg => (p01, (p00 + p11) * 0.5, p10),
            BayerPattern::None => unreachable!(),
        };

        r_row[rx] = r;
        g_row[rx] = g;
        b_row[rx] = b;
    }
}
