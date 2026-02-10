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

    for y in 0..out_h {
        for x in 0..out_w {
            let in_y = y * 2;
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

            let idx = y * out_w + x;
            output[idx] = r;
            output[idx + plane_size] = g;
            output[idx + plane_size * 2] = b;
        }
    }

    (output, out_w, out_h)
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

    #[cfg(target_arch = "aarch64")]
    {
        return debayer_f32_neon(input, width, height, pattern);
    }

    #[cfg(target_arch = "x86_64")]
    {
        return debayer_f32_sse2(input, width, height, pattern);
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        debayer_f32_scalar(input, width, height, pattern)
    }
}

#[allow(dead_code)]
fn debayer_f32_scalar(
    input: &[f32],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> (Vec<f32>, usize, usize) {
    let out_w = width / 2;
    let out_h = height / 2;
    let plane_size = out_w * out_h;
    let mut output = vec![0f32; plane_size * 3];

    for y in 0..out_h {
        for x in 0..out_w {
            let in_y = y * 2;
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

            let idx = y * out_w + x;
            output[idx] = r;
            output[idx + plane_size] = g;
            output[idx + plane_size * 2] = b;
        }
    }

    (output, out_w, out_h)
}

#[allow(dead_code)]
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

#[cfg(target_arch = "aarch64")]
fn debayer_f32_neon(
    input: &[f32],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> (Vec<f32>, usize, usize) {
    use std::arch::aarch64::*;

    let out_w = width / 2;
    let out_h = height / 2;
    let plane_size = out_w * out_h;
    let mut output = vec![0f32; plane_size * 3];

    let (ri, gai, gbi, bi) = pattern_indices(pattern);
    let half = unsafe { vdupq_n_f32(0.5) };

    for y in 0..out_h {
        let in_y = y * 2;
        let row0_base = in_y * width;
        let row1_base = (in_y + 1) * width;
        let out_base = y * out_w;
        let mut x = 0;

        unsafe {
            // Process 2 output pixels at a time
            // Each output pixel needs a 2x2 input block = 4 floats (2 per row)
            // For 2 output pixels: 4 floats from row0, 4 floats from row1
            while x + 2 <= out_w {
                let in_x = x * 2;

                // Load 4 consecutive f32 from each row
                let r0 = vld1q_f32(input.as_ptr().add(row0_base + in_x));
                let r1 = vld1q_f32(input.as_ptr().add(row1_base + in_x));

                // r0 = [p00_a, p01_a, p00_b, p01_b]  (2 blocks worth of row0)
                // r1 = [p10_a, p11_a, p10_b, p11_b]  (2 blocks worth of row1)

                // Build a 4-element array: [p00, p01, p10, p11] for each block
                // Deinterleave: evens = [p00_a, p00_b, p10_a, p10_b]
                //                odds = [p01_a, p01_b, p11_a, p11_b]
                let r0_deint = vuzp_f32(vget_low_f32(r0), vget_high_f32(r0));
                let r1_deint = vuzp_f32(vget_low_f32(r1), vget_high_f32(r1));

                // Now we have for each of 2 blocks:
                // p00 = r0_deint.0 lane 0,1  (evens of row0)
                // p01 = r0_deint.1 lane 0,1  (odds of row0)
                // p10 = r1_deint.0 lane 0,1  (evens of row1)
                // p11 = r1_deint.1 lane 0,1  (odds of row1)

                // Combine into arrays indexed by [p00, p01, p10, p11]
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

                // Store to planar output (2 values each)
                vst1_f32(output.as_mut_ptr().add(out_base + x), r_vals);
                vst1_f32(output.as_mut_ptr().add(out_base + x + plane_size), g_vals);
                vst1_f32(output.as_mut_ptr().add(out_base + x + plane_size * 2), b_vals);

                x += 2;
            }
        }

        // Scalar remainder
        for rx in x..out_w {
            let in_y = y * 2;
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

            let idx = y * out_w + rx;
            output[idx] = r;
            output[idx + plane_size] = g;
            output[idx + plane_size * 2] = b;
        }
    }

    (output, out_w, out_h)
}

#[cfg(target_arch = "x86_64")]
fn debayer_f32_sse2(
    input: &[f32],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> (Vec<f32>, usize, usize) {
    use std::arch::x86_64::*;

    let out_w = width / 2;
    let out_h = height / 2;
    let plane_size = out_w * out_h;
    let mut output = vec![0f32; plane_size * 3];

    let (ri, gai, gbi, bi) = pattern_indices(pattern);

    unsafe {
        let v_half = _mm_set1_ps(0.5);

        for y in 0..out_h {
            let in_y = y * 2;
            let row0_base = in_y * width;
            let row1_base = (in_y + 1) * width;
            let out_base = y * out_w;
            let mut x = 0;

            // Process 2 output pixels at a time
            while x + 2 <= out_w {
                let in_x = x * 2;

                // Load 4 floats from each row
                let r0 = _mm_loadu_ps(input.as_ptr().add(row0_base + in_x));
                let r1 = _mm_loadu_ps(input.as_ptr().add(row1_base + in_x));

                // r0 = [p00_a, p01_a, p00_b, p01_b]
                // r1 = [p10_a, p11_a, p10_b, p11_b]

                // Deinterleave: evens and odds
                // p00 = [p00_a, p00_b, ?, ?] = shuffle(r0, r0, 0,2,x,x)
                // p01 = [p01_a, p01_b, ?, ?] = shuffle(r0, r0, 1,3,x,x)
                // p10 = [p10_a, p10_b, ?, ?] = shuffle(r1, r1, 0,2,x,x)
                // p11 = [p11_a, p11_b, ?, ?] = shuffle(r1, r1, 1,3,x,x)
                let p00 = _mm_shuffle_ps(r0, r0, 0b00_00_10_00); // [0,2,0,0]
                let p01 = _mm_shuffle_ps(r0, r0, 0b00_00_11_01); // [1,3,0,0]
                let p10 = _mm_shuffle_ps(r1, r1, 0b00_00_10_00); // [0,2,0,0]
                let p11 = _mm_shuffle_ps(r1, r1, 0b00_00_11_01); // [1,3,0,0]

                let vals = [p00, p01, p10, p11];

                let r_vals = vals[ri];
                let ga_vals = vals[gai];
                let gb_vals = vals[gbi];
                let b_vals = vals[bi];

                let g_vals = _mm_mul_ps(_mm_add_ps(ga_vals, gb_vals), v_half);

                // Store low 2 floats (64 bits) to each plane
                // Cast to f64 and use _mm_store_sd to write 64 bits
                _mm_store_sd(
                    output.as_mut_ptr().add(out_base + x) as *mut f64,
                    _mm_castps_pd(r_vals),
                );
                _mm_store_sd(
                    output.as_mut_ptr().add(out_base + x + plane_size) as *mut f64,
                    _mm_castps_pd(g_vals),
                );
                _mm_store_sd(
                    output.as_mut_ptr().add(out_base + x + plane_size * 2) as *mut f64,
                    _mm_castps_pd(b_vals),
                );

                x += 2;
            }

            // Scalar remainder
            for rx in x..out_w {
                let in_y = y * 2;
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

                let idx = y * out_w + rx;
                output[idx] = r;
                output[idx + plane_size] = g;
                output[idx + plane_size * 2] = b;
            }
        }
    }

    (output, out_w, out_h)
}
