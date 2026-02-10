use rayon::prelude::*;

/// 2x2 averaging for mono preview mode.
pub fn bin_2x2_float(input: &[f32], width: usize, height: usize) -> (Vec<f32>, usize, usize) {
    let out_w = width / 2;
    let out_h = height / 2;
    let mut output = vec![0f32; out_w * out_h];

    output
        .par_chunks_mut(out_w)
        .enumerate()
        .for_each(|(y, out_row)| {
            bin_row(input, width, y, out_row);
        });

    (output, out_w, out_h)
}

fn bin_row(input: &[f32], width: usize, y: usize, out_row: &mut [f32]) {

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe { bin_row_avx2(input, width, y, out_row) };
        } else {
            bin_row_sse2(input, width, y, out_row);
        }
        return;
    }

    #[cfg(target_arch = "aarch64")]
    {
        bin_row_neon(input, width, y, out_row);
        return;
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        bin_row_scalar(input, width, y, out_row);
    }
}

#[allow(dead_code)]
fn bin_row_scalar(input: &[f32], width: usize, y: usize, out_row: &mut [f32]) {
    let out_w = out_row.len();
    let row0 = (y * 2) * width;
    let row1 = (y * 2 + 1) * width;
    for x in 0..out_w {
        let in_x = x * 2;
        out_row[x] = (input[row0 + in_x]
            + input[row0 + in_x + 1]
            + input[row1 + in_x]
            + input[row1 + in_x + 1])
            * 0.25;
    }
}

#[cfg(target_arch = "x86_64")]
fn bin_row_sse2(input: &[f32], width: usize, y: usize, out_row: &mut [f32]) {
    use std::arch::x86_64::*;

    let out_w = out_row.len();
    let row0 = (y * 2) * width;
    let row1 = (y * 2 + 1) * width;
    let mut x = 0;

    unsafe {
        let v_quarter = _mm_set1_ps(0.25);

        while x + 4 <= out_w {
            let in_x = x * 2;

            let r0a = _mm_loadu_ps(input.as_ptr().add(row0 + in_x));
            let r0b = _mm_loadu_ps(input.as_ptr().add(row0 + in_x + 4));
            let r1a = _mm_loadu_ps(input.as_ptr().add(row1 + in_x));
            let r1b = _mm_loadu_ps(input.as_ptr().add(row1 + in_x + 4));

            let r0_even = _mm_shuffle_ps(r0a, r0b, 0b10_00_10_00);
            let r0_odd = _mm_shuffle_ps(r0a, r0b, 0b11_01_11_01);
            let r1_even = _mm_shuffle_ps(r1a, r1b, 0b10_00_10_00);
            let r1_odd = _mm_shuffle_ps(r1a, r1b, 0b11_01_11_01);

            let sum = _mm_add_ps(
                _mm_add_ps(r0_even, r0_odd),
                _mm_add_ps(r1_even, r1_odd),
            );
            let avg = _mm_mul_ps(sum, v_quarter);

            _mm_storeu_ps(out_row.as_mut_ptr().add(x), avg);
            x += 4;
        }
    }

    // Scalar remainder
    for rx in x..out_w {
        let in_x = rx * 2;
        out_row[rx] = (input[row0 + in_x]
            + input[row0 + in_x + 1]
            + input[row1 + in_x]
            + input[row1 + in_x + 1])
            * 0.25;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn bin_row_avx2(input: &[f32], width: usize, y: usize, out_row: &mut [f32]) {
    use std::arch::x86_64::*;

    let out_w = out_row.len();
    let row0 = (y * 2) * width;
    let row1 = (y * 2 + 1) * width;
    let mut x = 0;

    let v_quarter = _mm256_set1_ps(0.25);

    while x + 8 <= out_w {
        let in_x = x * 2;

        let r0a = _mm256_loadu_ps(input.as_ptr().add(row0 + in_x));
        let r0b = _mm256_loadu_ps(input.as_ptr().add(row0 + in_x + 8));
        let r1a = _mm256_loadu_ps(input.as_ptr().add(row1 + in_x));
        let r1b = _mm256_loadu_ps(input.as_ptr().add(row1 + in_x + 8));

        let r0_even = _mm256_shuffle_ps(r0a, r0b, 0b10_00_10_00);
        let r0_odd = _mm256_shuffle_ps(r0a, r0b, 0b11_01_11_01);
        let r1_even = _mm256_shuffle_ps(r1a, r1b, 0b10_00_10_00);
        let r1_odd = _mm256_shuffle_ps(r1a, r1b, 0b11_01_11_01);

        let r0_even = _mm256_permute4x64_pd(_mm256_castps_pd(r0_even), 0b11_01_10_00);
        let r0_even = _mm256_castpd_ps(r0_even);
        let r0_odd = _mm256_permute4x64_pd(_mm256_castps_pd(r0_odd), 0b11_01_10_00);
        let r0_odd = _mm256_castpd_ps(r0_odd);
        let r1_even = _mm256_permute4x64_pd(_mm256_castps_pd(r1_even), 0b11_01_10_00);
        let r1_even = _mm256_castpd_ps(r1_even);
        let r1_odd = _mm256_permute4x64_pd(_mm256_castps_pd(r1_odd), 0b11_01_10_00);
        let r1_odd = _mm256_castpd_ps(r1_odd);

        let sum = _mm256_add_ps(
            _mm256_add_ps(r0_even, r0_odd),
            _mm256_add_ps(r1_even, r1_odd),
        );
        let avg = _mm256_mul_ps(sum, v_quarter);

        _mm256_storeu_ps(out_row.as_mut_ptr().add(x), avg);
        x += 8;
    }

    // Scalar remainder
    for rx in x..out_w {
        let in_x = rx * 2;
        *out_row.get_unchecked_mut(rx) = (*input.get_unchecked(row0 + in_x)
            + *input.get_unchecked(row0 + in_x + 1)
            + *input.get_unchecked(row1 + in_x)
            + *input.get_unchecked(row1 + in_x + 1))
            * 0.25;
    }
}

#[cfg(target_arch = "aarch64")]
fn bin_row_neon(input: &[f32], width: usize, y: usize, out_row: &mut [f32]) {
    use std::arch::aarch64::*;

    let out_w = out_row.len();
    let row0 = (y * 2) * width;
    let row1 = (y * 2 + 1) * width;
    let mut x = 0;

    unsafe {
        let v_quarter = vdupq_n_f32(0.25);

        while x + 4 <= out_w {
            let in_x = x * 2;

            let r0 = vld2q_f32(input.as_ptr().add(row0 + in_x));
            let r1 = vld2q_f32(input.as_ptr().add(row1 + in_x));

            let sum = vaddq_f32(
                vaddq_f32(r0.0, r0.1),
                vaddq_f32(r1.0, r1.1),
            );
            let avg = vmulq_f32(sum, v_quarter);

            vst1q_f32(out_row.as_mut_ptr().add(x), avg);
            x += 4;
        }
    }

    // Scalar remainder
    for rx in x..out_w {
        let in_x = rx * 2;
        out_row[rx] = (input[row0 + in_x]
            + input[row0 + in_x + 1]
            + input[row1 + in_x]
            + input[row1 + in_x + 1])
            * 0.25;
    }
}
