/// Quickselect: find k-th smallest element (in-place, modifies slice).
fn quickselect(arr: &mut [f32], k: usize) -> f32 {
    let mut left = 0usize;
    let mut right = arr.len() - 1;

    while left < right {
        let mid = left + (right - left) / 2;

        // Median-of-three pivot
        if arr[mid] < arr[left] {
            arr.swap(left, mid);
        }
        if arr[right] < arr[left] {
            arr.swap(left, right);
        }
        if arr[right] < arr[mid] {
            arr.swap(mid, right);
        }

        let pivot = arr[mid];

        // Move pivot to right-1
        arr.swap(mid, right - 1);

        // Partition
        let mut i = left;
        let mut j = right - 1;

        loop {
            i += 1;
            while arr[i] < pivot {
                i += 1;
            }
            j -= 1;
            while arr[j] > pivot {
                j -= 1;
            }
            if i >= j {
                break;
            }
            arr.swap(i, j);
        }

        // Restore pivot
        arr.swap(i, right - 1);

        if i == k {
            return arr[k];
        } else if i > k {
            right = i - 1;
        } else {
            left = i + 1;
        }
    }

    arr[k]
}

fn find_median(data: &mut [f32]) -> f32 {
    let k = data.len() / 2;
    quickselect(data, k)
}

pub struct StretchParams {
    pub shadows: f32,
    pub highlights: f32,
    pub midtones: f32,
}

pub fn compute_stretch_params(data: &[f32], max_input: f32) -> StretchParams {
    const MAX_SAMPLES: usize = 500_000;

    let num_samples = data.len().min(MAX_SAMPLES);
    let mut samples = if data.len() <= MAX_SAMPLES {
        data.to_vec()
    } else {
        let step = data.len() / MAX_SAMPLES;
        (0..num_samples).map(|i| data[i * step]).collect()
    };

    let median = find_median(&mut samples);

    // Compute MADN
    let mut deviations: Vec<f32> = samples.iter().map(|&v| (v - median).abs()).collect();
    let madn = 1.4826 * find_median(&mut deviations);

    // Normalize
    let norm_median = median / max_input;
    let norm_madn = madn / max_input;

    let upper_half = norm_median > 0.5;

    // Shadows
    let shadows = if upper_half || norm_madn == 0.0 {
        0.0
    } else {
        (norm_median + (-2.8 * norm_madn)).clamp(0.0, 1.0)
    };

    // Highlights
    let highlights = if !upper_half || norm_madn == 0.0 {
        1.0
    } else {
        (norm_median - (-2.8 * norm_madn)).clamp(0.0, 1.0)
    };

    // Midtones via QuickFits formula
    let b = 0.25f32;
    let (x, m) = if !upper_half {
        (norm_median - shadows, b)
    } else {
        (b, highlights - norm_median)
    };

    let midtones = if x == 0.0 {
        0.0
    } else if x == m {
        0.5
    } else if x == 1.0 {
        1.0
    } else {
        ((m - 1.0) * x) / ((2.0 * m - 1.0) * x - m)
    };

    StretchParams {
        shadows,
        highlights,
        midtones,
    }
}

/// Apply stretch to a channel, writing into `output` at positions `output_offset + i * stride`.
pub fn apply_stretch(
    channel_data: &[f32],
    output: &mut [u8],
    output_offset: usize,
    stride: usize,
    native_shadows: f32,
    native_highlights: f32,
    k1: f32,
    k2: f32,
    midtones: f32,
) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                apply_stretch_avx2(
                    channel_data,
                    output,
                    output_offset,
                    stride,
                    native_shadows,
                    native_highlights,
                    k1,
                    k2,
                    midtones,
                );
            }
        } else {
            apply_stretch_sse2(
                channel_data,
                output,
                output_offset,
                stride,
                native_shadows,
                native_highlights,
                k1,
                k2,
                midtones,
            );
        }
        return;
    }

    #[cfg(target_arch = "aarch64")]
    {
        apply_stretch_neon(
            channel_data,
            output,
            output_offset,
            stride,
            native_shadows,
            native_highlights,
            k1,
            k2,
            midtones,
        );
        return;
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        apply_stretch_scalar(
            channel_data,
            output,
            output_offset,
            stride,
            native_shadows,
            native_highlights,
            k1,
            k2,
            midtones,
        );
    }
}

fn stretch_pixel(input: f32, native_shadows: f32, native_highlights: f32, k1: f32, k2: f32, midtones: f32) -> u8 {
    let out = if input < native_shadows {
        0.0f32
    } else if input >= native_highlights {
        255.0f32
    } else {
        let input_floored = input - native_shadows;
        (input_floored * k1) / (input_floored * k2 - midtones)
    };
    out.clamp(0.0, 255.0) as u8
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
fn apply_stretch_scalar(
    channel_data: &[f32],
    output: &mut [u8],
    output_offset: usize,
    stride: usize,
    native_shadows: f32,
    native_highlights: f32,
    k1: f32,
    k2: f32,
    midtones: f32,
) {
    for (i, &input) in channel_data.iter().enumerate() {
        output[output_offset + i * stride] =
            stretch_pixel(input, native_shadows, native_highlights, k1, k2, midtones);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn apply_stretch_avx2(
    channel_data: &[f32],
    output: &mut [u8],
    output_offset: usize,
    stride: usize,
    native_shadows: f32,
    native_highlights: f32,
    k1: f32,
    k2: f32,
    midtones: f32,
) {
    use std::arch::x86_64::*;

    let count = channel_data.len();
    let mut i = 0;

    let v_shadows = _mm256_set1_ps(native_shadows);
    let v_highlights = _mm256_set1_ps(native_highlights);
    let v_k1 = _mm256_set1_ps(k1);
    let v_k2 = _mm256_set1_ps(k2);
    let v_midtones = _mm256_set1_ps(midtones);
    let v_zero = _mm256_setzero_ps();
    let v_255 = _mm256_set1_ps(255.0);

    while i + 8 <= count {
        let input = _mm256_loadu_ps(channel_data.as_ptr().add(i));

        // Masks: CMP returns all-ones or all-zeros per lane
        let mask_low = _mm256_cmp_ps(input, v_shadows, _CMP_LT_OS);
        let mask_high = _mm256_cmp_ps(input, v_highlights, _CMP_GE_OS);

        // Midtones transfer
        let input_floored = _mm256_sub_ps(input, v_shadows);
        let numerator = _mm256_mul_ps(input_floored, v_k1);
        let denominator = _mm256_sub_ps(_mm256_mul_ps(input_floored, v_k2), v_midtones);
        let output_val = _mm256_div_ps(numerator, denominator);

        // Blend: 0 if < shadows, 255 if >= highlights, computed otherwise
        let result = _mm256_blendv_ps(output_val, v_zero, mask_low);
        let result = _mm256_blendv_ps(result, v_255, mask_high);

        // Clamp to [0, 255]
        let clamped = _mm256_min_ps(_mm256_max_ps(result, v_zero), v_255);

        // Convert f32 → i32 → pack to u8
        let output_int = _mm256_cvtps_epi32(clamped);
        // Pack 8x i32 → 8x i16 (with saturation)
        // _mm256_packs_epi32 packs within 128-bit lanes: [0..3,4..7] → [0..3 as i16, 4..7 as i16] per lane
        let packed16 = _mm256_packs_epi32(output_int, output_int);
        // Pack 8x i16 → 8x u8 (with saturation)
        let packed8 = _mm256_packus_epi16(packed16, packed16);

        // Extract the 8 bytes: they're at bytes [0..3] and [16..19] due to lane interleaving
        // Use permute to bring lane1 data adjacent to lane0
        let shuffled = _mm256_permute4x64_epi64(packed8, 0b11_01_10_00);
        let lo128 = _mm256_castsi256_si128(shuffled);
        let packed = _mm_cvtsi128_si64(lo128) as u64;

        for j in 0..8 {
            output[output_offset + (i + j) * stride] = ((packed >> (j * 8)) & 0xFF) as u8;
        }

        i += 8;
    }

    // Scalar remainder
    for idx in i..count {
        output[output_offset + idx * stride] =
            stretch_pixel(channel_data[idx], native_shadows, native_highlights, k1, k2, midtones);
    }
}

#[cfg(target_arch = "x86_64")]
fn apply_stretch_sse2(
    channel_data: &[f32],
    output: &mut [u8],
    output_offset: usize,
    stride: usize,
    native_shadows: f32,
    native_highlights: f32,
    k1: f32,
    k2: f32,
    midtones: f32,
) {
    use std::arch::x86_64::*;

    let count = channel_data.len();
    let mut i = 0;

    unsafe {
        let v_shadows = _mm_set1_ps(native_shadows);
        let v_highlights = _mm_set1_ps(native_highlights);
        let v_k1 = _mm_set1_ps(k1);
        let v_k2 = _mm_set1_ps(k2);
        let v_midtones = _mm_set1_ps(midtones);
        let v_zero = _mm_setzero_ps();
        let v_255 = _mm_set1_ps(255.0);

        while i + 4 <= count {
            let input = _mm_loadu_ps(channel_data.as_ptr().add(i));

            let mask_low = _mm_cmplt_ps(input, v_shadows);
            let mask_high = _mm_cmpge_ps(input, v_highlights);

            let input_floored = _mm_sub_ps(input, v_shadows);
            let numerator = _mm_mul_ps(input_floored, v_k1);
            let denominator = _mm_sub_ps(_mm_mul_ps(input_floored, v_k2), v_midtones);
            let output_val = _mm_div_ps(numerator, denominator);

            // Blend: 0 if < shadows, 255 if >= highlights, computed otherwise
            let all_ones = _mm_castsi128_ps(_mm_set1_epi32(-1));
            let mid_mask = _mm_andnot_ps(mask_low, _mm_andnot_ps(mask_high, all_ones));
            let result = _mm_or_ps(
                _mm_and_ps(output_val, mid_mask),
                _mm_and_ps(v_255, mask_high),
            );

            let clamped = _mm_min_ps(_mm_max_ps(result, v_zero), v_255);
            let output_int = _mm_cvtps_epi32(clamped);
            let packed16 = _mm_packs_epi32(output_int, output_int);
            let packed8 = _mm_packus_epi16(packed16, packed16);
            let packed = _mm_cvtsi128_si32(packed8) as u32;

            for j in 0..4 {
                output[output_offset + (i + j) * stride] = ((packed >> (j * 8)) & 0xFF) as u8;
            }

            i += 4;
        }
    }

    // Scalar remainder
    for idx in i..count {
        output[output_offset + idx * stride] =
            stretch_pixel(channel_data[idx], native_shadows, native_highlights, k1, k2, midtones);
    }
}

#[cfg(target_arch = "aarch64")]
fn apply_stretch_neon(
    channel_data: &[f32],
    output: &mut [u8],
    output_offset: usize,
    stride: usize,
    native_shadows: f32,
    native_highlights: f32,
    k1: f32,
    k2: f32,
    midtones: f32,
) {
    use std::arch::aarch64::*;

    let count = channel_data.len();
    let mut i = 0;

    unsafe {
        let v_shadows = vdupq_n_f32(native_shadows);
        let v_highlights = vdupq_n_f32(native_highlights);
        let v_k1 = vdupq_n_f32(k1);
        let v_k2 = vdupq_n_f32(k2);
        let v_midtones = vdupq_n_f32(midtones);
        let v_zero = vdupq_n_f32(0.0);
        let v_255 = vdupq_n_f32(255.0);

        while i + 4 <= count {
            let input = vld1q_f32(channel_data.as_ptr().add(i));

            let mask_low = vcltq_f32(input, v_shadows);
            let mask_high = vcgeq_f32(input, v_highlights);

            let input_floored = vsubq_f32(input, v_shadows);
            let numerator = vmulq_f32(input_floored, v_k1);
            let denominator = vsubq_f32(vmulq_f32(input_floored, v_k2), v_midtones);

            // Newton-Raphson reciprocal
            let recip = vrecpeq_f32(denominator);
            let recip = vmulq_f32(vrecpsq_f32(denominator, recip), recip);
            let output_val = vmulq_f32(numerator, recip);

            let result = vbslq_f32(mask_low, v_zero, output_val);
            let result = vbslq_f32(mask_high, v_255, result);

            let clamped = vminq_f32(vmaxq_f32(result, v_zero), v_255);
            let output_int = vcvtq_u32_f32(clamped);
            let output_16 = vmovn_u32(output_int);
            let output_8 = vmovn_u16(vcombine_u16(output_16, output_16));

            let mut packed = [0u8; 8];
            vst1_u8(packed.as_mut_ptr(), output_8);

            for j in 0..4 {
                output[output_offset + (i + j) * stride] = packed[j];
            }

            i += 4;
        }
    }

    // Scalar remainder
    for idx in i..count {
        output[output_offset + idx * stride] =
            stretch_pixel(channel_data[idx], native_shadows, native_highlights, k1, k2, midtones);
    }
}
