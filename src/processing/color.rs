/// Replicate grayscale u8 to interleaved RGB u8.
pub fn replicate_gray_to_rgb(gray: &[u8]) -> Vec<u8> {
    let mut rgb = vec![0u8; gray.len() * 3];

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe { gray_to_rgb_avx2(gray, &mut rgb) };
        } else {
            gray_to_rgb_sse2(gray, &mut rgb);
        }
        return rgb;
    }

    #[cfg(target_arch = "aarch64")]
    {
        gray_to_rgb_neon(gray, &mut rgb);
        return rgb;
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        gray_to_rgb_scalar(gray, &mut rgb);
        rgb
    }
}

#[allow(dead_code)]
fn gray_to_rgb_scalar(gray: &[u8], rgb: &mut [u8]) {
    for (i, &val) in gray.iter().enumerate() {
        rgb[i * 3] = val;
        rgb[i * 3 + 1] = val;
        rgb[i * 3 + 2] = val;
    }
}

#[cfg(target_arch = "aarch64")]
fn gray_to_rgb_neon(gray: &[u8], rgb: &mut [u8]) {
    use std::arch::aarch64::*;

    let count = gray.len();
    let mut i = 0;

    unsafe {
        // Process 16 pixels at a time (reads 16 bytes, writes 48 bytes)
        while i + 16 <= count {
            let g = vld1q_u8(gray.as_ptr().add(i));
            let triplet = uint8x16x3_t(g, g, g);
            vst3q_u8(rgb.as_mut_ptr().add(i * 3), triplet);
            i += 16;
        }
    }

    // Scalar remainder
    for idx in i..count {
        rgb[idx * 3] = gray[idx];
        rgb[idx * 3 + 1] = gray[idx];
        rgb[idx * 3 + 2] = gray[idx];
    }
}

// SSE2 doesn't have _mm_shuffle_epi8 (that's SSSE3), so gray→RGB
// on baseline x86_64 uses scalar. The AVX2 path (which implies SSSE3)
// uses pshufb for the real SIMD acceleration.
#[cfg(target_arch = "x86_64")]
fn gray_to_rgb_sse2(gray: &[u8], rgb: &mut [u8]) {
    gray_to_rgb_scalar(gray, rgb);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn gray_to_rgb_avx2(gray: &[u8], rgb: &mut [u8]) {
    use std::arch::x86_64::*;

    let count = gray.len();
    let mut i = 0;

    // With AVX2 we have _mm256_shuffle_epi8 (SSSE3 + AVX2)
    // Process 32 input pixels → 96 output bytes at a time
    // Strategy: load 32 bytes, use shuffle to create 3x32 byte blocks

    // Shuffle masks to triplicate bytes within 128-bit lanes
    // For first 16 output bytes from first 5-6 input bytes:
    // We process 10 pixels at a time from a 128-bit lane (30 output bytes fits in ~2 stores)

    // Simpler approach: process 16 pixels per iteration using __m128i + pshufb
    let shuf0 = _mm_setr_epi8(0,0,0, 1,1,1, 2,2,2, 3,3,3, 4,4,4, 5);
    let shuf1 = _mm_setr_epi8(5,5, 6,6,6, 7,7,7, 8,8,8, 9,9,9, 10,10);
    let shuf2 = _mm_setr_epi8(10, 11,11,11, 12,12,12, 13,13,13, 14,14,14, 15,15,15);

    while i + 16 <= count {
        let g = _mm_loadu_si128(gray.as_ptr().add(i) as *const __m128i);

        // pshufb to create triplets
        let out0 = _mm_shuffle_epi8(g, shuf0); // bytes 0-15 of output (pixels 0-5)
        let out1 = _mm_shuffle_epi8(g, shuf1); // bytes 16-31 of output (pixels 5-10)
        let out2 = _mm_shuffle_epi8(g, shuf2); // bytes 32-47 of output (pixels 10-15)

        let base = i * 3;
        _mm_storeu_si128(rgb.as_mut_ptr().add(base) as *mut __m128i, out0);
        _mm_storeu_si128(rgb.as_mut_ptr().add(base + 16) as *mut __m128i, out1);
        _mm_storeu_si128(rgb.as_mut_ptr().add(base + 32) as *mut __m128i, out2);

        i += 16;
    }

    // Scalar remainder
    for idx in i..count {
        let val = gray[idx];
        let out = idx * 3;
        *rgb.get_unchecked_mut(out) = val;
        *rgb.get_unchecked_mut(out + 1) = val;
        *rgb.get_unchecked_mut(out + 2) = val;
    }
}

/// Vertical flip of an RGB image (3 bytes per pixel).
pub fn vertical_flip_rgb(data: &mut [u8], width: usize, height: usize) {
    let row_bytes = width * 3;
    let mut temp = vec![0u8; row_bytes];
    for y in 0..height / 2 {
        let top = y * row_bytes;
        let bot = (height - 1 - y) * row_bytes;
        temp.copy_from_slice(&data[top..top + row_bytes]);
        data.copy_within(bot..bot + row_bytes, top);
        data[bot..bot + row_bytes].copy_from_slice(&temp);
    }
}

/// Convert u16 slice to f32 (parallel with SIMD dispatch per chunk).
pub fn u16_to_f32(data: &[u16]) -> Vec<f32> {
    use rayon::prelude::*;

    let mut output = vec![0f32; data.len()];

    const CHUNK: usize = 65536;
    data.par_chunks(CHUNK)
        .zip(output.par_chunks_mut(CHUNK))
        .for_each(|(src, dst)| {
            u16_to_f32_chunk(src, dst);
        });

    output
}

/// SIMD-dispatched u16→f32 conversion for a chunk.
fn u16_to_f32_chunk(data: &[u16], output: &mut [f32]) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe { u16_to_f32_avx2(data, output) };
        } else {
            u16_to_f32_sse2(data, output);
        }
        return;
    }

    #[cfg(target_arch = "aarch64")]
    {
        u16_to_f32_neon(data, output);
        return;
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        for (i, &v) in data.iter().enumerate() {
            output[i] = v as f32;
        }
    }
}

#[cfg(target_arch = "aarch64")]
fn u16_to_f32_neon(data: &[u16], output: &mut [f32]) {
    use std::arch::aarch64::*;

    let count = data.len();
    let mut i = 0;

    unsafe {
        // Process 4 u16 at a time
        while i + 4 <= count {
            let u16x4 = vld1_u16(data.as_ptr().add(i));
            let u32x4 = vmovl_u16(u16x4);
            let f32x4 = vcvtq_f32_u32(u32x4);
            vst1q_f32(output.as_mut_ptr().add(i), f32x4);
            i += 4;
        }
    }

    // Scalar remainder
    for idx in i..count {
        output[idx] = data[idx] as f32;
    }
}

#[cfg(target_arch = "x86_64")]
fn u16_to_f32_sse2(data: &[u16], output: &mut [f32]) {
    use std::arch::x86_64::*;

    let count = data.len();
    let mut i = 0;

    unsafe {
        // Process 4 u16 at a time
        while i + 4 <= count {
            // Load 4 u16 (64 bits)
            let u16x4 = _mm_loadl_epi64(data.as_ptr().add(i) as *const __m128i);
            // Zero-extend u16 → u32
            let u32x4 = _mm_unpacklo_epi16(u16x4, _mm_setzero_si128());
            // Convert i32 → f32 (values are non-negative, so this is correct)
            let f32x4 = _mm_cvtepi32_ps(u32x4);
            _mm_storeu_ps(output.as_mut_ptr().add(i), f32x4);
            i += 4;
        }
    }

    // Scalar remainder
    for idx in i..count {
        output[idx] = data[idx] as f32;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn u16_to_f32_avx2(data: &[u16], output: &mut [f32]) {
    use std::arch::x86_64::*;

    let count = data.len();
    let mut i = 0;

    // Process 8 u16 at a time
    while i + 8 <= count {
        // Load 8 u16 (128 bits)
        let u16x8 = _mm_loadu_si128(data.as_ptr().add(i) as *const __m128i);
        // Zero-extend u16 → u32 using AVX2 (8 values at once)
        let u32x8 = _mm256_cvtepu16_epi32(u16x8);
        // Convert i32 → f32
        let f32x8 = _mm256_cvtepi32_ps(u32x8);
        _mm256_storeu_ps(output.as_mut_ptr().add(i), f32x8);
        i += 8;
    }

    // Scalar remainder
    for idx in i..count {
        *output.get_unchecked_mut(idx) = *data.get_unchecked(idx) as f32;
    }
}
