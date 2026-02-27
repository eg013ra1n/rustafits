use rayon::prelude::*;

pub fn downscale_u16(input: &[u16], width: usize, height: usize, factor: usize) -> (Vec<u16>, usize, usize) {
    let new_w = width / factor;
    let new_h = height / factor;
    let mut output = vec![0u16; new_w * new_h];

    output
        .par_chunks_mut(new_w)
        .enumerate()
        .for_each(|(y, row)| {
            let src_y = y * factor;
            for (x, dst) in row.iter_mut().enumerate() {
                *dst = input[src_y * width + x * factor];
            }
        });

    (output, new_w, new_h)
}

#[allow(dead_code)]
pub fn downscale_f32(input: &[f32], width: usize, height: usize, factor: usize) -> (Vec<f32>, usize, usize) {
    let new_w = width / factor;
    let new_h = height / factor;
    let mut output = vec![0f32; new_w * new_h];

    output
        .par_chunks_mut(new_w)
        .enumerate()
        .for_each(|(y, row)| {
            let src_y = y * factor;
            for (x, dst) in row.iter_mut().enumerate() {
                *dst = input[src_y * width + x * factor];
            }
        });

    (output, new_w, new_h)
}

/// Downscale planar multi-channel f32 data (each channel stored contiguously).
pub fn downscale_f32_planar(
    input: &[f32],
    width: usize,
    height: usize,
    channels: usize,
    factor: usize,
) -> (Vec<f32>, usize, usize) {
    let new_w = width / factor;
    let new_h = height / factor;
    let plane_in = width * height;
    let plane_out = new_w * new_h;
    let mut output = vec![0f32; plane_out * channels];

    output
        .par_chunks_mut(plane_out)
        .enumerate()
        .for_each(|(c, dst_plane)| {
            let src = &input[c * plane_in..];
            dst_plane
                .par_chunks_mut(new_w)
                .enumerate()
                .for_each(|(y, row)| {
                    let src_y = y * factor;
                    let src_row = &src[src_y * width..];
                    if factor == 2 {
                        downscale_row_2x(src_row, row);
                    } else {
                        for (x, dst) in row.iter_mut().enumerate() {
                            *dst = src_row[x * factor];
                        }
                    }
                });
        });

    (output, new_w, new_h)
}

/// SIMD-accelerated factor-2 nearest-neighbor downscale for a single row.
fn downscale_row_2x(src: &[f32], dst: &mut [f32]) {
    let n = dst.len();
    let mut x = 0;

    #[cfg(target_arch = "aarch64")]
    {
        use std::arch::aarch64::*;
        // vld2q_f32 deinterleaves: .0 = even indices, .1 = odd indices
        while x + 4 <= n {
            unsafe {
                let pair = vld2q_f32(src.as_ptr().add(x * 2));
                vst1q_f32(dst.as_mut_ptr().add(x), pair.0);
            }
            x += 4;
        }
    }

    #[cfg(target_arch = "x86_64")]
    {
        use std::arch::x86_64::*;
        // SSE: load 8 f32, shuffle to extract even-indexed elements (4 outputs)
        while x + 4 <= n {
            unsafe {
                let a = _mm_loadu_ps(src.as_ptr().add(x * 2));     // [0,1,2,3]
                let b = _mm_loadu_ps(src.as_ptr().add(x * 2 + 4)); // [4,5,6,7]
                let evens = _mm_shuffle_ps(a, b, 0b10_00_10_00);   // [0,2,4,6]
                _mm_storeu_ps(dst.as_mut_ptr().add(x), evens);
            }
            x += 4;
        }
    }

    // Scalar remainder
    for i in x..n {
        dst[i] = src[i * 2];
    }
}
