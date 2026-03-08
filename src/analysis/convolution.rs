/// Separable Gaussian convolution with SIMD dispatch.
///
/// A 2D Gaussian is separable into two 1D passes: horizontal then vertical.
/// For kernel half-width K, this is O(W×H×2K) vs O(W×H×K²) for naive 2D.
///
/// SIMD dispatch: AVX2 (8 px FMA), SSE2 (4 px FMA), NEON (4 px vfma).

use rayon::prelude::*;

/// Generate a zero-sum 1D Gaussian kernel.
///
/// Returns `(kernel, kernel_energy_sq)` where kernel has length `2*radius+1`,
/// is zero-mean (DAOFIND-style), and `kernel_energy_sq = Σ k²`.
pub fn generate_1d_kernel(sigma: f32) -> (Vec<f32>, f32) {
    let radius = (2.0 * sigma).ceil() as usize;
    let ksize = 2 * radius + 1;
    let inv_2s2 = 1.0 / (2.0 * sigma * sigma);

    let mut kernel = Vec::with_capacity(ksize);
    let mut ksum = 0.0_f32;
    for i in 0..ksize {
        let d = i as f32 - radius as f32;
        let g = (-inv_2s2 * d * d).exp();
        kernel.push(g);
        ksum += g;
    }

    // Subtract mean to make zero-sum
    let kmean = ksum / ksize as f32;
    let mut energy_sq = 0.0_f32;
    for v in kernel.iter_mut() {
        *v -= kmean;
        energy_sq += *v * *v;
    }

    (kernel, energy_sq)
}

/// Separable convolution: horizontal pass, then vertical pass.
///
/// Output is written to `output` (must be `width * height` elements).
/// The effective kernel_energy_sq for the 2D separable kernel is
/// `h_energy * v_energy` (returned by caller from two `generate_1d_kernel` calls,
/// but since both are the same kernel, it's `energy_sq²`).
///
/// Pixels within `radius` of the border are left as 0.0.
pub fn separable_convolve(
    data: &[f32],
    width: usize,
    height: usize,
    kernel: &[f32],
    output: &mut [f32],
) {
    let radius = kernel.len() / 2;

    // Temporary buffer for horizontal pass result
    let mut hpass = vec![0.0_f32; width * height];

    // Horizontal pass (parallelize over rows)
    hpass
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            let src_row = &data[y * width..(y + 1) * width];
            convolve_row_1d(src_row, kernel, radius, row);
        });

    // Vertical pass (parallelize over rows)
    output
        .par_chunks_mut(width) // We'll use a different parallel strategy
        .enumerate()
        .for_each(|(y, out_row)| {
            if y < radius || y >= height - radius {
                return;
            }
            convolve_vertical_row(&hpass, width, height, kernel, radius, y, out_row);
        });
}

/// Convolve a single row with a 1D kernel. SIMD dispatch with scalar fallback.
#[inline]
fn convolve_row_1d(src: &[f32], kernel: &[f32], radius: usize, dst: &mut [f32]) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
            unsafe { convolve_row_avx2(src, kernel, radius, dst) };
            return;
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        unsafe { convolve_row_neon(src, kernel, radius, dst) };
        return;
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        let width = src.len();
        let ksize = kernel.len();
        for x in radius..(width - radius) {
            let mut sum = 0.0_f32;
            for k in 0..ksize {
                sum += src[x + k - radius] * kernel[k];
            }
            dst[x] = sum;
        }
    }
}

/// Convolve a single output row vertically by gathering from hpass columns.
#[inline]
fn convolve_vertical_row(
    hpass: &[f32],
    width: usize,
    _height: usize,
    kernel: &[f32],
    radius: usize,
    y: usize,
    dst: &mut [f32],
) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
            unsafe { convolve_vertical_row_avx2(hpass, width, kernel, radius, y, dst) };
            return;
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        unsafe { convolve_vertical_row_neon(hpass, width, kernel, radius, y, dst) };
        return;
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        let ksize = kernel.len();
        for x in radius..(width - radius) {
            let mut sum = 0.0_f32;
            for k in 0..ksize {
                let iy = y + k - radius;
                sum += hpass[iy * width + x] * kernel[k];
            }
            dst[x] = sum;
        }
    }
}

// ── AVX2 implementations ────────────────────────────────────────────────────

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn convolve_row_avx2(src: &[f32], kernel: &[f32], radius: usize, dst: &mut [f32]) {
    use std::arch::x86_64::*;

    let width = src.len();
    let ksize = kernel.len();
    let end = width - radius;
    let mut x = radius;

    // Process 8 pixels at a time
    while x + 8 <= end {
        let mut acc = _mm256_setzero_ps();
        for k in 0..ksize {
            let kval = _mm256_set1_ps(kernel[k]);
            let src_ptr = src.as_ptr().add(x + k - radius);
            let pixels = _mm256_loadu_ps(src_ptr);
            acc = _mm256_fmadd_ps(pixels, kval, acc);
        }
        _mm256_storeu_ps(dst.as_mut_ptr().add(x), acc);
        x += 8;
    }

    // Scalar tail
    while x < end {
        let mut sum = 0.0_f32;
        for k in 0..ksize {
            sum += src[x + k - radius] * kernel[k];
        }
        dst[x] = sum;
        x += 1;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn convolve_vertical_row_avx2(
    hpass: &[f32],
    width: usize,
    kernel: &[f32],
    radius: usize,
    y: usize,
    dst: &mut [f32],
) {
    use std::arch::x86_64::*;

    let ksize = kernel.len();
    let end = width - radius;
    let mut x = radius;

    while x + 8 <= end {
        let mut acc = _mm256_setzero_ps();
        for k in 0..ksize {
            let iy = y + k - radius;
            let kval = _mm256_set1_ps(kernel[k]);
            let src_ptr = hpass.as_ptr().add(iy * width + x);
            let pixels = _mm256_loadu_ps(src_ptr);
            acc = _mm256_fmadd_ps(pixels, kval, acc);
        }
        _mm256_storeu_ps(dst.as_mut_ptr().add(x), acc);
        x += 8;
    }

    while x < end {
        let mut sum = 0.0_f32;
        for k in 0..ksize {
            let iy = y + k - radius;
            sum += hpass[iy * width + x] * kernel[k];
        }
        dst[x] = sum;
        x += 1;
    }
}

// ── NEON implementations ────────────────────────────────────────────────────

#[cfg(target_arch = "aarch64")]
unsafe fn convolve_row_neon(src: &[f32], kernel: &[f32], radius: usize, dst: &mut [f32]) {
    use std::arch::aarch64::*;

    let width = src.len();
    let ksize = kernel.len();
    let end = width - radius;
    let mut x = radius;

    while x + 4 <= end {
        let mut acc = vdupq_n_f32(0.0);
        for k in 0..ksize {
            let kval = vdupq_n_f32(kernel[k]);
            let pixels = vld1q_f32(src.as_ptr().add(x + k - radius));
            acc = vfmaq_f32(acc, pixels, kval);
        }
        vst1q_f32(dst.as_mut_ptr().add(x), acc);
        x += 4;
    }

    while x < end {
        let mut sum = 0.0_f32;
        for k in 0..ksize {
            sum += src[x + k - radius] * kernel[k];
        }
        dst[x] = sum;
        x += 1;
    }
}

#[cfg(target_arch = "aarch64")]
unsafe fn convolve_vertical_row_neon(
    hpass: &[f32],
    width: usize,
    kernel: &[f32],
    radius: usize,
    y: usize,
    dst: &mut [f32],
) {
    use std::arch::aarch64::*;

    let ksize = kernel.len();
    let end = width - radius;
    let mut x = radius;

    while x + 4 <= end {
        let mut acc = vdupq_n_f32(0.0);
        for k in 0..ksize {
            let iy = y + k - radius;
            let kval = vdupq_n_f32(kernel[k]);
            let pixels = vld1q_f32(hpass.as_ptr().add(iy * width + x));
            acc = vfmaq_f32(acc, pixels, kval);
        }
        vst1q_f32(dst.as_mut_ptr().add(x), acc);
        x += 4;
    }

    while x < end {
        let mut sum = 0.0_f32;
        for k in 0..ksize {
            let iy = y + k - radius;
            sum += hpass[iy * width + x] * kernel[k];
        }
        dst[x] = sum;
        x += 1;
    }
}

/// Separable B3-spline smoothing (DC-preserving).
///
/// Kernel: [1/16, 1/4, 3/8, 1/4, 1/16], radius=2.
/// Uses reflected boundary for border pixels.
/// This is the first layer of the à trous wavelet transform used for
/// MRS (Multiscale Residual Spectrum) noise estimation.
pub fn b3_spline_smooth(
    data: &[f32],
    width: usize,
    height: usize,
    output: &mut [f32],
) {
    const K: [f32; 5] = [1.0 / 16.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 4.0, 1.0 / 16.0];

    let mut hpass = vec![0.0_f32; width * height];

    // Horizontal pass with reflected boundary
    hpass
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            let src = &data[y * width..(y + 1) * width];
            for x in 0..width {
                let mut sum = 0.0_f32;
                for (ki, &kv) in K.iter().enumerate() {
                    let sx = x as i32 + ki as i32 - 2;
                    // Reflect at boundaries
                    let sx = if sx < 0 {
                        -sx
                    } else if sx >= width as i32 {
                        2 * width as i32 - 2 - sx
                    } else {
                        sx
                    } as usize;
                    sum += src[sx] * kv;
                }
                row[x] = sum;
            }
        });

    // Vertical pass with reflected boundary
    output
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            for x in 0..width {
                let mut sum = 0.0_f32;
                for (ki, &kv) in K.iter().enumerate() {
                    let sy = y as i32 + ki as i32 - 2;
                    let sy = if sy < 0 {
                        -sy
                    } else if sy >= height as i32 {
                        2 * height as i32 - 2 - sy
                    } else {
                        sy
                    } as usize;
                    sum += hpass[sy * width + x] * kv;
                }
                row[x] = sum;
            }
        });
}

/// Dilated B3-spline smoothing for à trous wavelet layers 2+.
///
/// Same kernel [1/16, 1/4, 3/8, 1/4, 1/16] but with spacing `2^(layer-1)`
/// between taps. Layer 1 has spacing 1 (use `b3_spline_smooth`).
/// Layer 2 has spacing 2, layer 3 has spacing 4, etc.
///
/// # Panics
///
/// Panics if `layer` is 0.
pub fn b3_spline_smooth_dilated(
    data: &[f32],
    width: usize,
    height: usize,
    output: &mut [f32],
    layer: usize,
) {
    use rayon::prelude::*;

    const K: [f32; 5] = [1.0 / 16.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 4.0, 1.0 / 16.0];
    let spacing = 1_i32 << (layer - 1); // 2^(layer-1)

    let mut hpass = vec![0.0_f32; width * height];

    // Horizontal pass with reflected boundary
    hpass
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            let src = &data[y * width..(y + 1) * width];
            for x in 0..width {
                let mut sum = 0.0_f32;
                for (ki, &kv) in K.iter().enumerate() {
                    let sx = x as i32 + (ki as i32 - 2) * spacing;
                    let sx = if sx < 0 {
                        -sx
                    } else if sx >= width as i32 {
                        2 * width as i32 - 2 - sx
                    } else {
                        sx
                    } as usize;
                    sum += src[sx.min(width - 1)] * kv;
                }
                row[x] = sum;
            }
        });

    // Vertical pass with reflected boundary
    output
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            for x in 0..width {
                let mut sum = 0.0_f32;
                for (ki, &kv) in K.iter().enumerate() {
                    let sy = y as i32 + (ki as i32 - 2) * spacing;
                    let sy = if sy < 0 {
                        -sy
                    } else if sy >= height as i32 {
                        2 * height as i32 - 2 - sy
                    } else {
                        sy
                    } as usize;
                    sum += hpass[sy.min(height - 1) * width + x] * kv;
                }
                row[x] = sum;
            }
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_1d_kernel_zero_sum() {
        let (kernel, energy) = generate_1d_kernel(1.5);
        let sum: f32 = kernel.iter().sum();
        assert!(sum.abs() < 1e-5, "Kernel sum {} should be ~0", sum);
        assert!(energy > 0.0, "Energy should be positive");
    }

    #[test]
    fn test_separable_detects_point_source() {
        // Verify separable convolution finds a point source at the correct location
        // and suppresses flat background (zero-sum property).
        let width = 100;
        let height = 100;
        let sigma = 1.5_f32;

        // Point source on flat background
        let mut data = vec![1000.0_f32; width * height];
        data[50 * width + 50] += 5000.0;

        let (kernel_1d, _) = generate_1d_kernel(sigma);
        let radius = kernel_1d.len() / 2;
        let mut output = vec![0.0_f32; width * height];
        separable_convolve(&data, width, height, &kernel_1d, &mut output);

        // Peak should be at (50, 50)
        let peak_val = output.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let peak_idx = output.iter().position(|&v| v == peak_val).unwrap();
        assert_eq!(peak_idx, 50 * width + 50, "Peak should be at (50, 50)");
        assert!(peak_val > 100.0, "Peak response should be strongly positive");

        // Away from the point source, response should be near zero (flat bg suppressed)
        let far_val = output[20 * width + 20];
        assert!(
            far_val.abs() < 1.0,
            "Far from source, response {} should be ~0",
            far_val
        );
    }

    #[test]
    fn test_b3_spline_flat_image_preserved() {
        // A flat image should be unchanged after B3 spline smoothing (DC-preserving)
        let width = 50;
        let height = 50;
        let val = 1000.0_f32;
        let data = vec![val; width * height];
        let mut output = vec![0.0_f32; width * height];
        b3_spline_smooth(&data, width, height, &mut output);

        for y in 0..height {
            for x in 0..width {
                assert!(
                    (output[y * width + x] - val).abs() < 0.01,
                    "Flat image at ({},{}) = {} should be ~{}",
                    x, y, output[y * width + x], val
                );
            }
        }
    }

    #[test]
    fn test_b3_spline_smooths_noise() {
        // B3 spline should reduce high-frequency noise while preserving large-scale structure
        let width = 100;
        let height = 100;
        let mut data = vec![1000.0_f32; width * height];

        // Add a point source
        data[50 * width + 50] += 5000.0;

        let mut output = vec![0.0_f32; width * height];
        b3_spline_smooth(&data, width, height, &mut output);

        // The point source should be spread out (peak lower than original)
        assert!(
            output[50 * width + 50] < data[50 * width + 50],
            "Peak should be smoothed: {} < {}",
            output[50 * width + 50], data[50 * width + 50]
        );
        // But still present
        assert!(
            output[50 * width + 50] > 1000.0,
            "Peak should still be above background"
        );
    }

    #[test]
    fn test_separable_flat_image_zero_output() {
        // A flat image convolved with a zero-sum kernel should produce ~0
        let width = 50;
        let height = 50;
        let data = vec![1000.0_f32; width * height];
        let (kernel, _) = generate_1d_kernel(1.5);
        let mut output = vec![0.0_f32; width * height];
        separable_convolve(&data, width, height, &kernel, &mut output);

        let radius = kernel.len() / 2;
        for y in radius..(height - radius) {
            for x in radius..(width - radius) {
                assert!(
                    output[y * width + x].abs() < 0.1,
                    "Flat image at ({},{}) = {} should be ~0",
                    x, y, output[y * width + x]
                );
            }
        }
    }

    #[test]
    fn test_b3_spline_dilated_layer2() {
        // Both layers produce the same center-peak for a point source (the center
        // tap always captures the spike), but layer 2 places its outer taps at
        // spacing 2, so the pixel immediately adjacent to the spike receives more
        // energy under layer 1 (tap at distance 1) than under layer 2 (nearest
        // non-center tap is at distance 2, skipping the adjacent pixel entirely).
        let w = 64;
        let h = 64;
        let mut data = vec![0.0_f32; w * h];
        data[32 * w + 32] = 1000.0;

        let mut out1 = vec![0.0_f32; w * h];
        b3_spline_smooth(&data, w, h, &mut out1);
        // Layer 1: pixel at distance 1 from the spike gets weight 1/4 * 3/8 from
        // the separable horizontal+vertical passes.
        let near1_h = out1[32 * w + 33]; // one pixel right of center (horizontal neighbour)

        let mut out2 = vec![0.0_f32; w * h];
        b3_spline_smooth_dilated(&data, w, h, &mut out2, 2);
        // Layer 2 spacing=2: the nearest non-center taps land at distance 2, so
        // the pixel at distance 1 from the spike receives zero from the horizontal
        // pass and thus zero overall.
        let near2_h = out2[32 * w + 33];

        assert!(
            near1_h > near2_h,
            "Layer 1 adjacent response {} should exceed layer 2 adjacent response {} \
             (layer 2 skips distance-1 pixels)",
            near1_h,
            near2_h
        );

        // Conversely, layer 2 places energy at distance 2 that layer 1 also reaches,
        // but layer 2 carries more weight there because spacing=2 means the distance-2
        // pixel is the *nearest* neighbour tap (weight 1/4), not the second one.
        let far1 = out1[32 * w + 34]; // distance 2 from centre
        let far2 = out2[32 * w + 34];
        assert!(
            far2 > far1,
            "Layer 2 response at distance 2 ({}) should exceed layer 1 ({})",
            far2,
            far1
        );
    }

    #[test]
    fn test_b3_spline_dilated_flat_preserved() {
        let w = 32;
        let h = 32;
        let data = vec![500.0_f32; w * h];
        let mut out = vec![0.0_f32; w * h];
        b3_spline_smooth_dilated(&data, w, h, &mut out, 2);
        for &v in &out {
            assert!(
                (v - 500.0).abs() < 0.01,
                "Flat image should be preserved, got {}",
                v
            );
        }
    }
}
