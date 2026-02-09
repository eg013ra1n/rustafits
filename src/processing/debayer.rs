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
