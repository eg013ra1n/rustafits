pub fn downscale_u16(input: &[u16], width: usize, height: usize, factor: usize) -> (Vec<u16>, usize, usize) {
    let new_w = width / factor;
    let new_h = height / factor;
    let mut output = vec![0u16; new_w * new_h];

    for y in 0..new_h {
        for x in 0..new_w {
            output[y * new_w + x] = input[(y * factor) * width + (x * factor)];
        }
    }

    (output, new_w, new_h)
}

#[allow(dead_code)]
pub fn downscale_f32(input: &[f32], width: usize, height: usize, factor: usize) -> (Vec<f32>, usize, usize) {
    let new_w = width / factor;
    let new_h = height / factor;
    let mut output = vec![0f32; new_w * new_h];

    for y in 0..new_h {
        for x in 0..new_w {
            output[y * new_w + x] = input[(y * factor) * width + (x * factor)];
        }
    }

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

    for c in 0..channels {
        let src = &input[c * plane_in..];
        let dst = &mut output[c * plane_out..];
        for y in 0..new_h {
            for x in 0..new_w {
                dst[y * new_w + x] = src[(y * factor) * width + (x * factor)];
            }
        }
    }

    (output, new_w, new_h)
}
