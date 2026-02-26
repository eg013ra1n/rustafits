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
                    for (x, dst) in row.iter_mut().enumerate() {
                        *dst = src[src_y * width + x * factor];
                    }
                });
        });

    (output, new_w, new_h)
}
