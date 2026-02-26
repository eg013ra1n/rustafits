use std::path::Path;

use anyhow::Result;
use rayon::prelude::*;

use crate::formats;
use crate::processing::{binning, color, debayer, downscale, stretch};
use crate::types::{BayerPattern, ImageMetadata, PixelData, ProcessConfig, ProcessedImage};

pub fn process_image(path: &Path, config: &ProcessConfig) -> Result<ProcessedImage> {
    let (meta, pixels) = formats::read_image(path)?;

    match pixels {
        PixelData::Uint16(data) => process_u16(data, meta, config),
        PixelData::Float32(data) => process_f32(data, meta, config),
    }
}

fn process_u16(
    mut data: Vec<u16>,
    meta: ImageMetadata,
    config: &ProcessConfig,
) -> Result<ProcessedImage> {
    let mut width = meta.width;
    let mut height = meta.height;

    let (float_data, is_color, num_channels);

    if config.apply_debayer && meta.bayer_pattern != BayerPattern::None {
        // Debayer first: u16 mono → f32 planar RGB, half dims (inherent 2x reduction)
        let (mut rgb, mut ow, mut oh) =
            debayer::super_pixel_debayer_u16(&data, width, height, meta.bayer_pattern);
        // Debayer counts as 2x, so only apply additional downscale for factor > 2
        let extra = config.downscale_factor / 2;
        if extra > 1 {
            let (d, nw, nh) = downscale::downscale_f32_planar(&rgb, ow, oh, 3, extra);
            rgb = d;
            ow = nw;
            oh = nh;
        }
        width = ow;
        height = oh;
        float_data = rgb;
        is_color = true;
        num_channels = 3;
    } else {
        // Non-Bayer: downscale raw data directly
        if config.downscale_factor > 1 {
            let (d, nw, nh) =
                downscale::downscale_u16(&data, width, height, config.downscale_factor);
            data = d;
            width = nw;
            height = nh;
        }

        // Convert to float
        let mut fdata = color::u16_to_f32(&data);

        // Preview binning for mono
        if config.preview_mode && meta.bayer_pattern == BayerPattern::None {
            let (binned, nw, nh) = binning::bin_2x2_float(&fdata, width, height);
            fdata = binned;
            width = nw;
            height = nh;
        }

        float_data = fdata;
        is_color = false;
        num_channels = 1;
    }

    apply_stretch_and_finalize(float_data, width, height, is_color, num_channels, &meta, config)
}

fn process_f32(
    mut data: Vec<f32>,
    meta: ImageMetadata,
    config: &ProcessConfig,
) -> Result<ProcessedImage> {
    let mut width = meta.width;
    let mut height = meta.height;

    let (float_data, is_color, num_channels);

    if meta.channels == 1 && config.apply_debayer && meta.bayer_pattern != BayerPattern::None {
        // Debayer first: f32 mono → f32 planar RGB, half dims (inherent 2x reduction)
        let (mut rgb, mut ow, mut oh) =
            debayer::super_pixel_debayer_f32(&data, width, height, meta.bayer_pattern);
        // Debayer counts as 2x, so only apply additional downscale for factor > 2
        let extra = config.downscale_factor / 2;
        if extra > 1 {
            let (d, nw, nh) = downscale::downscale_f32_planar(&rgb, ow, oh, 3, extra);
            rgb = d;
            ow = nw;
            oh = nh;
        }
        width = ow;
        height = oh;
        float_data = rgb;
        is_color = true;
        num_channels = 3;
    } else if meta.channels == 3 {
        // Already RGB: downscale planar data directly
        if config.downscale_factor > 1 {
            let (d, nw, nh) =
                downscale::downscale_f32_planar(&data, width, height, 3, config.downscale_factor);
            data = d;
            width = nw;
            height = nh;
        }
        float_data = data;
        is_color = true;
        num_channels = 3;
    } else {
        // Mono f32 (no Bayer): downscale then optional binning
        if config.downscale_factor > 1 {
            let (d, nw, nh) =
                downscale::downscale_f32_planar(&data, width, height, 1, config.downscale_factor);
            data = d;
            width = nw;
            height = nh;
        }
        if config.preview_mode && meta.channels == 1 && meta.bayer_pattern == BayerPattern::None {
            let (binned, nw, nh) = binning::bin_2x2_float(&data, width, height);
            data = binned;
            width = nw;
            height = nh;
        }
        float_data = data;
        is_color = false;
        num_channels = 1;
    }

    apply_stretch_and_finalize(float_data, width, height, is_color, num_channels, &meta, config)
}

fn compute_stretch_coefficients(channel_data: &[f32]) -> (f32, f32, f32, f32, f32) {
    let max_input = 65536.0f32;
    let params = stretch::compute_stretch_params(channel_data, max_input);
    let hs_range_factor = if params.highlights == params.shadows {
        1.0
    } else {
        1.0 / (params.highlights - params.shadows)
    };
    let native_shadows = params.shadows * max_input;
    let native_highlights = params.highlights * max_input;
    let k1 = (params.midtones - 1.0) * hs_range_factor * 255.0 / max_input;
    let k2 = (2.0 * params.midtones - 1.0) * hs_range_factor / max_input;
    (native_shadows, native_highlights, k1, k2, params.midtones)
}

fn apply_stretch_and_finalize(
    float_data: Vec<f32>,
    width: usize,
    height: usize,
    is_color: bool,
    num_channels: usize,
    meta: &ImageMetadata,
    config: &ProcessConfig,
) -> Result<ProcessedImage> {
    let channel_size = width * height;
    let mut rgb_data = vec![0u8; width * height * 3];

    if config.auto_stretch {
        if is_color {
            // Compute stretch params for each channel in parallel (includes quickselect)
            let coeffs: Vec<_> = (0..num_channels)
                .into_par_iter()
                .map(|c| {
                    let ch = &float_data[c * channel_size..(c + 1) * channel_size];
                    compute_stretch_coefficients(ch)
                })
                .collect();

            // Apply stretch with stride=3 directly into rgb_data (no intermediate buffers)
            for c in 0..num_channels {
                let ch = &float_data[c * channel_size..(c + 1) * channel_size];
                let (ns, nh, k1, k2, m) = coeffs[c];
                stretch::apply_stretch(ch, &mut rgb_data, c, 3, ns, nh, k1, k2, m);
            }
        } else {
            let channel_data = &float_data[0..channel_size];
            let (native_shadows, native_highlights, k1, k2, midtones) =
                compute_stretch_coefficients(channel_data);

            let mut temp = vec![0u8; channel_size];
            stretch::apply_stretch(
                channel_data,
                &mut temp,
                0,
                1,
                native_shadows,
                native_highlights,
                k1,
                k2,
                midtones,
            );
            rgb_data = color::replicate_gray_to_rgb(&temp);
        }
    }

    // Vertical flip
    if meta.flip_vertical {
        color::vertical_flip_rgb(&mut rgb_data, width, height);
    }

    Ok(ProcessedImage {
        data: rgb_data,
        width,
        height,
        is_color,
    })
}
