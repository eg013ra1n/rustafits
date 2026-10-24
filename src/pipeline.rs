use std::path::Path;

use anyhow::Result;
use rayon::prelude::*;

use crate::formats;
use crate::processing::{binning, color, debayer, downscale, stretch, vng};
use crate::types::{BayerPattern, ImageMetadata, PixelData, ProcessConfig, ProcessedImage};

pub fn process_image(path: &Path, config: &ProcessConfig) -> Result<ProcessedImage> {
    let (meta, pixels) = formats::read_image(path)?;
    process_image_data(meta, pixels, config)
}

pub fn process_image_data(
    meta: ImageMetadata,
    pixels: PixelData,
    config: &ProcessConfig,
) -> Result<ProcessedImage> {
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
        // Debayer first: u16 mono → f32 planar RGB.
        //
        // The two debayers reduce differently, and that is what `extra` below
        // is about. Super-pixel folds each 2x2 Bayer tile into one pixel, so it
        // *is* a 2x downscale and only the remainder of `downscale_factor` is
        // left to apply. VNG keeps the native grid, so the whole factor still
        // has to be applied.
        let (mut rgb, mut ow, mut oh) = if config.vng_debayer {
            let cfa = color::u16_to_f32(&data);
            // The u16 plane is dead from here on; release it before VNG
            // allocates its 12-bytes-per-pixel planar RGB output.
            drop(std::mem::take(&mut data));
            let rgb = vng::vng_debayer_f32(&cfa, width, height, meta.bayer_pattern);
            (rgb, width, height)
        } else {
            debayer::super_pixel_debayer_u16(&data, width, height, meta.bayer_pattern)
        };
        let extra = if config.vng_debayer {
            config.downscale_factor
        } else {
            config.downscale_factor / 2
        };
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

        // Preview binning for mono (needs f32 for binning)
        if config.preview_mode && meta.bayer_pattern == BayerPattern::None {
            let mut fdata = color::u16_to_f32(&data);
            let (binned, nw, nh) = binning::bin_2x2_float(&fdata, width, height);
            fdata = binned;
            width = nw;
            height = nh;
            float_data = fdata;
            is_color = false;
            num_channels = 1;
            return apply_stretch_and_finalize(float_data, width, height, is_color, num_channels, &meta, config);
        }

        // Fast path: stretch directly from u16, skip f32 intermediate
        return apply_stretch_u16_and_finalize(&data, width, height, &meta, config);
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
        // Debayer first: f32 mono → f32 planar RGB. See `process_u16` for why
        // `extra` differs between the two debayers.
        let (mut rgb, mut ow, mut oh) = if config.vng_debayer {
            let rgb = vng::vng_debayer_f32(&data, width, height, meta.bayer_pattern);
            // The CFA plane is dead from here on — see `process_u16`.
            drop(std::mem::take(&mut data));
            (rgb, width, height)
        } else {
            debayer::super_pixel_debayer_f32(&data, width, height, meta.bayer_pattern)
        };
        let extra = if config.vng_debayer {
            config.downscale_factor
        } else {
            config.downscale_factor / 2
        };
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

fn compute_stretch_coefficients_u16(channel_data: &[u16]) -> stretch::StretchCoeffs {
    let max_input = 65536.0f32;
    let params = stretch::compute_stretch_params_u16(channel_data, max_input);
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

/// Fast u16→u8 stretch for mono images, skipping the f32 intermediate.
fn apply_stretch_u16_and_finalize(
    data: &[u16],
    width: usize,
    height: usize,
    meta: &ImageMetadata,
    config: &ProcessConfig,
) -> Result<ProcessedImage> {
    let channel_size = width * height;
    let bpp: usize = if config.rgba_output { 4 } else { 3 };

    let mut out_data = vec![0u8; channel_size * bpp];
    if config.auto_stretch {
        let coeffs = compute_stretch_coefficients_u16(&data[..channel_size]);
        // Stretch and gray→RGB(A) in one parallel pass, no full-frame gray buffer.
        stretch::stretch_gray_u16_interleaved(&data[..channel_size], &mut out_data, bpp, coeffs);
    }

    let mut result = ProcessedImage {
        data: out_data,
        width,
        height,
        is_color: false,
        channels: bpp as u8,
        flip_vertical: meta.flip_vertical,
    };

    if meta.flip_vertical {
        color::vertical_flip(&mut result.data, width, height, bpp);
    }

    Ok(result)
}

fn compute_stretch_coefficients(channel_data: &[f32]) -> stretch::StretchCoeffs {
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
    let bpp: usize = if config.rgba_output { 4 } else { 3 };
    let mut out_data = vec![0u8; channel_size * bpp];

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

            // All channels (and alpha) interleaved in one parallel pass.
            stretch::stretch_planar_interleaved(&float_data, channel_size, &mut out_data, bpp, &coeffs);
        } else {
            let channel_data = &float_data[0..channel_size];
            let coeffs = compute_stretch_coefficients(channel_data);
            stretch::stretch_gray_f32_interleaved(channel_data, &mut out_data, bpp, coeffs);
        }
    }

    // Vertical flip
    if meta.flip_vertical {
        color::vertical_flip(&mut out_data, width, height, bpp);
    }

    Ok(ProcessedImage {
        data: out_data,
        width,
        height,
        is_color,
        channels: bpp as u8,
        flip_vertical: meta.flip_vertical,
    })
}

#[cfg(test)]
mod tests {
    use crate::types::{DataType, ImageMetadata, PixelData};
    use crate::{BayerPattern, ImageConverter};

    fn cfa_meta(w: usize, h: usize, dtype: DataType) -> ImageMetadata {
        ImageMetadata {
            width: w,
            height: h,
            channels: 1,
            dtype,
            bayer_pattern: BayerPattern::Rggb,
            flip_vertical: false,
        }
    }

    /// A CFA ramp: enough structure that the gradient method has something to
    /// follow, and never uniform (a flat frame stretches to nothing).
    fn ramp_u16(w: usize, h: usize) -> Vec<u16> {
        (0..w * h).map(|i| ((i * 37) % 4096) as u16).collect()
    }

    #[test]
    fn vng_keeps_native_dimensions() {
        let (w, h) = (16, 16);
        let out = ImageConverter::new()
            .with_vng_debayer()
            .process_data(
                cfa_meta(w, h, DataType::Uint16),
                PixelData::Uint16(ramp_u16(w, h)),
            )
            .unwrap();

        assert_eq!((out.width, out.height), (w, h));
        assert!(out.is_color);
        assert_eq!(out.channels, 3);
        assert_eq!(out.data.len(), w * h * 3);
    }

    #[test]
    fn superpixel_still_halves_without_the_flag() {
        let (w, h) = (16, 16);
        let out = ImageConverter::new()
            .process_data(
                cfa_meta(w, h, DataType::Uint16),
                PixelData::Uint16(ramp_u16(w, h)),
            )
            .unwrap();

        assert_eq!((out.width, out.height), (w / 2, h / 2));
        assert!(out.is_color);
    }

    /// Super-pixel debayer *is* a 2x reduction, so it consumes half of
    /// `downscale_factor`. VNG is not, so the whole factor must still apply —
    /// inheriting the `/2` compensation would silently halve the request.
    #[test]
    fn vng_applies_the_full_downscale_factor() {
        let (w, h) = (16, 16);
        let vng = ImageConverter::new()
            .with_vng_debayer()
            .with_downscale(2)
            .process_data(
                cfa_meta(w, h, DataType::Uint16),
                PixelData::Uint16(ramp_u16(w, h)),
            )
            .unwrap();
        assert_eq!((vng.width, vng.height), (w / 2, h / 2));

        let superpixel = ImageConverter::new()
            .with_downscale(2)
            .process_data(
                cfa_meta(w, h, DataType::Uint16),
                PixelData::Uint16(ramp_u16(w, h)),
            )
            .unwrap();
        assert_eq!((superpixel.width, superpixel.height), (w / 2, h / 2));
    }

    #[test]
    fn vng_covers_the_f32_input_path() {
        let (w, h) = (16, 16);
        let data: Vec<f32> = ramp_u16(w, h).into_iter().map(|v| v as f32).collect();
        let out = ImageConverter::new()
            .with_vng_debayer()
            .process_data(
                cfa_meta(w, h, DataType::Float32),
                PixelData::Float32(data),
            )
            .unwrap();

        assert_eq!((out.width, out.height), (w, h));
        assert!(out.is_color);
        assert_eq!(out.channels, 3);
    }

    /// The gradient method needs a full 5x5 window; a frame smaller than that
    /// has no interior at all and must fall through to the border pass rather
    /// than panic.
    #[test]
    fn vng_survives_a_frame_smaller_than_its_window() {
        let (w, h) = (4, 4);
        let out = ImageConverter::new()
            .with_vng_debayer()
            .process_data(
                cfa_meta(w, h, DataType::Uint16),
                PixelData::Uint16(ramp_u16(w, h)),
            )
            .unwrap();

        assert_eq!((out.width, out.height), (w, h));
    }
}
