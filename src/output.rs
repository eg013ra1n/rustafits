use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use image::codecs::png::PngEncoder;
use image::{ColorType, ImageEncoder};

use crate::types::ProcessedImage;

/// Encode an RGB/RGBA buffer to a baseline JPEG byte stream (4:2:0 chroma
/// subsampling), in memory.
///
/// `channels` must be 3 (RGB) or 4 (RGBA — the alpha byte is read but ignored
/// by the encoder, so callers must not pre-strip it).
pub fn encode_jpeg(
    data: &[u8],
    width: usize,
    height: usize,
    channels: usize,
    quality: u8,
) -> Result<Vec<u8>> {
    let format = match channels {
        3 => turbojpeg::PixelFormat::RGB,
        4 => turbojpeg::PixelFormat::RGBA,
        n => anyhow::bail!("unsupported channel count for JPEG encoding: {n}"),
    };

    let expected = width
        .checked_mul(height)
        .and_then(|n| n.checked_mul(channels))
        .context("image dimensions overflow")?;
    if data.len() < expected {
        anyhow::bail!(
            "pixel buffer too small: {} bytes, need {expected} for {width}x{height}x{channels}",
            data.len()
        );
    }

    let image = turbojpeg::Image {
        pixels: data,
        width,
        pitch: width * channels,
        height,
        format,
    };
    let jpeg = turbojpeg::compress(image, quality as i32, turbojpeg::Subsamp::Sub2x2)
        .context("JPEG encoding failed")?;
    Ok(jpeg.to_vec())
}

pub fn save_image(image: &ProcessedImage, path: &Path, quality: u8) -> Result<()> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    match ext.as_str() {
        "png" => {
            let file = File::create(path).context("Failed to create output file")?;
            let writer = BufWriter::new(file);
            let color_type = if image.channels == 4 {
                ColorType::Rgba8
            } else {
                ColorType::Rgb8
            };
            let encoder = PngEncoder::new(writer);
            encoder
                .write_image(
                    &image.data,
                    image.width as u32,
                    image.height as u32,
                    color_type.into(),
                )
                .context("PNG encoding failed")?;
        }
        _ => {
            // JPEG via libjpeg-turbo (SIMD-accelerated DCT on NEON/AVX2)
            // Strip alpha if RGBA; otherwise pass RGB data by reference (zero-copy).
            let rgb_tmp;
            let pixels: &[u8] = if image.channels == 4 {
                rgb_tmp = image
                    .data
                    .chunks_exact(4)
                    .flat_map(|px| [px[0], px[1], px[2]])
                    .collect::<Vec<u8>>();
                &rgb_tmp
            } else {
                &image.data
            };

            let tj_image = turbojpeg::Image {
                pixels,
                width: image.width,
                pitch: image.width * 3,
                height: image.height,
                format: turbojpeg::PixelFormat::RGB,
            };

            let jpeg_data = turbojpeg::compress(tj_image, quality as i32, turbojpeg::Subsamp::Sub2x2)
                .context("JPEG encoding failed")?;

            let file = File::create(path).context("Failed to create output file")?;
            let mut writer = BufWriter::new(file);
            writer.write_all(&jpeg_data).context("Failed to write JPEG data")?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synth_rgb(w: usize, h: usize) -> Vec<u8> {
        (0..w * h * 3).map(|i| (i % 251) as u8).collect()
    }

    #[test]
    fn encode_jpeg_rgb_produces_a_jpeg_stream() {
        let jpeg = encode_jpeg(&synth_rgb(32, 24), 32, 24, 3, 90).unwrap();
        assert!(jpeg.len() > 2);
        assert_eq!(&jpeg[..2], &[0xFF, 0xD8], "SOI marker");
        assert_eq!(&jpeg[jpeg.len() - 2..], &[0xFF, 0xD9], "EOI marker");
    }

    /// The RGBA contract: alpha is read but ignored, so an RGBA buffer with the
    /// same colors as an RGB one must produce byte-identical output.
    #[test]
    fn encode_jpeg_rgba_ignores_alpha() {
        let rgb = synth_rgb(32, 24);
        let rgba: Vec<u8> = rgb
            .chunks_exact(3)
            .flat_map(|px| [px[0], px[1], px[2], 0xAB])
            .collect();
        let from_rgb = encode_jpeg(&rgb, 32, 24, 3, 90).unwrap();
        let from_rgba = encode_jpeg(&rgba, 32, 24, 4, 90).unwrap();
        assert_eq!(from_rgb, from_rgba);
    }

    #[test]
    fn encode_jpeg_rejects_bad_channel_count_and_short_buffer() {
        let rgb = synth_rgb(8, 8);
        assert!(encode_jpeg(&rgb, 8, 8, 2, 90).is_err());
        assert!(encode_jpeg(&rgb[..10], 8, 8, 3, 90).is_err());
    }
}
