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
/// Pure-Rust encoder with NEON (aarch64) and AVX2/SSE2 (x86_64) code paths and a
/// scalar fallback — no C toolchain, no cmake/nasm, so cross-compilation and
/// distro packaging need nothing beyond the Rust toolchain.
///
/// `channels` must be 3 (RGB) or 4 (RGBA — the alpha byte is read but ignored
/// by the encoder, so callers must not pre-strip it).
///
/// `quality` is clamped to 1..=100, matching `ImageConverter::with_quality`.
/// The encoder clamps internally too; doing it here keeps the documented range
/// a property of this function rather than of whichever backend is underneath.
/// Note this differs from the old `turbojpeg` binding, which returned an error
/// for an out-of-range quality instead of clamping.
pub fn encode_jpeg(
    data: &[u8],
    width: usize,
    height: usize,
    channels: usize,
    quality: u8,
) -> Result<Vec<u8>> {
    let quality = quality.clamp(1, 100);
    let format = match channels {
        3 => libjpeg_turbo_rs::PixelFormat::Rgb,
        4 => libjpeg_turbo_rs::PixelFormat::Rgba,
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

    libjpeg_turbo_rs::compress(
        data,
        width,
        height,
        format,
        quality,
        libjpeg_turbo_rs::Subsampling::S420,
    )
    .map_err(|e| anyhow::anyhow!("JPEG encoding failed: {e}"))
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
            // RGBA is handed to the encoder as-is — it reads 4 bytes per pixel
            // and discards alpha itself, so the old de-interleaving Vec
            // allocation (W*H*3 bytes per frame) is gone. Byte-equality of the
            // two paths is pinned by
            // `tests/jpeg_encoder.rs::rgba_input_encodes_identically_to_hand_stripped_rgb`.
            let jpeg_data = encode_jpeg(
                &image.data,
                image.width,
                image.height,
                image.channels as usize,
                quality,
            )?;

            let file = File::create(path).context("Failed to create output file")?;
            let mut writer = BufWriter::new(file);
            writer
                .write_all(&jpeg_data)
                .context("Failed to write JPEG data")?;
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

    /// Out-of-range quality clamps to 1..=100 instead of erroring (the old
    /// `turbojpeg` binding rejected it) or reaching the encoder unbounded.
    #[test]
    fn encode_jpeg_clamps_out_of_range_quality() {
        let rgb = synth_rgb(32, 24);
        assert_eq!(
            encode_jpeg(&rgb, 32, 24, 3, 0).unwrap(),
            encode_jpeg(&rgb, 32, 24, 3, 1).unwrap(),
            "quality 0 should behave as 1"
        );
        assert_eq!(
            encode_jpeg(&rgb, 32, 24, 3, 200).unwrap(),
            encode_jpeg(&rgb, 32, 24, 3, 100).unwrap(),
            "quality above 100 should behave as 100"
        );
    }

    #[test]
    fn encode_jpeg_rejects_bad_channel_count_and_short_buffer() {
        let rgb = synth_rgb(8, 8);
        assert!(encode_jpeg(&rgb, 8, 8, 2, 90).is_err());
        assert!(encode_jpeg(&rgb[..10], 8, 8, 3, 90).is_err());
    }
}
