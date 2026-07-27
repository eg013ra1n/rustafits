use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use image::codecs::png::PngEncoder;
use image::{ColorType, ImageEncoder};

use crate::types::ProcessedImage;

/// Encode an RGB/RGBA buffer to a baseline JPEG byte stream (4:2:0 chroma
/// subsampling).
///
/// Pure-Rust encoder with NEON (aarch64) and AVX2/SSE2 (x86_64) code paths and a
/// scalar fallback — no C toolchain, no cmake/nasm, so cross-compilation and
/// distro packaging need nothing beyond the Rust toolchain.
///
/// `channels` must be 3 (RGB) or 4 (RGBA — the alpha byte is read but ignored,
/// so no de-interleaving copy is needed).
pub fn encode_jpeg(data: &[u8], width: usize, height: usize, channels: usize, quality: u8) -> Result<Vec<u8>> {
    let pixel_format = match channels {
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
        pixel_format,
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
            // RGBA is fed to the encoder as-is — it reads 4 bytes per pixel and
            // drops alpha itself, so the old de-interleaving Vec allocation
            // (W*H*3 bytes per frame) is gone.
            let jpeg_data = encode_jpeg(
                &image.data,
                image.width,
                image.height,
                image.channels as usize,
                quality,
            )?;

            let file = File::create(path).context("Failed to create output file")?;
            let mut writer = BufWriter::new(file);
            writer.write_all(&jpeg_data).context("Failed to write JPEG data")?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::encode_jpeg;

    /// Synthetic starfield: dark background + a few gaussian blobs, the shape of
    /// data this encoder actually sees.
    fn starfield(w: usize, h: usize, channels: usize) -> Vec<u8> {
        let mut v = vec![0u8; w * h * channels];
        let mut s: u32 = 0x9e3779b9;
        let mut rng = || {
            s ^= s << 13;
            s ^= s >> 17;
            s ^= s << 5;
            s
        };
        for px in v.chunks_exact_mut(channels) {
            let n = (rng() % 12) as u8 + 8;
            px[0] = n;
            px[1] = n;
            px[2] = n;
            if channels == 4 {
                px[3] = 255;
            }
        }
        for _ in 0..40 {
            let cx = (rng() as usize) % w;
            let cy = (rng() as usize) % h;
            let peak = 150 + (rng() % 100);
            for dy in 0..7usize {
                for dx in 0..7usize {
                    let (x, y) = (cx + dx, cy + dy);
                    if x >= w || y >= h {
                        continue;
                    }
                    let (fx, fy) = (dx as f32 - 3.0, dy as f32 - 3.0);
                    let g = (-(fx * fx + fy * fy) / 3.0).exp();
                    let val = (peak as f32 * g) as u32;
                    let i = (y * w + x) * channels;
                    for c in 0..3 {
                        v[i + c] = v[i + c].saturating_add(val.min(255) as u8);
                    }
                }
            }
        }
        v
    }

    /// Round-trip guard for the encoder swap: the JPEG must be structurally
    /// valid (SOI/EOI + decodable) and lossy-close to the source. Catches a
    /// wrong pixel format / subsampling / channel-order regression, which a
    /// "does it compile" check would not.
    fn round_trip(channels: usize) {
        let (w, h) = (192usize, 128usize);
        let src = starfield(w, h, channels);

        let jpeg = encode_jpeg(&src, w, h, channels, 92).expect("encode failed");
        assert_eq!(&jpeg[..2], &[0xFF, 0xD8], "missing SOI marker");
        assert_eq!(&jpeg[jpeg.len() - 2..], &[0xFF, 0xD9], "missing EOI marker");

        let decoded = libjpeg_turbo_rs::decompress_to(&jpeg, libjpeg_turbo_rs::PixelFormat::Rgb)
            .expect("decode failed");
        assert_eq!(decoded.width, w);
        assert_eq!(decoded.height, h);

        // Compare RGB triplets, skipping the alpha byte on the RGBA input.
        let mut sse = 0f64;
        for i in 0..(w * h) {
            for c in 0..3 {
                let a = src[i * channels + c] as f64;
                let b = decoded.data[i * 3 + c] as f64;
                sse += (a - b) * (a - b);
            }
        }
        let mse = sse / (w * h * 3) as f64;
        let psnr = 10.0 * (255.0f64 * 255.0 / mse).log10();
        assert!(
            psnr > 35.0,
            "round-trip PSNR too low for {channels}-channel input: {psnr:.1} dB"
        );
    }

    #[test]
    fn jpeg_round_trip_rgb() {
        round_trip(3);
    }

    #[test]
    fn jpeg_round_trip_rgba_ignores_alpha() {
        round_trip(4);
    }

    #[test]
    fn jpeg_rejects_short_buffer() {
        let err = encode_jpeg(&[0u8; 10], 64, 64, 3, 90).unwrap_err();
        assert!(err.to_string().contains("too small"), "got: {err}");
    }
}
