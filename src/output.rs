use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use image::codecs::png::{CompressionType, FilterType, PngEncoder};
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

    // Output byte-matches the old turbojpeg build not because of the arguments
    // below but because of the four defaults this shim fixes internally
    // (`CompressParams::new` + `high_level::compress`): IsLow DCT, Annex-K
    // quantization tables scaled by quality, standard Annex-K Huffman tables
    // (no two-pass optimization), and no restart interval — the same choices C
    // libjpeg-turbo makes by default. None of them is part of this function's
    // signature, so the `=0.8.0` pin in Cargo.toml is what protects the parity;
    // a version bump can change any of them silently.
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

/// Output file format, chosen from the output path's extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Jpeg,
    Png,
}

impl OutputFormat {
    /// `.png` → PNG; `.jpg`, `.jpeg` or no extension at all → JPEG (case-insensitive).
    /// Any other extension is an error rather than a JPEG stream under a misleading
    /// name (`out.tif` used to receive JPEG bytes silently).
    pub fn from_path(path: &Path) -> Result<Self> {
        let Some(ext) = path.extension() else {
            return Ok(OutputFormat::Jpeg);
        };
        match ext.to_str().map(str::to_ascii_lowercase).as_deref() {
            Some("png") => Ok(OutputFormat::Png),
            Some("jpg") | Some("jpeg") => Ok(OutputFormat::Jpeg),
            _ => anyhow::bail!(
                "unsupported output format '.{}' (expected .jpg, .jpeg or .png)",
                ext.to_string_lossy()
            ),
        }
    }
}

/// PNG DEFLATE effort. Stretched astro frames are mostly noise, so the levels
/// differ a lot. Measured on a 6252x4176 RGB render (78 MB raw, 4 cores):
/// `None` ~0.1 s / 78 MB, `Fast` ~0.27 s / 77 MB, `Balanced` ~8 s / 30 MB,
/// `High` ~12 s / 28.5 MB.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PngCompression {
    /// Stored, unfiltered: fastest, largest.
    None,
    /// The encoder's default (fast DEFLATE, adaptive filter); byte-identical to
    /// what `save_processed` has always written.
    #[default]
    Fast,
    /// zlib-default DEFLATE with adaptive filtering.
    Balanced,
    /// Maximum DEFLATE effort with adaptive filtering.
    High,
}

impl PngCompression {
    fn encoder_settings(self) -> (CompressionType, FilterType) {
        match self {
            // Filtering only helps DEFLATE; without it, it is wasted work.
            PngCompression::None => (CompressionType::Uncompressed, FilterType::NoFilter),
            PngCompression::Fast => (CompressionType::Fast, FilterType::Adaptive),
            PngCompression::Balanced => (CompressionType::Default, FilterType::Adaptive),
            PngCompression::High => (CompressionType::Best, FilterType::Adaptive),
        }
    }
}

impl std::str::FromStr for PngCompression {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_ascii_lowercase().as_str() {
            "none" => Ok(PngCompression::None),
            "fast" => Ok(PngCompression::Fast),
            "balanced" => Ok(PngCompression::Balanced),
            "high" => Ok(PngCompression::High),
            _ => anyhow::bail!("unknown PNG compression '{s}' (expected none, fast, balanced or high)"),
        }
    }
}

pub fn save_image(image: &ProcessedImage, path: &Path, quality: u8) -> Result<()> {
    save_image_with(image, path, quality, PngCompression::default())
}

/// [`save_image`] with an explicit PNG compression level (ignored for JPEG).
pub fn save_image_with(
    image: &ProcessedImage,
    path: &Path,
    quality: u8,
    png_compression: PngCompression,
) -> Result<()> {
    match OutputFormat::from_path(path)? {
        OutputFormat::Png => {
            let file = File::create(path).context("Failed to create output file")?;
            let writer = BufWriter::new(file);
            let color_type = if image.channels == 4 {
                ColorType::Rgba8
            } else {
                ColorType::Rgb8
            };
            let (compression, filter) = png_compression.encoder_settings();
            let encoder = PngEncoder::new_with_quality(writer, compression, filter);
            encoder
                .write_image(
                    &image.data,
                    image.width as u32,
                    image.height as u32,
                    color_type.into(),
                )
                .context("PNG encoding failed")?;
        }
        OutputFormat::Jpeg => {
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

    fn encode_png(image: &ProcessedImage, level: Option<PngCompression>) -> Vec<u8> {
        let mut buf = Vec::new();
        let color = if image.channels == 4 { ColorType::Rgba8 } else { ColorType::Rgb8 };
        let encoder = match level {
            None => PngEncoder::new(&mut buf),
            Some(l) => {
                let (c, f) = l.encoder_settings();
                PngEncoder::new_with_quality(&mut buf, c, f)
            }
        };
        encoder
            .write_image(&image.data, image.width as u32, image.height as u32, color.into())
            .unwrap();
        buf
    }

    fn synth_image(w: usize, h: usize) -> ProcessedImage {
        ProcessedImage { data: synth_rgb(w, h), width: w, height: h, is_color: true, channels: 3, flip_vertical: false }
    }

    /// The default level must keep writing exactly what `PngEncoder::new`
    /// (the pre-option behaviour) writes.
    #[test]
    fn default_png_compression_is_the_old_encoder_output() {
        let img = synth_image(37, 23);
        assert_eq!(encode_png(&img, Some(PngCompression::default())), encode_png(&img, None));
    }

    #[test]
    fn every_png_level_round_trips_the_pixels() {
        let img = synth_image(37, 23);
        for level in [PngCompression::None, PngCompression::Fast, PngCompression::Balanced, PngCompression::High] {
            let png = encode_png(&img, Some(level));
            let decoded = image::load_from_memory_with_format(&png, image::ImageFormat::Png).unwrap();
            assert_eq!(decoded.into_rgb8().into_raw(), img.data, "{level:?}");
        }
    }

    #[test]
    fn png_compression_parses_its_names() {
        assert_eq!("none".parse::<PngCompression>().unwrap(), PngCompression::None);
        assert_eq!("FAST".parse::<PngCompression>().unwrap(), PngCompression::Fast);
        assert_eq!("balanced".parse::<PngCompression>().unwrap(), PngCompression::Balanced);
        assert_eq!("high".parse::<PngCompression>().unwrap(), PngCompression::High);
        assert!("best".parse::<PngCompression>().is_err());
    }

    #[test]
    fn output_format_follows_the_extension() {
        let f = |p: &str| OutputFormat::from_path(Path::new(p));
        assert_eq!(f("a.png").unwrap(), OutputFormat::Png);
        assert_eq!(f("a.PNG").unwrap(), OutputFormat::Png);
        assert_eq!(f("a.jpg").unwrap(), OutputFormat::Jpeg);
        assert_eq!(f("a.JPEG").unwrap(), OutputFormat::Jpeg);
        assert_eq!(f("dir/a").unwrap(), OutputFormat::Jpeg, "no extension stays JPEG");
        assert!(f("a.tif").is_err());
        assert!(f("a.jpg.bak").is_err());
    }

    #[test]
    fn save_image_refuses_an_unknown_extension_without_writing() {
        let dir = std::env::temp_dir().join(format!("rustafits-ext-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("out.tif");
        let image = ProcessedImage {
            data: synth_rgb(8, 8),
            width: 8,
            height: 8,
            is_color: true,
            channels: 3,
            flip_vertical: false,
        };
        assert!(save_image(&image, &path, 90).is_err());
        assert!(!path.exists(), "nothing may be written under the wrong name");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn encode_jpeg_rejects_bad_channel_count_and_short_buffer() {
        let rgb = synth_rgb(8, 8);
        assert!(encode_jpeg(&rgb, 8, 8, 2, 90).is_err());
        assert!(encode_jpeg(&rgb[..10], 8, 8, 3, 90).is_err());
    }
}
