use anyhow::{Context, Result};
use image::{ImageBuffer, RgbImage, Rgb};
use std::path::Path;
use crate::processor::ProcessedImage;
use ndarray::Array3;

pub struct JpegConfig {
    /// JPEG quality (1-100)
    pub quality: u8,
}

impl Default for JpegConfig {
    fn default() -> Self {
        JpegConfig { quality: 95 }
    }
}

/// Save 16-bit RGB data as TIFF (before stretching)
pub fn save_tiff_16bit<P: AsRef<Path>>(
    rgb_data: &Array3<f32>,
    path: P,
) -> Result<()> {
    let (channels, height, width) = rgb_data.dim();

    if channels != 3 {
        return Err(anyhow::anyhow!("TIFF output requires RGB data (3 channels)"));
    }

    // Convert to u16 and interleaved format
    let mut rgb_bytes = Vec::with_capacity(height * width * 3);

    for y in 0..height {
        for x in 0..width {
            let r = rgb_data[[0, y, x]].round().clamp(0.0, 65535.0) as u16;
            let g = rgb_data[[1, y, x]].round().clamp(0.0, 65535.0) as u16;
            let b = rgb_data[[2, y, x]].round().clamp(0.0, 65535.0) as u16;
            rgb_bytes.push(Rgb([r, g, b]));
        }
    }

    // Create image buffer
    let img_buffer: ImageBuffer<Rgb<u16>, Vec<u16>> = ImageBuffer::from_vec(
        width as u32,
        height as u32,
        rgb_bytes.into_iter().flat_map(|p| vec![p.0[0], p.0[1], p.0[2]]).collect()
    ).context("Failed to create image buffer")?;

    // Save as TIFF
    img_buffer.save(path.as_ref())
        .context("Failed to save TIFF")?;

    Ok(())
}

pub fn save_jpeg<P: AsRef<Path>>(
    image: &ProcessedImage,
    path: P,
    config: JpegConfig,
) -> Result<()> {
    let width = image.width() as u32;
    let height = image.height() as u32;

    // Convert to RGB bytes (interleaved format)
    let rgb_bytes = image.to_rgb_bytes();

    // Create image buffer
    let img_buffer: RgbImage = ImageBuffer::from_raw(width, height, rgb_bytes)
        .context("Failed to create image buffer from processed data")?;

    // Save as JPEG
    let mut output = std::fs::File::create(path.as_ref())
        .context("Failed to create output file")?;

    let encoder = image::codecs::jpeg::JpegEncoder::new_with_quality(&mut output, config.quality);

    img_buffer.write_with_encoder(encoder)
        .context("Failed to encode JPEG")?;

    Ok(())
}
