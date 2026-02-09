use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};
use image::codecs::jpeg::JpegEncoder;
use image::codecs::png::PngEncoder;
use image::{ColorType, ImageEncoder};

use crate::types::ProcessedImage;

pub fn save_image(image: &ProcessedImage, path: &Path, quality: u8) -> Result<()> {
    let file = File::create(path).context("Failed to create output file")?;
    let writer = BufWriter::new(file);

    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    match ext.as_str() {
        "png" => {
            let encoder = PngEncoder::new(writer);
            encoder
                .write_image(
                    &image.data,
                    image.width as u32,
                    image.height as u32,
                    ColorType::Rgb8.into(),
                )
                .context("PNG encoding failed")?;
        }
        _ => {
            let encoder = JpegEncoder::new_with_quality(writer, quality);
            encoder
                .write_image(
                    &image.data,
                    image.width as u32,
                    image.height as u32,
                    ColorType::Rgb8.into(),
                )
                .context("JPEG encoding failed")?;
        }
    }

    Ok(())
}
