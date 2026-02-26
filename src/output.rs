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

    let color_type = if image.channels == 4 {
        ColorType::Rgba8
    } else {
        ColorType::Rgb8
    };

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
                    color_type.into(),
                )
                .context("PNG encoding failed")?;
        }
        _ => {
            // JPEG doesn't support alpha — strip to RGB if needed
            if image.channels == 4 {
                let rgb: Vec<u8> = image
                    .data
                    .chunks_exact(4)
                    .flat_map(|px| [px[0], px[1], px[2]])
                    .collect();
                let encoder = JpegEncoder::new_with_quality(writer, quality);
                encoder
                    .write_image(
                        &rgb,
                        image.width as u32,
                        image.height as u32,
                        ColorType::Rgb8.into(),
                    )
                    .context("JPEG encoding failed")?;
            } else {
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
    }

    Ok(())
}
