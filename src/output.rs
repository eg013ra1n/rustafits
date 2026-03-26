use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use image::codecs::png::PngEncoder;
use image::{ColorType, ImageEncoder};

use crate::types::ProcessedImage;

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
            let (pixels, format) = if image.channels == 4 {
                // Strip alpha for JPEG
                let rgb: Vec<u8> = image
                    .data
                    .chunks_exact(4)
                    .flat_map(|px| [px[0], px[1], px[2]])
                    .collect();
                (rgb, turbojpeg::PixelFormat::RGB)
            } else {
                (image.data.clone(), turbojpeg::PixelFormat::RGB)
            };

            let tj_image = turbojpeg::Image {
                pixels: pixels.as_slice(),
                width: image.width,
                pitch: image.width * 3,
                height: image.height,
                format,
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
