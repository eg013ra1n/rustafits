use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use image::codecs::png::PngEncoder;
use image::{ColorType, ImageEncoder};
use img_parts::jpeg::Jpeg;
use img_parts::ImageEXIF;
use little_exif::exif_tag::ExifTag;
use little_exif::metadata::Metadata;
use little_exif::rational::uR64;

use crate::types::{FitsMetadata, ProcessedImage};

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

            let jpeg_data =
                turbojpeg::compress(tj_image, quality as i32, turbojpeg::Subsamp::Sub2x2)
                    .context("JPEG encoding failed")?;

            if ext == "jpg" || ext == "jpeg" {
                embed_exif(&jpeg_data, path, &image.observational)?;
            } else {
                let file = File::create(path).context("Failed to create output file")?;
                let mut writer = BufWriter::new(file);
                writer
                    .write_all(&jpeg_data)
                    .context("Failed to write JPEG data")?;
            }
        }
    }

    Ok(())
}

fn embed_exif(jpeg_data: &[u8], path: &Path, meta: &FitsMetadata) -> Result<()> {
    let mut exif = Metadata::new();

    if let Some(date) = &meta.date_obs {
        let exif_dt = convert_fits_datetime_to_exif(date);
        if !exif_dt.is_empty() {
            exif.set_tag(ExifTag::DateTimeOriginal(exif_dt.clone()));
            exif.set_tag(ExifTag::CreateDate(exif_dt));
        }
    }

    if let Some(exp) = meta.exposure {
        let nom = exp as u32;
        exif.set_tag(ExifTag::ExposureTime(vec![uR64 {
            nominator: nom,
            denominator: 1,
        }]));
    }

    if let Some(iso) = meta.iso {
        exif.set_tag(ExifTag::ISO(vec![iso as u16]));
    }

    if let Some(apt) = meta.aperture {
        exif.set_tag(ExifTag::FNumber(vec![uR64 {
            nominator: (apt * 100.0) as u32,
            denominator: 100,
        }]));
    }

    if let Some(fl) = meta.aperture {
        exif.set_tag(ExifTag::FocalLength(vec![uR64 {
            nominator: (fl * 100.0) as u32,
            denominator: 100,
        }]));
    }

    if let Some(inst) = &meta.instrument {
        exif.set_tag(ExifTag::Model(inst.clone()));
    }

    if let Some(scope) = &meta.telescope {
        exif.set_tag(ExifTag::Make(scope.clone()));
    }

    if let Some(obj) = &meta.object {
        exif.set_tag(ExifTag::ImageDescription(obj.clone()));
    }

    if let Some(obs) = &meta.observer {
        exif.set_tag(ExifTag::Artist(obs.clone()));
    }

    exif.set_tag(ExifTag::Software(format!(
        "rustafits {}",
        env!("CARGO_PKG_VERSION")
    )));

    let mut jpeg = Jpeg::from_bytes(jpeg_data.to_vec().into()).context("Failed to parse JPEG")?;

    let exif_bytes = exif.encode().context("Failed to encode EXIF")?;

    jpeg.set_exif(Some(exif_bytes.into()));

    let mut file = File::create(path).context("Failed to create output file")?;
    jpeg.encoder()
        .write_to(&mut file)
        .context("Failed to write JPEG with EXIF")?;

    Ok(())
}

fn convert_fits_datetime_to_exif(fits: &str) -> String {
    if fits.is_empty() {
        return String::new();
    }
    let normalized = fits.replace("T", " ").replace("/", " ");
    let parts: Vec<&str> = normalized.split_whitespace().collect();
    if parts.is_empty() {
        return String::new();
    }
    let date_part = parts[0];
    let time_part = if parts.len() > 1 { parts[1] } else { "" };

    let date_fixed: String = date_part
        .chars()
        .map(|c| if c == '-' { ':' } else { c })
        .collect();

    let date_trimmed = date_fixed.trim_end_matches(':');
    if time_part.is_empty() {
        date_trimmed.to_string()
    } else {
        format!("{} {}", date_trimmed, time_part)
    }
}

fn parse_ra_dec_to_degrees(ra: &str, dec: &str) -> Option<(f64, f64)> {
    let ra_deg = parse_sexagesimal_ra(ra)?;
    let dec_deg = parse_sexagesimal_dec(dec)?;
    Some((dec_deg, ra_deg))
}

fn parse_sexagesimal_ra(s: &str) -> Option<f64> {
    let s = s.trim();
    let parts: Vec<&str> = s
        .split(|c| c == 'h' || c == 'H' || c == 'd' || c == 'D' || c == ':' || c == ' ')
        .collect();
    if parts.len() < 3 {
        return None;
    }
    let hours: f64 = parts[0].parse().ok()?;
    let mins: f64 = parts[1].parse().ok()?;
    let secs: f64 = parts[2].parse().ok()?;
    let is_negative = s.contains('-');
    let total = hours + mins / 60.0 + secs / 3600.0;
    Some(if is_negative { -total } else { total } * 15.0)
}

fn parse_sexagesimal_dec(s: &str) -> Option<f64> {
    let s = s.trim();
    let parts: Vec<&str> = s
        .split(|c| c == 'd' || c == 'D' || c == ':' || c == ' ')
        .collect();
    if parts.len() < 3 {
        return None;
    }
    let deg: f64 = parts[0].parse().ok()?;
    let mins: f64 = parts[1].parse().ok()?;
    let secs: f64 = parts[2].parse().ok()?;
    let is_negative = s.starts_with('-');
    let total = deg.abs() + mins / 60.0 + secs / 3600.0;
    Some(if is_negative { -total } else { total })
}
