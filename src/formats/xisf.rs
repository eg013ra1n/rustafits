use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{bail, Context, Result};
use base64::Engine;
use quick_xml::events::Event;
use quick_xml::reader::Reader;
use rayon::prelude::*;

use crate::types::{BayerPattern, DataType, ImageMetadata, PixelData};

const XISF_SIGNATURE: &[u8; 8] = b"XISF0100";

#[derive(Debug, Clone, Copy, PartialEq)]
enum XisfCompression {
    None,
    Zlib,
    Lz4,
    Lz4hc,
    Zstd,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum XisfSampleFormat {
    Uint8,
    Uint16,
    Uint32,
    Float32,
    Float64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum XisfLocation {
    Attachment,
    Embedded,
}

struct XisfImageInfo {
    width: usize,
    height: usize,
    channels: usize,
    sample_format: XisfSampleFormat,
    is_planar: bool,
    location: XisfLocation,
    attachment_pos: u64,
    block_size: u64,
    compression: XisfCompression,
    uncompressed_size: usize,
    byte_shuffled: bool,
    shuffle_item_size: usize,
    bayer_pattern: BayerPattern,
    flip_vertical: bool,
}

fn parse_sample_format(s: &str) -> XisfSampleFormat {
    match s {
        "UInt8" => XisfSampleFormat::Uint8,
        "UInt16" => XisfSampleFormat::Uint16,
        "UInt32" => XisfSampleFormat::Uint32,
        "Float32" => XisfSampleFormat::Float32,
        "Float64" => XisfSampleFormat::Float64,
        _ => XisfSampleFormat::Float32,
    }
}

fn parse_geometry(s: &str) -> Result<(usize, usize, usize)> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() < 2 {
        bail!("Invalid geometry format: {}", s);
    }
    let width: usize = parts[0].parse().context("Invalid width")?;
    let height: usize = parts[1].parse().context("Invalid height")?;
    let channels: usize = if parts.len() >= 3 {
        parts[2].parse().context("Invalid channels")?
    } else {
        1
    };
    Ok((width, height, channels))
}

fn parse_location(s: &str) -> Result<(XisfLocation, u64, u64)> {
    if s.starts_with("attachment:") {
        let rest = &s[11..];
        let parts: Vec<&str> = rest.split(':').collect();
        if parts.len() < 2 {
            bail!("Invalid attachment location: {}", s);
        }
        let pos: u64 = parts[0].parse().context("Invalid attachment position")?;
        let size: u64 = parts[1].parse().context("Invalid attachment size")?;
        Ok((XisfLocation::Attachment, pos, size))
    } else if s == "embedded" {
        Ok((XisfLocation::Embedded, 0, 0))
    } else {
        bail!("Unsupported XISF location: {}", s);
    }
}

fn parse_compression(s: &str) -> Result<(XisfCompression, usize, bool, usize)> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() < 2 {
        bail!("Invalid compression format: {}", s);
    }

    let codec_str = parts[0];
    let (codec_name, byte_shuffled) = if let Some(stripped) = codec_str.strip_suffix("+sh") {
        (stripped, true)
    } else {
        (codec_str, false)
    };

    let compression = match codec_name.to_ascii_lowercase().as_str() {
        "zlib" => XisfCompression::Zlib,
        "lz4" => XisfCompression::Lz4,
        "lz4hc" | "lz4+hc" => XisfCompression::Lz4hc,
        "zstd" => XisfCompression::Zstd,
        _ => bail!("Unknown compression codec: {}", codec_name),
    };

    let uncompressed_size: usize = parts[1].parse().context("Invalid uncompressed size")?;
    let shuffle_item_size: usize = if parts.len() >= 3 {
        parts[2].parse().unwrap_or(1)
    } else {
        1
    };

    Ok((compression, uncompressed_size, byte_shuffled, shuffle_item_size))
}

fn parse_xisf_xml(xml: &str) -> Result<XisfImageInfo> {
    let mut info = XisfImageInfo {
        width: 0,
        height: 0,
        channels: 1,
        sample_format: XisfSampleFormat::Float32,
        is_planar: true,
        location: XisfLocation::Attachment,
        attachment_pos: 0,
        block_size: 0,
        compression: XisfCompression::None,
        uncompressed_size: 0,
        byte_shuffled: false,
        shuffle_item_size: 1,
        bayer_pattern: BayerPattern::None,
        flip_vertical: false,
    };

    let mut reader = Reader::from_str(xml);
    let mut found_image = false;

    loop {
        match reader.read_event() {
            Ok(Event::Empty(ref e)) | Ok(Event::Start(ref e)) if e.name().as_ref() == b"Image" => {
                found_image = true;
                for attr in e.attributes().flatten() {
                    let key = std::str::from_utf8(attr.key.as_ref()).unwrap_or("");
                    let val = attr.unescape_value().unwrap_or_default();

                    match key {
                        "geometry" => {
                            let (w, h, c) = parse_geometry(&val)?;
                            info.width = w;
                            info.height = h;
                            info.channels = c;
                        }
                        "sampleFormat" => {
                            info.sample_format = parse_sample_format(&val);
                        }
                        "pixelStorage" => {
                            info.is_planar = val.eq_ignore_ascii_case("planar");
                        }
                        "location" => {
                            let (loc, pos, size) = parse_location(&val)?;
                            info.location = loc;
                            info.attachment_pos = pos;
                            info.block_size = size;
                        }
                        "compression" => {
                            let (comp, uncomp_size, shuffled, item_size) =
                                parse_compression(&val)?;
                            info.compression = comp;
                            info.uncompressed_size = uncomp_size;
                            info.byte_shuffled = shuffled;
                            info.shuffle_item_size = item_size;
                        }
                        _ => {}
                    }
                }
                break;
            }
            Ok(Event::Eof) => break,
            Err(e) => bail!("XML parse error: {}", e),
            _ => {}
        }
    }

    if !found_image {
        bail!("No <Image> element found in XISF header");
    }

    if info.width == 0 || info.height == 0 {
        bail!("Invalid XISF image dimensions");
    }

    // If no compression, compute expected size
    if info.compression == XisfCompression::None {
        let bps = match info.sample_format {
            XisfSampleFormat::Uint8 => 1,
            XisfSampleFormat::Uint16 => 2,
            XisfSampleFormat::Uint32 => 4,
            XisfSampleFormat::Float32 => 4,
            XisfSampleFormat::Float64 => 8,
        };
        info.uncompressed_size = info.width * info.height * info.channels * bps;
    }

    // Look for Bayer pattern in XML
    if xml.contains("RGGB") {
        info.bayer_pattern = BayerPattern::Rggb;
    } else if xml.contains("BGGR") {
        info.bayer_pattern = BayerPattern::Bggr;
    } else if xml.contains("GBRG") {
        info.bayer_pattern = BayerPattern::Gbrg;
    } else if xml.contains("GRBG") {
        info.bayer_pattern = BayerPattern::Grbg;
    }

    Ok(info)
}

fn unshuffle_bytes(data: &mut [u8], item_size: usize) {
    if item_size <= 1 || data.is_empty() {
        return;
    }

    let num_items = data.len() / item_size;
    if num_items == 0 {
        return;
    }

    let temp = data.to_vec();

    for i in 0..num_items {
        for b in 0..item_size {
            data[i * item_size + b] = temp[b * num_items + i];
        }
    }
}

fn decompress_block(compressed: &[u8], uncompressed_size: usize, codec: XisfCompression) -> Result<Vec<u8>> {
    match codec {
        XisfCompression::None => Ok(compressed.to_vec()),
        XisfCompression::Zlib => {
            let mut decompressor = flate2::Decompress::new(true);
            let mut output = vec![0u8; uncompressed_size];
            let status = decompressor
                .decompress(compressed, &mut output, flate2::FlushDecompress::Finish)
                .context("zlib decompression failed")?;
            if status != flate2::Status::StreamEnd {
                // Try raw deflate (no zlib header)
                let mut decompressor = flate2::Decompress::new(false);
                decompressor
                    .decompress(compressed, &mut output, flate2::FlushDecompress::Finish)
                    .context("zlib decompression failed (raw)")?;
            }
            Ok(output)
        }
        XisfCompression::Lz4 | XisfCompression::Lz4hc => {
            let output = lz4_flex::decompress(compressed, uncompressed_size)
                .map_err(|e| anyhow::anyhow!("LZ4 decompression failed: {}", e))?;
            Ok(output)
        }
        XisfCompression::Zstd => {
            let mut output = vec![0u8; uncompressed_size];
            let mut decoder = ruzstd::frame_decoder::FrameDecoder::new();
            decoder.reset(&compressed[..]).context("zstd decoder reset failed")?;
            let mut cursor = std::io::Cursor::new(&mut output[..]);
            std::io::copy(&mut decoder, &mut cursor)
                .context("zstd decompression failed")?;
            Ok(output)
        }
    }
}

fn convert_chunk(src: &[u8], dst: &mut [f32], format: XisfSampleFormat) {
    match format {
        XisfSampleFormat::Uint8 => {
            for i in 0..dst.len() {
                dst[i] = src[i] as f32 * 256.0;
            }
        }
        XisfSampleFormat::Uint16 => {
            for i in 0..dst.len() {
                let val = u16::from_le_bytes([src[i * 2], src[i * 2 + 1]]);
                dst[i] = val as f32;
            }
        }
        XisfSampleFormat::Uint32 => {
            for i in 0..dst.len() {
                let off = i * 4;
                let val = u32::from_le_bytes([
                    src[off], src[off + 1], src[off + 2], src[off + 3],
                ]);
                dst[i] = (val >> 16) as f32;
            }
        }
        XisfSampleFormat::Float32 => {
            for i in 0..dst.len() {
                let off = i * 4;
                let val = f32::from_le_bytes([
                    src[off], src[off + 1], src[off + 2], src[off + 3],
                ]);
                dst[i] = val * 65535.0;
            }
        }
        XisfSampleFormat::Float64 => {
            for i in 0..dst.len() {
                let off = i * 8;
                let val = f64::from_le_bytes([
                    src[off], src[off + 1], src[off + 2], src[off + 3],
                    src[off + 4], src[off + 5], src[off + 6], src[off + 7],
                ]);
                dst[i] = (val * 65535.0) as f32;
            }
        }
    }
}

fn convert_to_float32(raw_data: &[u8], info: &XisfImageInfo) -> Vec<f32> {
    let num_samples = info.width * info.height * info.channels;
    let mut float_data = vec![0f32; num_samples];
    let bytes_per_sample = match info.sample_format {
        XisfSampleFormat::Uint8 => 1,
        XisfSampleFormat::Uint16 => 2,
        XisfSampleFormat::Uint32 | XisfSampleFormat::Float32 => 4,
        XisfSampleFormat::Float64 => 8,
    };
    let format = info.sample_format;
    const CHUNK: usize = 65536;
    const PAR_THRESHOLD: usize = CHUNK * 2;

    if num_samples >= PAR_THRESHOLD {
        raw_data.par_chunks(CHUNK * bytes_per_sample)
            .zip(float_data.par_chunks_mut(CHUNK))
            .for_each(|(src, dst)| {
                convert_chunk(src, dst, format);
            });
    } else {
        convert_chunk(raw_data, &mut float_data, format);
    }

    float_data
}

fn convert_normal_to_planar(data: &mut Vec<f32>, width: usize, height: usize, channels: usize) {
    if channels != 3 {
        return;
    }

    let plane_size = width * height;
    let mut temp = vec![0f32; plane_size * 3];
    let (r_plane, rest) = temp.split_at_mut(plane_size);
    let (g_plane, b_plane) = rest.split_at_mut(plane_size);

    const PAR_ROW_THRESHOLD: usize = 128;

    if height >= PAR_ROW_THRESHOLD {
        r_plane.par_chunks_mut(width)
            .zip(g_plane.par_chunks_mut(width))
            .zip(b_plane.par_chunks_mut(width))
            .enumerate()
            .for_each(|(row, ((r_row, g_row), b_row))| {
                let base = row * width * 3;
                for x in 0..r_row.len() {
                    r_row[x] = data[base + x * 3];
                    g_row[x] = data[base + x * 3 + 1];
                    b_row[x] = data[base + x * 3 + 2];
                }
            });
    } else {
        for row in 0..height {
            let base = row * width * 3;
            for x in 0..width {
                r_plane[row * width + x] = data[base + x * 3];
                g_plane[row * width + x] = data[base + x * 3 + 1];
                b_plane[row * width + x] = data[base + x * 3 + 2];
            }
        }
    }

    *data = temp;
}

pub fn read_xisf_image(path: &Path) -> Result<(ImageMetadata, PixelData)> {
    let file = File::open(path).context("Failed to open XISF file")?;
    let mut reader = BufReader::new(file);

    // Read and validate signature
    let mut header = [0u8; 16];
    reader
        .read_exact(&mut header)
        .context("Failed to read XISF header")?;

    if &header[..8] != XISF_SIGNATURE {
        bail!("Invalid XISF signature");
    }

    // Bytes 8-11: XML length (little-endian u32)
    let xml_length =
        u32::from_le_bytes([header[8], header[9], header[10], header[11]]) as usize;

    if xml_length == 0 || xml_length > 100 * 1024 * 1024 {
        bail!("Invalid XISF header length");
    }

    // Read XML header
    let mut xml_bytes = vec![0u8; xml_length];
    reader
        .read_exact(&mut xml_bytes)
        .context("Failed to read XISF XML header")?;
    let xml = String::from_utf8_lossy(&xml_bytes);

    // Parse XML
    let info = parse_xisf_xml(&xml)?;

    // Read pixel data
    let compressed_data = match info.location {
        XisfLocation::Attachment => {
            reader
                .seek(SeekFrom::Start(info.attachment_pos))
                .context("Failed to seek to attachment")?;
            let mut data = vec![0u8; info.block_size as usize];
            reader
                .read_exact(&mut data)
                .context("Failed to read attachment data")?;
            data
        }
        XisfLocation::Embedded => {
            // Find base64 content in XML between <Data ...>...</Data> or after Image element
            let b64_content = extract_embedded_data(&xml)?;
            base64::engine::general_purpose::STANDARD
                .decode(b64_content.as_bytes())
                .context("Base64 decode failed")?
        }
    };

    // Decompress if needed
    let raw_data = if info.compression != XisfCompression::None {
        let mut decompressed =
            decompress_block(&compressed_data, info.uncompressed_size, info.compression)?;

        if info.byte_shuffled && info.shuffle_item_size > 1 {
            unshuffle_bytes(&mut decompressed, info.shuffle_item_size);
        }

        decompressed
    } else {
        compressed_data
    };

    // Convert to float32
    let mut float_data = convert_to_float32(&raw_data, &info);

    // Convert interleaved to planar if needed
    if !info.is_planar && info.channels > 1 {
        convert_normal_to_planar(&mut float_data, info.width, info.height, info.channels);
    }

    let meta = ImageMetadata {
        width: info.width,
        height: info.height,
        channels: info.channels,
        dtype: DataType::Float32,
        bayer_pattern: info.bayer_pattern,
        flip_vertical: info.flip_vertical,
    };

    Ok((meta, PixelData::Float32(float_data)))
}

fn extract_embedded_data(xml: &str) -> Result<String> {
    // Look for <Data> element content
    if let Some(start) = xml.find("<Data") {
        if let Some(gt) = xml[start..].find('>') {
            let content_start = start + gt + 1;
            if let Some(end) = xml[content_start..].find("</") {
                let data = &xml[content_start..content_start + end];
                // Strip whitespace
                let clean: String = data.chars().filter(|c| !c.is_whitespace()).collect();
                return Ok(clean);
            }
        }
    }
    bail!("Failed to find embedded data in XISF")
}
