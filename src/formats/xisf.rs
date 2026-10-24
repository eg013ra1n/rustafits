use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{bail, Context, Result};
use base64::Engine;
use quick_xml::events::Event;
use quick_xml::reader::Reader;
use rayon::prelude::*;

use crate::types::{BayerPattern, DataType, FileFormat, HeaderKeyword, ImageHeader, ImageMetadata, PixelData};

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

/// XISF 1.0 §7.2.2 allows a sample's stored byte order to be declared per
/// `<Image>` via `byteOrder="big"|"little"`; absent, it defaults to little.
#[derive(Debug, Clone, Copy, PartialEq)]
enum XisfByteOrder {
    Little,
    Big,
}

#[derive(Debug, Clone)]
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
    byte_order: XisfByteOrder,
    /// The numeric range `[bounds_lo, bounds_hi]` the stored samples
    /// represent (XISF 1.0 §7.2.2's `bounds` attribute on `<Image>`,
    /// float sample formats only). Default `0:1` — the overwhelmingly
    /// common case, and the identity mapping in `convert_to_float32`.
    bounds_lo: f32,
    bounds_hi: f32,
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

/// `bounds="lo:hi"` (XISF 1.0 §7.2.2) — the value range the stored float
/// samples represent. Malformed input falls back to the default `0:1`
/// rather than failing the whole header parse.
fn parse_bounds(s: &str) -> (f32, f32) {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 2 {
        return (0.0, 1.0);
    }
    match (parts[0].parse::<f32>(), parts[1].parse::<f32>()) {
        (Ok(lo), Ok(hi)) => (lo, hi),
        _ => (0.0, 1.0),
    }
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

fn default_image_info() -> XisfImageInfo {
    XisfImageInfo {
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
        byte_order: XisfByteOrder::Little,
        bounds_lo: 0.0,
        bounds_hi: 1.0,
    }
}

/// Parse one `<Image>` element's attributes into a fresh `XisfImageInfo`.
fn parse_image_attrs(e: &quick_xml::events::BytesStart) -> Result<XisfImageInfo> {
    let mut info = default_image_info();
    for attr in e.attributes().flatten() {
        let key = std::str::from_utf8(attr.key.as_ref()).unwrap_or("");
        // XISF headers are XML 1.0; 1.0 and 1.1 normalization differ
        // only for \x85/\x2028, which cannot appear in these values.
        let val = attr
            .normalized_value(quick_xml::XmlVersion::Implicit1_0)
            .unwrap_or_default();

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
                let (comp, uncomp_size, shuffled, item_size) = parse_compression(&val)?;
                info.compression = comp;
                info.uncompressed_size = uncomp_size;
                info.byte_shuffled = shuffled;
                info.shuffle_item_size = item_size;
            }
            "byteOrder" => {
                info.byte_order = if val.eq_ignore_ascii_case("big") {
                    XisfByteOrder::Big
                } else {
                    XisfByteOrder::Little
                };
            }
            "bounds" => {
                let (lo, hi) = parse_bounds(&val);
                info.bounds_lo = lo;
                info.bounds_hi = hi;
            }
            _ => {}
        }
    }
    Ok(info)
}

fn parse_xisf_xml(xml: &str) -> Result<XisfImageInfo> {
    let mut info: Option<XisfImageInfo> = None;

    let mut reader = Reader::from_str(xml);
    // Permissive dangling-`&` retained deliberately (0.36-era behavior; XISF
    // headers are machine-written but a malformed one should still render —
    // cycle decision 2026-07-29).
    reader.config_mut().allow_dangling_amp = true;

    // A header can carry more than one `<Image>` element — a thumbnail
    // before the real image, or (real-world master files written by the
    // external tool) a full-size weight map AFTER the actual pixel data,
    // same geometry. Scan every `<Image>` and keep the one with the
    // largest `width * height * channels`; a tie keeps whichever was seen
    // FIRST (not last) — a same-size weight map always follows the real
    // data in the files this reader has to handle, so "last wins" would
    // silently swap the light data for its own weight map on every
    // drizzled/integration master.
    // A candidate `<Image>` element can itself fail to parse (an
    // unsupported `location` scheme, an unknown compression codec, …) —
    // that must not fail the WHOLE header read when another `<Image>`
    // element is fine; skip it and remember the first error only in case
    // nothing usable turns up at all.
    let mut first_err: Option<anyhow::Error> = None;
    loop {
        match reader.read_event() {
            Ok(Event::Empty(ref e)) | Ok(Event::Start(ref e)) if e.name().as_ref() == b"Image" => {
                let candidate = match parse_image_attrs(e) {
                    Ok(c) => c,
                    Err(err) => {
                        first_err.get_or_insert(err);
                        continue;
                    }
                };
                let candidate_size = candidate.width * candidate.height * candidate.channels;
                let keep = match &info {
                    None => true,
                    Some(current) => {
                        candidate_size > current.width * current.height * current.channels
                    }
                };
                if keep {
                    info = Some(candidate);
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => bail!("XML parse error: {}", e),
            _ => {}
        }
    }

    let mut info = match info {
        Some(info) => info,
        None => match first_err {
            Some(err) => return Err(err),
            None => bail!("No <Image> element found in XISF header"),
        },
    };

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
            // `lz4_flex::decompress` (the crate-root re-export) is deprecated and
            // gated behind `alloc`, which this crate does not enable. Decompress
            // into a buffer of the size the header declares instead — that size is
            // authoritative for XISF, so a short read means a corrupt block.
            let mut output = vec![0u8; uncompressed_size];
            let written = lz4_flex::block::decompress_into(compressed, &mut output)
                .map_err(|e| anyhow::anyhow!("LZ4 decompression failed: {}", e))?;
            if written != uncompressed_size {
                bail!(
                    "LZ4 decompression size mismatch: got {} bytes, header declares {}",
                    written,
                    uncompressed_size
                );
            }
            Ok(output)
        }
        XisfCompression::Zstd => {
            // Must be `StreamingDecoder`, not a bare `FrameDecoder`: the latter's
            // `Read` impl only drains an already-decoded internal buffer and never
            // pulls from the source, so reading it yields 0 bytes and leaves the
            // output zero-filled — a silently black image. `StreamingDecoder` owns
            // the source and drives block decoding on each read.
            let mut output = vec![0u8; uncompressed_size];
            let mut decoder = ruzstd::decoding::StreamingDecoder::new(compressed)
                .map_err(|e| anyhow::anyhow!("zstd decoder init failed: {}", e))?;
            let mut cursor = std::io::Cursor::new(&mut output[..]);
            let written = std::io::copy(&mut decoder, &mut cursor)
                .context("zstd decompression failed")? as usize;
            if written != uncompressed_size {
                bail!(
                    "zstd decompression size mismatch: got {} bytes, header declares {}",
                    written,
                    uncompressed_size
                );
            }
            Ok(output)
        }
    }
}

/// Normalizes a stored float sample against `bounds="lo:hi"` (XISF 1.0
/// §7.2.2) into `[0, 1]` — identity for the default `0:1`, which is what
/// every real file this reader has been checked against uses — THEN
/// applies the u16-like ADU-domain scale every sample format in this
/// reader has always produced (see the "Convention" note on
/// `read_xisf_image` below for why that second step is load-bearing and
/// must not be removed). Degenerate bounds (`hi == lo`) skip
/// normalization rather than dividing by zero.
fn bounds_scale(raw: f32, lo: f32, hi: f32) -> f32 {
    let normalized = if hi != lo && (lo, hi) != (0.0, 1.0) {
        (raw - lo) / (hi - lo)
    } else {
        raw
    };
    normalized * 65535.0
}

/// As `bounds_scale`, but for `Float64` samples: normalize and scale in
/// f64 throughout, casting to f32 only once at the very end. Casting the
/// raw f64 sample down to f32 BEFORE normalizing/scaling (two roundings)
/// is not bit-identical to the pre-M4a reader, which always did this
/// arithmetic in f64 — a subnormal value in particular can round to
/// exactly `0.0` after an extra intermediate f32 cast where a single
/// final cast would not.
fn bounds_scale_f64(raw: f64, lo: f64, hi: f64) -> f32 {
    let normalized = if hi != lo && (lo, hi) != (0.0, 1.0) {
        (raw - lo) / (hi - lo)
    } else {
        raw
    };
    (normalized * 65535.0) as f32
}

fn convert_chunk(src: &[u8], dst: &mut [f32], format: XisfSampleFormat, byte_order: XisfByteOrder, bounds: (f32, f32)) {
    let big = byte_order == XisfByteOrder::Big;
    let (bounds_lo, bounds_hi) = bounds;
    match format {
        XisfSampleFormat::Uint8 => {
            // ×257 (not ×256) so 255 maps to 65535, matching the u16 domain
            // exactly — ×256 tops out at 65280 and renders 8-bit data ~0.4%
            // darker than the same image stored as u16.
            for i in 0..dst.len() {
                dst[i] = src[i] as f32 * 257.0;
            }
        }
        XisfSampleFormat::Uint16 => {
            for i in 0..dst.len() {
                let b = [src[i * 2], src[i * 2 + 1]];
                let val = if big { u16::from_be_bytes(b) } else { u16::from_le_bytes(b) };
                dst[i] = val as f32;
            }
        }
        XisfSampleFormat::Uint32 => {
            for i in 0..dst.len() {
                let off = i * 4;
                let b = [src[off], src[off + 1], src[off + 2], src[off + 3]];
                let val = if big { u32::from_be_bytes(b) } else { u32::from_le_bytes(b) };
                dst[i] = (val >> 16) as f32;
            }
        }
        XisfSampleFormat::Float32 => {
            for i in 0..dst.len() {
                let off = i * 4;
                let b = [src[off], src[off + 1], src[off + 2], src[off + 3]];
                let val = if big { f32::from_be_bytes(b) } else { f32::from_le_bytes(b) };
                dst[i] = bounds_scale(val, bounds_lo, bounds_hi);
            }
        }
        XisfSampleFormat::Float64 => {
            for i in 0..dst.len() {
                let off = i * 8;
                let b = [
                    src[off], src[off + 1], src[off + 2], src[off + 3],
                    src[off + 4], src[off + 5], src[off + 6], src[off + 7],
                ];
                let val = if big { f64::from_be_bytes(b) } else { f64::from_le_bytes(b) };
                dst[i] = bounds_scale_f64(val, bounds_lo as f64, bounds_hi as f64);
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
    let byte_order = info.byte_order;
    let bounds = (info.bounds_lo, info.bounds_hi);
    const CHUNK: usize = 65536;
    const PAR_THRESHOLD: usize = CHUNK * 2;

    if num_samples >= PAR_THRESHOLD {
        raw_data.par_chunks(CHUNK * bytes_per_sample)
            .zip(float_data.par_chunks_mut(CHUNK))
            .for_each(|(src, dst)| {
                convert_chunk(src, dst, format, byte_order, bounds);
            });
    } else {
        convert_chunk(raw_data, &mut float_data, format, byte_order, bounds);
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

/// Validate the signature and return the XML header (reader left after it).
fn read_xml_header(reader: &mut BufReader<File>) -> Result<String> {
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
    Ok(String::from_utf8_lossy(&xml_bytes).into_owned())
}

/// Header-only read: layout of the image `read_xisf_image` would pick, and
/// every `<FITSKeyword>` and `<Property>` in document order. No pixel data.
pub fn read_xisf_header_info(path: &Path) -> Result<ImageHeader> {
    let file = File::open(path).context("Failed to open XISF file")?;
    let mut reader = BufReader::new(file);
    let xml = read_xml_header(&mut reader)?;
    let info = parse_xisf_xml(&xml)?;
    let sample_type = match info.sample_format {
        XisfSampleFormat::Uint8 => "u8",
        XisfSampleFormat::Uint16 => "u16",
        XisfSampleFormat::Uint32 => "u32",
        XisfSampleFormat::Float32 => "f32",
        XisfSampleFormat::Float64 => "f64",
    };
    Ok(ImageHeader {
        format: FileFormat::Xisf,
        width: info.width,
        height: info.height,
        channels: info.channels,
        sample_type,
        bayer_pattern: info.bayer_pattern,
        flip_vertical: info.flip_vertical,
        keywords: parse_xisf_keywords(&xml)?,
    })
}

/// `<FITSKeyword name value comment/>` (value unquoted) and `<Property id
/// value/>` / `<Property id>text</Property>` → keywords, in document order.
fn parse_xisf_keywords(xml: &str) -> Result<Vec<HeaderKeyword>> {
    let mut reader = Reader::from_str(xml);
    reader.config_mut().allow_dangling_amp = true;
    let mut out = Vec::new();
    // A `<Property>` whose value is its text content, awaiting that text.
    let mut open_property: Option<HeaderKeyword> = None;
    loop {
        let event = reader.read_event().map_err(|e| anyhow::anyhow!("XML parse error: {e}"))?;
        match &event {
            Event::Start(e) | Event::Empty(e) => {
                let attr = |name: &[u8]| {
                    e.attributes()
                        .flatten()
                        .find(|a| a.key.as_ref() == name)
                        .and_then(|a| a.normalized_value(quick_xml::XmlVersion::Implicit1_0).ok())
                        .map(|v| v.into_owned())
                };
                match e.name().as_ref() {
                    b"FITSKeyword" => {
                        let value = attr(b"value").unwrap_or_default();
                        let value = value.trim();
                        let value = value
                            .strip_prefix('\'')
                            .and_then(|v| v.strip_suffix('\''))
                            .map_or(value.to_string(), |v| v.replace("''", "'").trim_end().to_string());
                        out.push(HeaderKeyword {
                            key: attr(b"name").unwrap_or_default(),
                            value,
                            comment: attr(b"comment").unwrap_or_default(),
                        });
                    }
                    b"Property" => {
                        let kw = HeaderKeyword {
                            key: attr(b"id").unwrap_or_default(),
                            value: attr(b"value").unwrap_or_default(),
                            comment: attr(b"comment").unwrap_or_default(),
                        };
                        if matches!(event, Event::Start(_)) && kw.value.is_empty() {
                            open_property = Some(kw);
                        } else {
                            out.push(kw);
                        }
                    }
                    _ => {}
                }
            }
            Event::Text(t) => {
                if let Some(kw) = open_property.as_mut() {
                    kw.value.push_str(&t.xml10_content().unwrap_or_default());
                }
            }
            Event::End(e) if e.name().as_ref() == b"Property" => {
                if let Some(mut kw) = open_property.take() {
                    kw.value = kw.value.trim().to_string();
                    out.push(kw);
                }
            }
            Event::Eof => break,
            _ => {}
        }
    }
    Ok(out)
}

/// Convention (M4a Task 1, controller ruling R-M4a-11): every Float32 or
/// Float64 sample this reader returns is in the u16-like ADU domain —
/// `bounds`-normalized to `[0, 1]` THEN multiplied by `65535.0` — the same
/// domain the FITS u16 reader's raw counts are naturally in. This is a
/// PRODUCTION CONTRACT, not just a rendering convenience: `athenaeum-core`'s
/// `integration::banded::spill_via_read_raw` (master builds, light
/// calibration) spills `PixelData::Float32` straight into its ADU-domain
/// band scratch, and `analysis::analyzer::analyze_frame` hands it straight
/// to a detector whose thresholds (`star_levels`, `saturation_limit`) are
/// ADU-domain by construction — either consumer fed a native-`[0, 1]` XISF
/// float would be silently wrong (a 65535×-too-small master; a detector
/// that finds nothing). A one-time attempt to remove this scale (M4a Task 1,
/// commit `b7d1d306`) broke nothing in THIS crate's own tests — its
/// `pipeline.rs`'s auto-stretch is provably invariant to a uniform rescale
/// of its float input (see `processing/stretch.rs`'s
/// `unit_range_floats_stretch_like_u16` test) — but was wrong at the
/// `athenaeum-core` call sites above, so it was reverted; DO NOT remove it
/// again without auditing every `ImageConverter::read_raw`/`read_xisf_image`
/// caller across both crates.
///
/// What WAS the actual bug behind the external tool's own masters measuring
/// at `fwhmPx: 4.3464` (`docs/superpowers/research/2026-09-10-m3-acceptance-run.md`
/// finding 1, "not trustworthy"): `athenaeum_core::stacking::measure`'s
/// estimators assume calibrated frames are float32 in `[0, 1]`
/// (`measure::ADU_SCALE`) and apply their OWN ×65535 internally — a caller
/// (`measure_probe.rs`'s XISF branch, and the ORIGINAL `weight_audit.rs`
/// skeleton) that passes this reader's ADU-domain Float32 straight into
/// `measure_plane`/`measure_plane_with_seeds` without dividing by 65535
/// first double-scales every sample. The fix belongs at THOSE call sites
/// (divide by 65535 before measuring, mirroring how they already divide a
/// `Uint16` result), not in this reader.
pub fn read_xisf_image(path: &Path) -> Result<(ImageMetadata, PixelData)> {
    let file = File::open(path).context("Failed to open XISF file")?;
    let mut reader = BufReader::new(file);
    let xml = read_xml_header(&mut reader)?;

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

#[cfg(test)]
mod tests {
    use super::*;

    /// A 128-byte payload compressed by the `zstd` CLI (an encoder independent
    /// of ruzstd, so a symmetric round-trip bug cannot hide here).
    const ZSTD_FRAME: &[u8] = &[
        0x28, 0xb5, 0x2f, 0xfd, 0x24, 0x80, 0x01, 0x04, 0x00, 0x00, 0x00, 0xf3, 0x6e, 0xe6, 0xdd,
        0xda, 0x4c, 0xcd, 0xbb, 0xc0, 0x2a, 0xb4, 0x99, 0xa7, 0x08, 0x9b, 0x77, 0x8e, 0xe6, 0x81,
        0x55, 0x75, 0xc4, 0x68, 0x33, 0x5b, 0xa2, 0x4f, 0x11, 0x42, 0x80, 0x36, 0xef, 0x29, 0x5e,
        0x1c, 0xcd, 0x10, 0x3c, 0x03, 0xab, 0xf7, 0x19, 0xea, 0x88, 0xdd, 0xf7, 0xd1, 0x66, 0xc4,
        0xd5, 0xb7, 0x44, 0xab, 0xb3, 0x9e, 0x22, 0x92, 0x91, 0x85, 0x00, 0x78, 0x6f, 0x6c, 0xde,
        0x5f, 0x4d, 0x53, 0xbc, 0x46, 0x2b, 0x39, 0x9a, 0x2d, 0x09, 0x20, 0x78, 0x13, 0xe7, 0x07,
        0x56, 0xfa, 0xc4, 0xee, 0x33, 0xe1, 0xa2, 0xd4, 0x11, 0xc8, 0x80, 0xbb, 0xef, 0xae, 0x5e,
        0xa2, 0xcd, 0x95, 0x3c, 0x89, 0xab, 0x7c, 0x1a, 0x6f, 0x89, 0x63, 0xf8, 0x56, 0x67, 0x4a,
        0xd6, 0x3d, 0x45, 0x30, 0xb4, 0x24, 0x23, 0x17, 0x92, 0x0a, 0x01, 0xfe, 0x6f, 0xf1, 0xde,
        0xe5, 0x4d, 0x96, 0x04, 0x31, 0x7d,
    ];

    /// The bytes ZSTD_FRAME encodes: 64 u16 values, little-endian.
    fn expected_payload() -> Vec<u8> {
        (0..64u64)
            .flat_map(|i| (((i * 2_654_435_761) >> 7) as u16).to_le_bytes())
            .collect()
    }

    /// Regression: a bare `FrameDecoder` never pulls from its source, so the
    /// Zstd arm used to return a zero-filled buffer and report success — every
    /// Zstd-compressed XISF silently decoded to a black frame.
    #[test]
    fn zstd_block_decompresses_to_the_real_bytes() {
        let expected = expected_payload();
        let got = decompress_block(ZSTD_FRAME, expected.len(), XisfCompression::Zstd)
            .expect("zstd decompression should succeed");

        assert_eq!(got.len(), expected.len(), "decompressed length");
        assert!(
            got.iter().any(|&b| b != 0),
            "decompressed to all zeros — the decoder is not being driven"
        );
        assert_eq!(got, expected, "decompressed payload");
    }

    /// quick-xml 0.41 rejects a lone `&` in text content by default, where 0.37
    /// read through. We opt back into the permissive behavior, so a header
    /// carrying one must still parse — and still yield correct geometry.
    #[test]
    fn dangling_ampersand_in_header_text_still_parses() {
        let xml = concat!(
            r#"<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0">"#,
            r#"<Property id="note">Dark & Flat calibration</Property>"#,
            r#"<Image geometry="32:16:1" sampleFormat="UInt16" colorSpace="Gray""#,
            r#" pixelStorage="Planar" location="attachment:4096:1024"/>"#,
            r#"</xisf>"#,
        );

        let info = parse_xisf_xml(xml)
            .expect("a lone `&` in text content must not fail the header parse");

        // Not just "no error": the element *after* the dangling `&` must still
        // have been reached and read correctly.
        assert_eq!(info.width, 32, "width");
        assert_eq!(info.height, 16, "height");
        assert_eq!(info.channels, 1, "channels");
        assert_eq!(info.sample_format, XisfSampleFormat::Uint16, "sampleFormat");
        assert_eq!(info.attachment_pos, 4096, "attachment position");
        assert_eq!(info.block_size, 1024, "block size");
    }

    /// A block that does not yield exactly the declared size must be an error,
    /// never a partially-filled (silently zero-padded) buffer.
    #[test]
    fn zstd_block_size_mismatch_is_reported() {
        let declared = expected_payload().len() + 64;
        let err = decompress_block(ZSTD_FRAME, declared, XisfCompression::Zstd)
            .expect_err("a short block must not be reported as success");
        assert!(
            err.to_string().contains("size mismatch"),
            "unexpected error: {err}"
        );
    }

    /// A unique scratch path per test process (mirrors `formats/fits.rs`'s
    /// existing test pattern — no `tempfile` dependency in this crate).
    fn scratch_path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "rustafits_test_{name}_{}_{}.xisf",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0)
        ))
    }

    fn write_xisf_header_and_attachments(
        path: &std::path::Path,
        image_tags: &[String],
        attachments: &[(u64, &[u8])],
    ) {
        let xml = format!(
            r#"<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0">{}</xisf>"#,
            image_tags.join(""),
        );
        let mut bytes = Vec::new();
        bytes.extend_from_slice(XISF_SIGNATURE);
        bytes.extend_from_slice(&(xml.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&[0u8; 4]);
        bytes.extend_from_slice(xml.as_bytes());
        for &(pos, data) in attachments {
            let pos = pos as usize;
            assert!(
                bytes.len() <= pos,
                "fixture attachment offset {pos} overlaps the header/XML (currently {} bytes) — bump the offset",
                bytes.len()
            );
            bytes.resize(pos, 0);
            bytes.extend_from_slice(data);
        }
        std::fs::write(path, bytes).expect("failed to write XISF fixture");
    }

    fn f32_le_bytes(values: &[f32]) -> Vec<u8> {
        values.iter().flat_map(|v| v.to_le_bytes()).collect()
    }

    /// The reader's u16-like ADU-domain convention (R-M4a-11): a stored
    /// float sample with the default `bounds="0:1"` comes back multiplied
    /// by 65535 — computed here, not written as literals, so the fixture
    /// assertions below stay correct if the scale constant ever changes.
    fn adu_scaled(values: &[f32]) -> Vec<f32> {
        values.iter().map(|v| v * 65535.0).collect()
    }

    /// header: `<Image geometry="4:2:1" .../>` (a 4x2 Gray thumbnail) then
    /// `<Image geometry="8:4:3" sampleFormat="Float32" colorSpace="RGB"
    /// pixelStorage="Planar" .../>` (the real data). Real multi-`<Image>`
    /// masters from the external tool put the weight map (same geometry as
    /// the data) AFTER the real image, not a smaller thumbnail before it —
    /// this fixture is the harder case (`read_xisf_image` choosing on size
    /// alone would need the true winner to be second, matching neither
    /// "always take the first" nor "always take the last").
    fn write_two_image_fixture(path: &std::path::Path, main: &[f32]) {
        assert_eq!(main.len(), 96, "fixture expects 96 main samples (8x4x3)");
        let thumb_pos: u64 = 4096;
        let thumb_size: u64 = 32; // 8 samples * 4 bytes
        let main_pos: u64 = 8192;
        let main_size: u64 = 384; // 96 samples * 4 bytes
        let tags = vec![
            format!(
                r#"<Image geometry="4:2:1" sampleFormat="Float32" colorSpace="Gray" location="attachment:{thumb_pos}:{thumb_size}"/>"#
            ),
            format!(
                r#"<Image geometry="8:4:3" sampleFormat="Float32" colorSpace="RGB" pixelStorage="Planar" location="attachment:{main_pos}:{main_size}"/>"#
            ),
        ];
        let thumb_data = f32_le_bytes(&[0.0f32; 8]);
        let main_data = f32_le_bytes(main);
        write_xisf_header_and_attachments(
            path,
            &tags,
            &[(thumb_pos, &thumb_data), (main_pos, &main_data)],
        );
    }

    #[test]
    fn reads_the_main_image_attachment_bit_exact_when_a_smaller_image_precedes_it() {
        let path = scratch_path("two_images");
        let main: Vec<f32> = (0..96).map(|i| i as f32 * 0.001).collect();
        write_two_image_fixture(&path, &main);
        let result = read_xisf_image(&path);
        let _ = std::fs::remove_file(&path);
        let (meta, pixels) = result.unwrap();
        assert_eq!((meta.width, meta.height, meta.channels), (8, 4, 3));
        match pixels {
            PixelData::Float32(v) => assert_eq!(
                v,
                adu_scaled(&main),
                "the main image's samples must come back bit-exact"
            ),
            _ => panic!("Float32 expected"),
        }
    }

    /// The real shape found in the external tool's own drizzled masters:
    /// the light data first, a WEIGHT MAP of the IDENTICAL geometry second
    /// (`masterLight_..._mono_drizzle_2x.xisf` carries `drizzle_integration`
    /// then `drizzle_weights`, both `12448:8336:1`). "Largest, ties → last"
    /// would silently swap the light data for its own weight map on every
    /// such file — this pins "ties → first" instead.
    #[test]
    fn a_same_size_image_following_the_real_data_does_not_win_the_tie() {
        let path = scratch_path("tied_images");
        let data_pos: u64 = 4096;
        let data_size: u64 = 32; // 4x2x1 samples * 4 bytes
        let weights_pos: u64 = 8192;
        let weights_size: u64 = 32;
        let tags = vec![
            format!(
                r#"<Image id="drizzle_integration" geometry="4:2:1" sampleFormat="Float32" colorSpace="Gray" location="attachment:{data_pos}:{data_size}"/>"#
            ),
            format!(
                r#"<Image id="drizzle_weights" geometry="4:2:1" sampleFormat="Float32" colorSpace="Gray" location="attachment:{weights_pos}:{weights_size}"/>"#
            ),
        ];
        let data: Vec<f32> = (0..8).map(|i| i as f32 * 0.01).collect();
        let weights = f32_le_bytes(&[1.0f32; 8]);
        let data_bytes = f32_le_bytes(&data);
        write_xisf_header_and_attachments(
            &path,
            &tags,
            &[(data_pos, &data_bytes), (weights_pos, &weights)],
        );
        let result = read_xisf_image(&path);
        let _ = std::fs::remove_file(&path);
        let (_, pixels) = result.unwrap();
        match pixels {
            PixelData::Float32(v) => assert_eq!(
                v,
                adu_scaled(&data),
                "a same-size image following the real data must not replace it"
            ),
            _ => panic!("Float32 expected"),
        }
    }

    /// XISF 1.0 §7.2.2: `byteOrder="big"` may override the little-endian
    /// default per `<Image>`.
    #[test]
    fn byte_order_big_is_honoured_for_float32_samples() {
        let path = scratch_path("byte_order_big");
        let pos: u64 = 4096;
        let size: u64 = 16; // 4 samples * 4 bytes
        let tags = vec![format!(
            r#"<Image geometry="4:1:1" sampleFormat="Float32" colorSpace="Gray" byteOrder="big" location="attachment:{pos}:{size}"/>"#
        )];
        let values = [0.1f32, 0.25, 0.5, 0.75];
        let be_bytes: Vec<u8> = values.iter().flat_map(|v| v.to_be_bytes()).collect();
        write_xisf_header_and_attachments(&path, &tags, &[(pos, &be_bytes)]);
        let result = read_xisf_image(&path);
        let _ = std::fs::remove_file(&path);
        let (_, pixels) = result.unwrap();
        match pixels {
            PixelData::Float32(v) => assert_eq!(v, adu_scaled(&values), "big-endian samples must be honoured"),
            _ => panic!("Float32 expected"),
        }
    }

    /// XISF 1.0 §7.2.2: `bounds="lo:hi"` is the value range the stored
    /// samples represent; a reader must map into `[0, 1]` before any
    /// further domain scale, not treat the raw stored value as already
    /// normalized.
    #[test]
    fn non_default_bounds_are_rescaled_into_unit_range() {
        let path = scratch_path("bounds");
        let pos: u64 = 4096;
        let size: u64 = 16; // 4 samples * 4 bytes
        let tags = vec![format!(
            r#"<Image geometry="4:1:1" sampleFormat="Float32" colorSpace="Gray" bounds="-1:1" location="attachment:{pos}:{size}"/>"#
        )];
        // Stored in [-1, 1]; expected normalized values are (stored+1)/2.
        let stored = [-1.0f32, -0.5, 0.0, 1.0];
        let stored_bytes = f32_le_bytes(&stored);
        write_xisf_header_and_attachments(&path, &tags, &[(pos, &stored_bytes)]);
        let result = read_xisf_image(&path);
        let _ = std::fs::remove_file(&path);
        let (_, pixels) = result.unwrap();
        let expected = adu_scaled(&[0.0f32, 0.25, 0.5, 1.0]);
        match pixels {
            PixelData::Float32(v) => {
                for (got, want) in v.iter().zip(expected.iter()) {
                    assert!(
                        (got - want).abs() < 1e-2,
                        "bounds-normalized value mismatch: got {got}, want {want}"
                    );
                }
            }
            _ => panic!("Float32 expected"),
        }
    }

    /// A candidate `<Image>` that fails to parse (an unsupported
    /// `location` scheme here) must not fail the whole header when
    /// another `<Image>` element parses fine.
    #[test]
    fn an_unparsable_image_element_is_skipped_not_fatal() {
        let xml = concat!(
            r#"<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0">"#,
            r#"<Image geometry="4:2:1" sampleFormat="Float32" colorSpace="Gray" location="attachment:4096:32"/>"#,
            r#"<Image geometry="4:2:1" sampleFormat="Float32" location="path:foo.dat"/>"#,
            r#"</xisf>"#,
        );
        let info = parse_xisf_xml(xml)
            .expect("a later un-parsable <Image> must not fail the whole header parse");
        assert_eq!((info.width, info.height, info.channels), (4, 2, 1));
        assert_eq!(info.attachment_pos, 4096, "must have picked the parsable image");
    }

    /// When the ONLY `<Image>` element fails to parse, the header read
    /// must still error (not silently succeed with no image).
    #[test]
    fn a_header_whose_only_image_is_unparsable_errors() {
        let xml = concat!(
            r#"<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0">"#,
            r#"<Image geometry="4:2:1" sampleFormat="Float32" location="path:foo.dat"/>"#,
            r#"</xisf>"#,
        );
        let err = parse_xisf_xml(xml)
            .expect_err("a header whose only <Image> can't be parsed must error");
        assert!(
            err.to_string().contains("Unsupported XISF location"),
            "unexpected error: {err}"
        );
    }

    /// `Float64` samples must be normalized/scaled in f64 throughout and
    /// cast to f32 only once at the end — casting to f32 before the
    /// arithmetic (two roundings) is not bit-identical to the
    /// pre-existing reader.
    #[test]
    fn float64_samples_keep_f64_arithmetic_until_the_final_cast() {
        let path = scratch_path("float64_precision");
        let pos: u64 = 4096;
        let size: u64 = 8; // 1 sample * 8 bytes
        let tags = vec![format!(
            r#"<Image geometry="1:1:1" sampleFormat="Float64" colorSpace="Gray" location="attachment:{pos}:{size}"/>"#
        )];
        let stored: f64 = 0.123456789012345;
        let bytes = stored.to_le_bytes().to_vec();
        write_xisf_header_and_attachments(&path, &tags, &[(pos, &bytes)]);
        let result = read_xisf_image(&path);
        let _ = std::fs::remove_file(&path);
        let (_, pixels) = result.unwrap();
        let expected = (stored * 65535.0) as f32;
        match pixels {
            PixelData::Float32(v) => assert_eq!(
                v[0], expected,
                "Float64 sample must scale in f64 arithmetic, one final cast"
            ),
            _ => panic!("Float32 expected"),
        }
    }

    /// A header with zero `<Image>` elements must error (the pin the
    /// M4a Task 1 review flagged as untested).
    #[test]
    fn a_header_with_no_image_elements_errors() {
        let xml = concat!(
            r#"<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0">"#,
            r#"</xisf>"#,
        );
        let err = parse_xisf_xml(xml).expect_err("a header with zero <Image> elements must error");
        assert!(
            err.to_string().contains("No <Image> element found"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn xisf_keywords_cover_fits_keywords_and_both_property_forms() {
        let xml = r#"<xisf><Image geometry="4:3:1" sampleFormat="UInt16">
            <FITSKeyword name="FILTER" value="'Ha      '" comment="Filter name"/>
            <FITSKeyword name="EXPTIME" value="300." comment="[s]"/>
            <Property id="Instrument:Camera:Gain" type="Float32" value="100"/>
            <Property id="Observation:Object:Name" type="String">M 31</Property>
        </Image></xisf>"#;
        let kws = parse_xisf_keywords(xml).unwrap();
        let get = |k: &str| kws.iter().find(|kw| kw.key == k).map(|kw| kw.value.as_str());
        assert_eq!(get("FILTER"), Some("Ha"));
        assert_eq!(get("EXPTIME"), Some("300."));
        assert_eq!(get("Instrument:Camera:Gain"), Some("100"));
        assert_eq!(get("Observation:Object:Name"), Some("M 31"));
        assert_eq!(kws[0].comment, "Filter name");
    }

    #[test]
    fn xisf_header_info_agrees_with_the_pixel_reader_on_the_fixture() {
        let path = std::path::Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/test.xisf"));
        if !path.exists() {
            return;
        }
        let info = read_xisf_header_info(path).unwrap();
        let (meta, _) = read_xisf_image(path).unwrap();
        assert_eq!((info.width, info.height, info.channels), (meta.width, meta.height, meta.channels));
        assert_eq!((info.bayer_pattern, info.flip_vertical), (meta.bayer_pattern, meta.flip_vertical));
        assert_eq!(info.sample_type, "u16");
        assert!(!info.keywords.is_empty());
    }
}
