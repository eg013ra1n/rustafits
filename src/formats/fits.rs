use std::borrow::Cow;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::fs::File;
use std::path::Path;

use anyhow::{bail, Context, Result};
use fitsio_pure::hdu::{FitsData, Hdu, HduInfo};
use rayon::prelude::*;

use crate::types::{BayerPattern, DataType, FileFormat, HeaderKeyword, ImageHeader, ImageMetadata, PixelData};

const FITS_BLOCK_SIZE: usize = 2880;
const FITS_CARD_SIZE: usize = 80;

struct FitsHeader {
    bitpix: i32,
    naxis: i32,
    naxis1: usize,
    naxis2: usize,
    naxis3: usize,
    bzero: f64,
    bscale: f64,
    bayerpat: String,
    xbayroff: i32,
    ybayroff: i32,
    roworder: String,
}

fn get_keyword_value(card: &str, keyword: &str) -> Option<String> {
    if !card.starts_with(keyword) {
        return None;
    }
    let eq_pos = card.find('=')?;
    let val = card[eq_pos + 1..].trim_start();
    Some(val.to_string())
}

fn parse_int_keyword(card: &str, keyword: &str) -> Option<i32> {
    let val = get_keyword_value(card, keyword)?;
    // Take chars until non-digit (or minus sign)
    let num_str: String = val
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '+')
        .collect();
    num_str.parse().ok()
}

fn parse_float_keyword(card: &str, keyword: &str) -> Option<f64> {
    let val = get_keyword_value(card, keyword)?;
    let num_str: String = val
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '+' || *c == '.' || *c == 'E' || *c == 'e')
        .collect();
    num_str.parse().ok()
}

fn parse_string_keyword(card: &str, keyword: &str) -> Option<String> {
    let val = get_keyword_value(card, keyword)?;
    let val = val.trim_start_matches('\'');
    let end = val.find('\'')?;
    let s = val[..end].trim_end().to_string();
    Some(s)
}

/// Read a header up to END. With `cards_out`, every card before END is also
/// collected (as the parser saw it) for [`read_fits_header_info`].
fn read_fits_header<R: Read>(reader: &mut R, mut cards_out: Option<&mut Vec<String>>) -> Result<FitsHeader> {
    let mut hdr = FitsHeader {
        bitpix: 0,
        naxis: 0,
        naxis1: 0,
        naxis2: 0,
        naxis3: 0,
        bzero: 0.0,
        bscale: 1.0,
        bayerpat: String::new(),
        xbayroff: 0,
        ybayroff: 0,
        roworder: String::new(),
    };

    let mut block = [0u8; FITS_BLOCK_SIZE];
    let mut found_end = false;

    while !found_end {
        reader
            .read_exact(&mut block)
            .context("Failed to read FITS header block")?;

        for i in 0..(FITS_BLOCK_SIZE / FITS_CARD_SIZE) {
            let card_bytes = &block[i * FITS_CARD_SIZE..(i + 1) * FITS_CARD_SIZE];
            let card = std::str::from_utf8(card_bytes).unwrap_or("");

            if card.starts_with("END") && card.as_bytes().get(3).map_or(true, |&b| b == b' ') {
                found_end = true;
                break;
            }
            if let Some(cards) = cards_out.as_deref_mut() {
                cards.push(card.to_string());
            }

            if let Some(v) = parse_int_keyword(card, "BITPIX  ") {
                hdr.bitpix = v;
            } else if let Some(v) = parse_int_keyword(card, "NAXIS   ") {
                hdr.naxis = v;
            } else if let Some(v) = parse_int_keyword(card, "NAXIS1") {
                hdr.naxis1 = v as usize;
            } else if let Some(v) = parse_int_keyword(card, "NAXIS2") {
                hdr.naxis2 = v as usize;
            } else if let Some(v) = parse_int_keyword(card, "NAXIS3") {
                hdr.naxis3 = v as usize;
            } else if let Some(v) = parse_float_keyword(card, "BZERO") {
                hdr.bzero = v;
            } else if let Some(v) = parse_float_keyword(card, "BSCALE") {
                hdr.bscale = v;
            } else if let Some(v) = parse_string_keyword(card, "BAYERPAT") {
                hdr.bayerpat = v;
            } else if let Some(v) = parse_int_keyword(card, "XBAYROFF") {
                hdr.xbayroff = v;
            } else if let Some(v) = parse_int_keyword(card, "YBAYROFF") {
                hdr.ybayroff = v;
            } else if let Some(v) = parse_string_keyword(card, "ROWORDER") {
                hdr.roworder = v;
            }
        }
    }

    if hdr.bitpix == 0 {
        bail!("Missing BITPIX keyword in FITS header");
    }

    Ok(hdr)
}

/// The image dimensions a header must carry for the pixel and header readers.
fn check_dims(hdr: &FitsHeader) -> Result<()> {
    if hdr.naxis < 2 {
        bail!("FITS image must have at least 2 dimensions");
    }
    if hdr.naxis1 == 0 || hdr.naxis2 == 0 {
        bail!("Invalid FITS image dimensions");
    }
    Ok(())
}

/// The primary header, when the file is plain FITS (not gzip) with its image
/// in the primary HDU: the in-place reader's case. `None` for everything left
/// to fitsio-pure ([`image_hdu`]): gzip, and an empty primary (`NAXIS = 0`)
/// with the image, plain or tile-compressed, in an extension.
fn plain_primary_header(
    reader: &mut BufReader<File>,
    cards_out: Option<&mut Vec<String>>,
) -> Result<Option<FitsHeader>> {
    if fitsio_pure::gzip::is_gzip(reader.fill_buf().context("Failed to read FITS file")?) {
        return Ok(None);
    }
    let hdr = read_fits_header(reader, cards_out)?;
    Ok((hdr.naxis != 0).then_some(hdr))
}

/// Channels, Bayer pattern and display flip implied by a header — shared by
/// the pixel reader and the header-only reader so they cannot disagree.
fn layout(hdr: &FitsHeader) -> Result<(usize, BayerPattern, bool)> {
    let channels = if hdr.naxis >= 3 && hdr.naxis3 > 0 {
        if hdr.naxis3 != 1 && hdr.naxis3 != 3 {
            bail!(
                "unsupported NAXIS3={}: expected 1 (mono) or 3 (RGB)",
                hdr.naxis3
            );
        }
        hdr.naxis3
    } else {
        1
    };

    let bayer_pattern = match hdr.bayerpat.as_str() {
        "RGGB" => BayerPattern::Rggb,
        "BGGR" => BayerPattern::Bggr,
        "GBRG" => BayerPattern::Gbrg,
        "GRBG" => BayerPattern::Grbg,
        _ => BayerPattern::None,
    };
    // BAYERPAT names the cell that starts XBAYROFF / YBAYROFF sites into the
    // stored frame; fold the offsets in so the pattern addresses pixel (0, 0).
    let bayer_pattern = shift_bayer_pattern(bayer_pattern, hdr.xbayroff % 2 != 0, hdr.ybayroff % 2 != 0);

    // Astronomical convention (PixInsight, ds9, Siril, AstroImageJ):
    //   ROWORDER = 'TOP-DOWN'  → row 0 is the top of the displayed image,
    //                            screen-row-0 = file-row-0, no flip needed.
    //   ROWORDER = 'BOTTOM-UP' → row 0 is the bottom; flip for display.
    //   missing                → astronomical default = bottom-up = flip.
    //
    // (Earlier versions of this code had the comparison inverted, which
    // displayed N.I.N.A. files — explicit TOP-DOWN — flipped relative to
    // PixInsight's display of the same file.)
    let flip_vertical = !hdr.roworder.eq_ignore_ascii_case("TOP-DOWN");

    Ok((channels, bayer_pattern, flip_vertical))
}

/// `pattern` as seen from one site further right (`odd_x`: the 2×2 cell's
/// columns swap, RGGB → GRBG) and/or one site further down (`odd_y`: its rows
/// swap, RGGB → GBRG).
fn shift_bayer_pattern(pattern: BayerPattern, odd_x: bool, odd_y: bool) -> BayerPattern {
    use BayerPattern as P;
    let pattern = if odd_x {
        match pattern {
            P::Rggb => P::Grbg,
            P::Grbg => P::Rggb,
            P::Gbrg => P::Bggr,
            P::Bggr => P::Gbrg,
            P::None => P::None,
        }
    } else {
        pattern
    };
    if odd_y {
        match pattern {
            P::Rggb => P::Gbrg,
            P::Gbrg => P::Rggb,
            P::Grbg => P::Bggr,
            P::Bggr => P::Grbg,
            P::None => P::None,
        }
    } else {
        pattern
    }
}

/// Header-only read: layout and every keyword, no pixel data.
pub fn read_fits_header_info(path: &Path) -> Result<ImageHeader> {
    let file = File::open(path).context("Failed to open FITS file")?;
    let mut reader = BufReader::new(file);
    let mut cards = Vec::new();
    let hdr = match plain_primary_header(&mut reader, Some(&mut cards))? {
        Some(hdr) => hdr,
        None => {
            cards.clear();
            read_fits_header(&mut image_hdu_header(path)?.as_slice(), Some(&mut cards))?
        }
    };
    check_dims(&hdr)?;
    let (channels, bayer_pattern, flip_vertical) = layout(&hdr)?;
    let sample_type = match hdr.bitpix {
        8 => "u8",
        16 if hdr.bzero == 32768.0 && hdr.bscale == 1.0 => "u16",
        16 => "i16",
        32 if hdr.bzero == 2147483648.0 && hdr.bscale == 1.0 => "u32",
        32 => "i32",
        -32 => "f32",
        -64 => "f64",
        other => bail!("Unsupported BITPIX value: {}", other),
    };
    Ok(ImageHeader {
        format: FileFormat::Fits,
        width: hdr.naxis1,
        height: hdr.naxis2,
        channels,
        sample_type,
        bayer_pattern,
        flip_vertical,
        keywords: cards.iter().filter_map(|c| parse_card(c)).collect(),
    })
}

/// One 80-column card → keyword. `KEY     = value / comment`; string values
/// are unquoted (`''` → `'`); cards without a value indicator (COMMENT,
/// HISTORY, …) keep their text as the value; blank cards are dropped.
fn parse_card(card: &str) -> Option<HeaderKeyword> {
    let key = card.get(..8).unwrap_or(card).trim_end();
    if key.is_empty() {
        return None;
    }
    let rest = card.get(8..).unwrap_or("");
    let (value, comment) = match rest.strip_prefix('=') {
        Some(v) => {
            let v = v.trim_start();
            if let Some(quoted) = v.strip_prefix('\'') {
                let mut value = String::new();
                let mut chars = quoted.chars();
                while let Some(c) = chars.next() {
                    if c != '\'' {
                        value.push(c);
                        continue;
                    }
                    let tail = chars.as_str();
                    if let Some(after) = tail.strip_prefix('\'') {
                        value.push('\'');
                        chars = after.chars();
                    } else {
                        break;
                    }
                }
                let comment = chars.as_str().split_once('/').map_or("", |(_, c)| c.trim());
                (value.trim_end().to_string(), comment.to_string())
            } else {
                let (val, com) = v.split_once('/').unwrap_or((v, ""));
                (val.trim().to_string(), com.trim().to_string())
            }
        }
        None => (rest.trim().to_string(), String::new()),
    };
    Some(HeaderKeyword { key: key.to_string(), value, comment })
}

/// Read a FITS image. Plain FITS with the image in the primary HDU is read in
/// place; gzip (`.fits.gz`), tile-compressed (`.fz`) and extension-HDU images
/// go through fitsio-pure, then the same decoding.
pub fn read_fits_image(path: &Path) -> Result<(ImageMetadata, PixelData)> {
    let file = File::open(path).context("Failed to open FITS file")?;
    let mut reader = BufReader::new(file);

    match plain_primary_header(&mut reader, None)? {
        Some(hdr) => read_pixels(path, &hdr, &mut reader),
        None => {
            let bytes = read_fits_bytes(path)?;
            let (header, data) = image_hdu(&bytes)?;
            let hdr = read_fits_header(&mut header.as_slice(), None)?;
            read_pixels(path, &hdr, &mut &*data)
        }
    }
}

/// Decode the data unit `reader` is positioned at (big-endian, as stored in
/// a plain FITS file) per `hdr`.
fn read_pixels<R: Read>(path: &Path, hdr: &FitsHeader, reader: &mut R) -> Result<(ImageMetadata, PixelData)> {
    check_dims(hdr)?;
    let (channels, bayer_pattern, flip_vertical) = layout(hdr)?;

    let num_pixels = hdr.naxis1 * hdr.naxis2 * channels;
    let bytes_per_pixel = (hdr.bitpix.unsigned_abs() as usize) / 8;
    let data_size = num_pixels * bytes_per_pixel;

    const CHUNK: usize = 65536;
    const PAR_THRESHOLD: usize = CHUNK * 2;
    let use_par = num_pixels >= PAR_THRESHOLD;

    // Same-width samples (BITPIX 16 → u16, ±32 → f32) are read straight into
    // the destination buffer and converted in place, so the raw bytes and the
    // decoded pixels never coexist (52 MB instead of 104 MB transient for a
    // 26 MP 16-bit frame). BITPIX 8 and -64 change width and keep a raw buffer.
    let (dtype, pixels) = match hdr.bitpix {
        16 => {
            let mut u16_data = vec![0u16; num_pixels];
            reader
                .read_exact(as_bytes_mut(&mut u16_data))
                .context("Failed to read FITS data")?;

            // `u16::from_be` on the in-memory sample is the file's big-endian
            // value on any host; the loops below auto-vectorize.
            let bzero = hdr.bzero;
            let bscale = hdr.bscale;
            // BITPIX 16 without BZERO: the stored value is signed. Some writers
            // store unsigned data shifted by -32768 and omit BZERO; more than 1 %
            // negative samples marks that case, read it as BZERO 32768. Otherwise
            // the file is genuinely signed and negatives clamp to 0 (never wrap).
            let offset_data = bzero == 0.0 && bscale == 1.0 && {
                // The sign bit of the big-endian sample, before any conversion.
                let count_neg =
                    |d: &[u16]| -> usize { d.iter().filter(|v| u16::from_be(**v) & 0x8000 != 0).count() };
                let negatives: usize = if use_par {
                    u16_data.par_chunks(CHUNK).map(count_neg).sum()
                } else {
                    count_neg(&u16_data)
                };
                let offset_data = negatives * 100 > num_pixels;
                tracing::debug!(
                    path = %path.display(),
                    negatives,
                    samples = num_pixels,
                    offset_data,
                    "fits bitpix16 bzero0 sign rule"
                );
                offset_data
            };
            let convert: &(dyn Fn(&mut [u16]) + Sync) = if (bzero == 32768.0 && bscale == 1.0) || offset_data {
                // Fast path: signed→unsigned via byte-swap + XOR 0x8000
                &|d: &mut [u16]| d.iter_mut().for_each(|v| *v = u16::from_be(*v) ^ 0x8000)
            } else if bzero == 0.0 && bscale == 1.0 {
                &|d: &mut [u16]| d.iter_mut().for_each(|v| *v = (u16::from_be(*v) as i16).max(0) as u16)
            } else {
                &move |d: &mut [u16]| {
                    for v in d.iter_mut() {
                        let val = u16::from_be(*v) as i16;
                        let scaled = bzero + bscale * val as f64;
                        *v = scaled.clamp(0.0, 65535.0) as u16;
                    }
                }
            };
            if use_par {
                u16_data.par_chunks_mut(CHUNK).for_each(convert);
            } else {
                convert(&mut u16_data);
            }

            (DataType::Uint16, PixelData::Uint16(u16_data))
        }
        -32 | 32 => {
            let mut f32_data = vec![0f32; num_pixels];
            reader
                .read_exact(as_bytes_mut(&mut f32_data))
                .context("Failed to read FITS data")?;

            let bzero = hdr.bzero;
            let bscale = hdr.bscale;
            let is_float = hdr.bitpix == -32;
            let convert = |d: &mut [f32]| {
                for v in d.iter_mut() {
                    let bits = u32::from_be(v.to_bits());
                    let val = if is_float { f32::from_bits(bits) as f64 } else { bits as i32 as f64 };
                    *v = (bzero + bscale * val) as f32;
                }
            };
            if use_par {
                f32_data.par_chunks_mut(CHUNK).for_each(convert);
            } else {
                convert(&mut f32_data);
            }

            (DataType::Float32, PixelData::Float32(f32_data))
        }
        8 => {
            let raw_data = read_raw_bytes(reader, data_size)?;
            let mut u16_data = vec![0u16; num_pixels];
            let bzero = hdr.bzero;
            let bscale = hdr.bscale;

            let convert = |s: &[u8], d: &mut [u16]| {
                for i in 0..d.len() {
                    let scaled = bzero + bscale * s[i] as f64;
                    d[i] = (scaled * 256.0) as u16;
                }
            };
            if use_par {
                raw_data.par_chunks(CHUNK).zip(u16_data.par_chunks_mut(CHUNK)).for_each(|(s, d)| convert(s, d));
            } else {
                convert(&raw_data, &mut u16_data);
            }

            (DataType::Uint16, PixelData::Uint16(u16_data))
        }
        -64 => {
            let raw_data = read_raw_bytes(reader, data_size)?;
            let mut f32_data = vec![0f32; num_pixels];
            let src = &raw_data;
            let bzero = hdr.bzero;
            let bscale = hdr.bscale;

            let convert = |s: &[u8], d: &mut [f32]| {
                for i in 0..d.len() {
                    let off = i * 8;
                    let val = f64::from_be_bytes([
                        s[off], s[off + 1], s[off + 2], s[off + 3],
                        s[off + 4], s[off + 5], s[off + 6], s[off + 7],
                    ]);
                    d[i] = (bzero + bscale * val) as f32;
                }
            };
            if use_par {
                src.par_chunks(CHUNK * 8).zip(f32_data.par_chunks_mut(CHUNK)).for_each(|(s, d)| convert(s, d));
            } else {
                convert(src, &mut f32_data);
            }

            (DataType::Float32, PixelData::Float32(f32_data))
        }
        other => bail!("Unsupported BITPIX value: {}", other),
    };

    let meta = ImageMetadata {
        width: hdr.naxis1,
        height: hdr.naxis2,
        channels,
        dtype,
        bayer_pattern,
        flip_vertical,
    };

    Ok((meta, pixels))
}

fn read_raw_bytes<R: Read>(reader: &mut R, len: usize) -> Result<Vec<u8>> {
    let mut raw = vec![0u8; len];
    reader.read_exact(&mut raw).context("Failed to read FITS data")?;
    Ok(raw)
}

// Files the in-place reader does not take go through fitsio-pure, which finds
// the image HDU and, for a tile-compressed image, decompresses it. Its header
// and its data unit (big-endian, as in a plain file) are then decoded by
// `read_fits_header` and `read_pixels` — the same code, so the same keyword
// parsing, BZERO / BSCALE rules and BITPIX → u16 / f32 conversions.

/// The file's bytes, inflated if it is gzip-compressed (`.fits.gz`).
fn read_fits_bytes(path: &Path) -> Result<Vec<u8>> {
    let bytes = std::fs::read(path).context("Failed to read FITS file")?;
    if fitsio_pure::gzip::is_gzip(&bytes) {
        return fitsio_pure::gzip::decompress(&bytes).context("Failed to decompress gzipped FITS file");
    }
    Ok(bytes)
}

/// The HDU holding the image: the primary if it has data axes, else the first
/// image or tile-compressed image extension.
fn find_image_hdu(fits: &FitsData) -> Result<&Hdu> {
    fits.iter()
        .find(|hdu| match &hdu.info {
            HduInfo::Primary { naxes, .. } | HduInfo::Image { naxes, .. } => !naxes.is_empty(),
            HduInfo::CompressedImage { .. } => true,
            _ => false,
        })
        .context("No image HDU in FITS file")
}

/// The header blocks [`read_fits_header`] reads for `hdu`, given its raw
/// header bytes: those bytes, or a tile-compressed image's logical header.
fn image_header(hdu: &Hdu, raw: &[u8]) -> Vec<u8> {
    if matches!(hdu.info, HduInfo::CompressedImage { .. }) {
        logical_header(raw)
    } else {
        raw.to_vec()
    }
}

/// The image HDU of a (gunzipped) FITS file: its header blocks and its data
/// unit, uncompressed and big-endian.
fn image_hdu(bytes: &[u8]) -> Result<(Vec<u8>, Cow<'_, [u8]>)> {
    let fits = fitsio_pure::hdu::parse_fits(bytes).context("Failed to parse FITS file")?;
    let hdu = find_image_hdu(&fits)?;
    let header = image_header(hdu, &bytes[hdu.header_start..hdu.data_start]);
    let data = if matches!(hdu.info, HduInfo::CompressedImage { .. }) {
        let pixels = fitsio_pure::image::read_image_data(bytes, hdu).context("Failed to decompress FITS image")?;
        Cow::Owned(fitsio_pure::image::serialize_image(&pixels))
    } else {
        Cow::Borrowed(&bytes[hdu.data_start..hdu.data_start + hdu.data_len])
    };
    Ok((header, data))
}

/// [`image_hdu`]'s header alone. Only a gzipped file is read whole; otherwise
/// just the header blocks are.
fn image_hdu_header(path: &Path) -> Result<Vec<u8>> {
    let mut file = File::open(path).context("Failed to open FITS file")?;
    let mut magic = [0u8; 2];
    let gzip = file.read_exact(&mut magic).is_ok() && fitsio_pure::gzip::is_gzip(&magic);
    if gzip {
        let bytes = read_fits_bytes(path)?;
        let fits = fitsio_pure::hdu::parse_fits(&bytes).context("Failed to parse FITS file")?;
        let hdu = find_image_hdu(&fits)?;
        return Ok(image_header(hdu, &bytes[hdu.header_start..hdu.data_start]));
    }
    let len = file.metadata().context("Failed to read FITS file")?.len() as usize;
    let mut read_at = |offset: usize, buf: &mut [u8]| -> std::io::Result<()> {
        file.seek(SeekFrom::Start(offset as u64))?;
        file.read_exact(buf)
    };
    let fits = fitsio_pure::hdu::parse_fits_headers(len, |offset, buf| Ok(read_at(offset, buf)?))
        .context("Failed to parse FITS file")?;
    let hdu = find_image_hdu(&fits)?;
    let mut raw = vec![0u8; hdu.data_start - hdu.header_start];
    read_at(hdu.header_start, &mut raw).context("Failed to read FITS header block")?;
    Ok(image_header(hdu, &raw))
}

/// A tile-compressed image's header as the image's own, the way cfitsio and
/// astropy present it: the `Z` keywords that keep the image's structural
/// keywords (`ZBITPIX`, `ZNAXISn`, `ZSIMPLE`, …) renamed back and put first in
/// standard order, the binary table's structure and the compression
/// parameters dropped, and every other card kept in place.
fn logical_header(raw: &[u8]) -> Vec<u8> {
    let mut lead: Vec<(u8, Vec<u8>)> = Vec::new();
    let mut rest: Vec<u8> = Vec::new();
    for card in raw.as_chunks::<FITS_CARD_SIZE>().0 {
        let key = std::str::from_utf8(&card[..8]).unwrap_or("").trim_end();
        if key == "END" {
            break;
        }
        let renamed = |name: &str| {
            let mut c = card.to_vec();
            c[..8].copy_from_slice(format!("{name:<8}").as_bytes());
            c
        };
        match key {
            "ZSIMPLE" => lead.push((0, renamed("SIMPLE"))),
            "ZTENSION" => lead.push((0, renamed("XTENSION"))),
            "ZBITPIX" => lead.push((1, renamed("BITPIX"))),
            "ZNAXIS" => lead.push((2, renamed("NAXIS"))),
            "ZPCOUNT" => lead.push((4, renamed("PCOUNT"))),
            "ZGCOUNT" => lead.push((5, renamed("GCOUNT"))),
            "ZEXTEND" => lead.push((6, renamed("EXTEND"))),
            "ZHECKSUM" => rest.extend(renamed("CHECKSUM")),
            "ZDATASUM" => rest.extend(renamed("DATASUM")),
            _ if indexed(key, "ZNAXIS") => lead.push((3, renamed(&key[1..]))),
            _ if is_table_or_compression_keyword(key) => {}
            _ => rest.extend_from_slice(card),
        }
    }
    // Stable: NAXIS1, NAXIS2, … keep their order.
    lead.sort_by_key(|(rank, _)| *rank);
    let mut out: Vec<u8> = lead.into_iter().flat_map(|(_, card)| card).collect();
    out.extend(rest);
    out.extend(format!("{:<width$}", "END", width = FITS_CARD_SIZE).bytes());
    out.resize(out.len().div_ceil(FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE, b' ');
    out
}

/// `key` is `prefix` followed by an axis / column number.
fn indexed(key: &str, prefix: &str) -> bool {
    key.strip_prefix(prefix).is_some_and(|n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()))
}

/// Keywords of a compressed image's binary table, or of its compression, that
/// describe nothing in the image itself.
fn is_table_or_compression_keyword(key: &str) -> bool {
    matches!(
        key,
        "XTENSION" | "BITPIX" | "NAXIS" | "PCOUNT" | "GCOUNT" | "TFIELDS" | "THEAP" | "CHECKSUM" | "DATASUM"
            | "ZIMAGE" | "ZCMPTYPE" | "ZQUANTIZ" | "ZDITHER0" | "ZMASKCMP" | "ZBLANK" | "ZSCALE" | "ZZERO"
    ) || ["NAXIS", "TTYPE", "TFORM", "TUNIT", "TNULL", "TSCAL", "TZERO", "TDIM", "ZTILE", "ZNAME", "ZVAL"]
        .iter()
        .any(|prefix| indexed(key, prefix))
}

/// Marker for plain sample types that may be filled byte-wise from a file.
trait PlainSample: Copy {}
impl PlainSample for u16 {}
impl PlainSample for f32 {}

/// View a sample buffer as its raw bytes, to `read_exact` into it directly.
fn as_bytes_mut<T: PlainSample>(data: &mut [T]) -> &mut [u8] {
    // SAFETY: u16 and f32 have no padding and every bit pattern is a valid
    // value, u8 has alignment 1, and the length covers exactly the same bytes.
    unsafe { std::slice::from_raw_parts_mut(data.as_mut_ptr().cast::<u8>(), std::mem::size_of_val(data)) }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Formats a single 80-byte FITS header card: an 8-char keyword field,
    /// `= `, then the value, right-padded with spaces to `FITS_CARD_SIZE`.
    fn make_card(keyword: &str, value: impl std::fmt::Display) -> String {
        let mut card = format!("{:<8}= {}", keyword, value);
        card.truncate(FITS_CARD_SIZE);
        format!("{:<width$}", card, width = FITS_CARD_SIZE)
    }

    /// Concatenates cards, appends an END card, then pads with spaces to a
    /// multiple of `FITS_BLOCK_SIZE` (the reader always reads full blocks).
    fn build_header_block(cards: &[String]) -> Vec<u8> {
        let mut bytes: Vec<u8> = cards.iter().flat_map(|c| c.bytes()).collect();
        bytes.extend(format!("{:<width$}", "END", width = FITS_CARD_SIZE).bytes());
        let rem = bytes.len() % FITS_BLOCK_SIZE;
        if rem != 0 {
            bytes.extend(std::iter::repeat(b' ').take(FITS_BLOCK_SIZE - rem));
        }
        bytes
    }

    /// Writes a minimal synthetic 16-bit FITS file (2x2xNAXIS3, zero-filled
    /// pixel data) with the given NAXIS/NAXIS3 to a temp path for testing.
    fn write_synthetic_fits(path: &Path, naxis: i32, naxis3: i32) {
        let cards = vec![
            make_card("SIMPLE", "T"),
            make_card("BITPIX", 16),
            make_card("NAXIS", naxis),
            make_card("NAXIS1", 2),
            make_card("NAXIS2", 2),
            make_card("NAXIS3", naxis3),
        ];
        let mut bytes = build_header_block(&cards);
        let channels = if naxis >= 3 { naxis3.max(1) as usize } else { 1 };
        let pixel_bytes = 2 * 2 * channels * 2; // BITPIX=16 -> 2 bytes/pixel
        bytes.extend(std::iter::repeat(0u8).take(pixel_bytes));
        std::fs::write(path, bytes).expect("failed to write synthetic FITS test file");
    }

    #[test]
    fn read_fits_image_rejects_invalid_naxis3() {
        let path = std::env::temp_dir()
            .join(format!("rustafits_test_naxis3_invalid_{}.fits", std::process::id()));
        write_synthetic_fits(&path, 3, 5);

        let result = read_fits_image(&path);
        let _ = std::fs::remove_file(&path);

        match result {
            Ok(_) => panic!("NAXIS3=5 must be rejected, not silently accepted"),
            Err(err) => assert!(
                err.to_string().contains("NAXIS3=5"),
                "unexpected error message: {err}"
            ),
        }
    }

    #[test]
    fn read_fits_image_accepts_naxis3_rgb() {
        let path = std::env::temp_dir()
            .join(format!("rustafits_test_naxis3_rgb_{}.fits", std::process::id()));
        write_synthetic_fits(&path, 3, 3);

        let result = read_fits_image(&path);
        let _ = std::fs::remove_file(&path);

        assert!(result.is_ok(), "NAXIS3=3 must still parse: {:?}", result.err());
    }

    #[test]
    fn parse_card_reads_values_strings_and_commentary() {
        let card = |s: &str| format!("{s:<80}");
        let kw = parse_card(&card("EXPTIME =                300.0 / [s] Exposure time")).unwrap();
        assert_eq!((kw.key.as_str(), kw.value.as_str(), kw.comment.as_str()), ("EXPTIME", "300.0", "[s] Exposure time"));
        let kw = parse_card(&card("FILTER  = 'Ha      '           / Filter name")).unwrap();
        assert_eq!((kw.value.as_str(), kw.comment.as_str()), ("Ha", "Filter name"));
        let kw = parse_card(&card("OBJECT  = 'O''Neil / A'")).unwrap();
        assert_eq!((kw.value.as_str(), kw.comment.as_str()), ("O'Neil / A", ""), "'' escape and a slash inside the string");
        let kw = parse_card(&card("HISTORY Calibrated with darks")).unwrap();
        assert_eq!((kw.key.as_str(), kw.value.as_str()), ("HISTORY", "Calibrated with darks"));
        assert!(parse_card(&card("")).is_none(), "blank card");
    }

    #[test]
    fn header_info_agrees_with_the_pixel_reader() {
        let mut bytes = build_header_block(&[
            make_card("SIMPLE", "T"),
            make_card("BITPIX", 16),
            make_card("NAXIS", 2),
            make_card("NAXIS1", 4),
            make_card("NAXIS2", 3),
            make_card("BZERO", 32768),
            make_card("BAYERPAT", "'GRBG'"),
            make_card("EXPTIME", "120.0 / seconds"),
        ]);
        bytes.extend(std::iter::repeat(0u8).take(2880));
        let path = std::env::temp_dir().join(format!("rustafits-hdr-{}.fits", std::process::id()));
        std::fs::write(&path, &bytes).unwrap();

        let info = read_fits_header_info(&path).unwrap();
        let (meta, _) = read_fits_image(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!((info.width, info.height, info.channels), (meta.width, meta.height, meta.channels));
        assert_eq!((info.bayer_pattern, info.flip_vertical), (meta.bayer_pattern, meta.flip_vertical));
        assert_eq!(info.sample_type, "u16");
        assert_eq!(info.get("exptime"), Some("120.0"));
        assert_eq!(info.get("BAYERPAT"), Some("GRBG"));
        assert_eq!(info.keywords.len(), 8);
    }

    /// Writes an unsigned 16-bit (BZERO 32768) top-down CFA FITS with the given
    /// BAYERPAT and, when `offsets` is set, XBAYROFF / YBAYROFF.
    fn write_cfa_fits(path: &Path, width: usize, height: usize, bayerpat: &str, offsets: Option<(i32, i32)>, data: &[u16]) {
        assert_eq!(data.len(), width * height);
        let mut cards = vec![
            make_card("SIMPLE", "T"),
            make_card("BITPIX", 16),
            make_card("NAXIS", 2),
            make_card("NAXIS1", width),
            make_card("NAXIS2", height),
            make_card("BZERO", 32768),
            make_card("BAYERPAT", format!("'{bayerpat}'")),
            make_card("ROWORDER", "'TOP-DOWN'"),
        ];
        if let Some((x, y)) = offsets {
            cards.push(make_card("XBAYROFF", format!("{x} / Bayer pattern X axis offset")));
            cards.push(make_card("YBAYROFF", format!("{y} / Bayer pattern Y axis offset")));
        }
        let mut bytes = build_header_block(&cards);
        for v in data {
            bytes.extend(((*v as i32 - 32768) as i16).to_be_bytes());
        }
        bytes.resize(bytes.len().div_ceil(FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE, 0);
        std::fs::write(path, bytes).expect("failed to write synthetic FITS test file");
    }

    #[test]
    fn bayer_offsets_fold_into_the_pattern() {
        use BayerPattern as P;
        // (BAYERPAT, [pattern at (0,0), (1,0), (0,1), (1,1) offsets]): an odd
        // XBAYROFF swaps the cell's columns, an odd YBAYROFF its rows.
        let table = [
            ("RGGB", [P::Rggb, P::Grbg, P::Gbrg, P::Bggr]),
            ("GRBG", [P::Grbg, P::Rggb, P::Bggr, P::Gbrg]),
            ("GBRG", [P::Gbrg, P::Bggr, P::Rggb, P::Grbg]),
            ("BGGR", [P::Bggr, P::Gbrg, P::Grbg, P::Rggb]),
        ];
        let path = std::env::temp_dir().join(format!("rustafits_test_bayroff_{}.fits", std::process::id()));
        for (bayerpat, expected) in table {
            let base = expected[0];
            // Only the parity counts: even and negative offsets behave like 0 / 1.
            let cases = [
                (None, base),
                (Some((0, 0)), expected[0]),
                (Some((1, 0)), expected[1]),
                (Some((0, 1)), expected[2]),
                (Some((1, 1)), expected[3]),
                (Some((2, 4)), expected[0]),
                (Some((3, -1)), expected[3]),
            ];
            for (offsets, want) in cases {
                write_cfa_fits(&path, 4, 2, bayerpat, offsets, &[0; 8]);
                let info = read_fits_header_info(&path).unwrap();
                let (meta, _) = read_fits_image(&path).unwrap();
                assert_eq!(meta.bayer_pattern, want, "BAYERPAT {bayerpat}, offsets {offsets:?}");
                assert_eq!(info.bayer_pattern, want, "header-only read, BAYERPAT {bayerpat}, offsets {offsets:?}");
            }
        }
        std::fs::remove_file(&path).ok();
    }

    /// A star lit only on the red sites of the sensor renders red, whichever
    /// (BAYERPAT, XBAYROFF, YBAYROFF) spelling describes that sensor — with
    /// the super-pixel and the VNG debayer.
    #[test]
    fn red_star_stays_red_for_every_bayer_offset() {
        let (w, h) = (32usize, 32usize);
        let path = std::env::temp_dir().join(format!("rustafits_test_red_star_{}.fits", std::process::id()));
        for bayerpat in ["RGGB", "GRBG", "GBRG", "BGGR"] {
            // Colour (0 = R, 1 = G, 2 = B) of each site of BAYERPAT's cell, by (y & 1) * 2 + (x & 1).
            let cell: Vec<usize> = bayerpat.chars().map(|c| "RGB".find(c).unwrap()).collect();
            for (xoff, yoff) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                // The sensor: pixel (x, y) sits at site (x + XBAYROFF, y + YBAYROFF) of the cell.
                let color_at = |x: usize, y: usize| cell[((y + yoff) & 1) * 2 + ((x + xoff) & 1)];
                let data: Vec<u16> = (0..w * h)
                    .map(|i| {
                        let (x, y) = (i % w, i / w);
                        let star = x.abs_diff(16) <= 3 && y.abs_diff(16) <= 3 && color_at(x, y) == 0;
                        if star { 50000 } else { 1000 + ((x / 2 * 7 + y / 2 * 13) % 50) as u16 }
                    })
                    .collect();
                write_cfa_fits(&path, w, h, bayerpat, Some((xoff as i32, yoff as i32)), &data);

                for (name, converter, centre) in [
                    ("super-pixel", crate::ImageConverter::new(), (8, 8)),
                    ("VNG", crate::ImageConverter::new().with_vng_debayer(), (16, 16)),
                ] {
                    let img = converter.process(&path).unwrap();
                    assert!(img.is_color);
                    let at = (centre.1 * img.width + centre.0) * img.channels as usize;
                    let (r, g, b) = (img.data[at], img.data[at + 1], img.data[at + 2]);
                    assert!(
                        r > g.saturating_add(100) && r > b.saturating_add(100),
                        "{name}, BAYERPAT {bayerpat}, offsets ({xoff}, {yoff}): star renders as RGB ({r}, {g}, {b})"
                    );
                }
            }
        }
        std::fs::remove_file(&path).ok();
    }

    /// Writes a BITPIX 16, BZERO 0 FITS (width x height, mono) with the given
    /// stored samples (big-endian signed).
    fn write_bitpix16_bzero0(path: &Path, width: usize, height: usize, stored: &[i16]) {
        assert_eq!(stored.len(), width * height);
        let cards = vec![
            make_card("SIMPLE", "T"),
            make_card("BITPIX", 16),
            make_card("NAXIS", 2),
            make_card("NAXIS1", width),
            make_card("NAXIS2", height),
            make_card("BZERO", 0),
            make_card("BSCALE", 1),
            make_card("ROWORDER", "'TOP-DOWN'"),
        ];
        let mut bytes = build_header_block(&cards);
        for v in stored {
            bytes.extend(v.to_be_bytes());
        }
        let rem = bytes.len() % FITS_BLOCK_SIZE;
        if rem != 0 {
            bytes.extend(std::iter::repeat(0u8).take(FITS_BLOCK_SIZE - rem));
        }
        std::fs::write(path, bytes).expect("failed to write synthetic FITS test file");
    }

    fn read_u16(path: &Path) -> Vec<u16> {
        let result = read_fits_image(path);
        let _ = std::fs::remove_file(path);
        match result.expect("read failed").1 {
            PixelData::Uint16(v) => v,
            _ => panic!("expected Uint16 pixels"),
        }
    }

    #[test]
    fn bitpix16_bzero0_offset_data_reads_with_32768() {
        let path = std::env::temp_dir()
            .join(format!("rustafits_test_b16z0_offset_{}.fits", std::process::id()));
        // 1000 samples, mostly negative (offset-stored unsigned data).
        let mut stored: Vec<i16> = (0..1000).map(|i| -30000 + i as i16).collect();
        stored[0] = -30000;
        stored[1] = 0;
        stored[2] = 32767;
        write_bitpix16_bzero0(&path, 100, 10, &stored);
        let out = read_u16(&path);
        assert_eq!(out[0], 2768);
        assert_eq!(out[1], 32768);
        assert_eq!(out[2], 65535);
        assert_eq!(out[3], (-30000i32 + 3 + 32768) as u16);
    }

    #[test]
    fn bitpix16_bzero0_signed_data_clamps_negatives() {
        let path = std::env::temp_dir()
            .join(format!("rustafits_test_b16z0_signed_{}.fits", std::process::id()));
        let mut stored: Vec<i16> = (0..1000).map(|i| 100 + (i % 901) as i16).collect();
        for i in [3usize, 100, 400, 700, 999] {
            stored[i] = -5;
        }
        write_bitpix16_bzero0(&path, 100, 10, &stored);
        let out = read_u16(&path);
        for (i, v) in stored.iter().enumerate() {
            let want = (*v).max(0) as u16;
            assert_eq!(out[i], want, "sample {i}");
        }
        assert_eq!(out[3], 0);
    }

    /// Stored samples of a synthetic image, as a plain FITS data unit holds them.
    enum Stored {
        I16(Vec<i16>),
        F32(Vec<f32>),
    }

    impl Stored {
        fn bitpix(&self) -> i32 {
            match self {
                Stored::I16(_) => 16,
                Stored::F32(_) => -32,
            }
        }

        fn be_bytes(&self) -> Vec<u8> {
            match self {
                Stored::I16(v) => v.iter().flat_map(|s| s.to_be_bytes()).collect(),
                Stored::F32(v) => v.iter().flat_map(|s| s.to_be_bytes()).collect(),
            }
        }

        /// One tile-compressed image extension: RICE_1 for integers, lossless
        /// GZIP_2 for floats.
        fn compress(&self, naxes: &[usize], extra: &[fitsio_pure::header::Card]) -> Vec<u8> {
            use fitsio_pure::compress::{compress_image_hdu, TileCompression};
            match self {
                Stored::I16(v) => compress_image_hdu(naxes, v, &TileCompression::rice(), extra),
                Stored::F32(v) => compress_image_hdu(naxes, v, &TileCompression::gzip2().lossless(), extra),
            }
            .expect("tile compression failed")
        }
    }

    /// A synthetic image: axes (NAXIS1 first), the cards after the structural
    /// ones, and the stored samples.
    struct Synthetic {
        naxes: Vec<usize>,
        cards: Vec<String>,
        stored: Stored,
    }

    impl Synthetic {
        /// `first` (SIMPLE or XTENSION), BITPIX, NAXIS, NAXISn.
        fn structural(&self, first: String) -> Vec<String> {
            let mut cards = vec![first, make_card("BITPIX", self.stored.bitpix()), make_card("NAXIS", self.naxes.len())];
            cards.extend(self.naxes.iter().enumerate().map(|(i, n)| make_card(&format!("NAXIS{}", i + 1), n)));
            cards
        }

        fn data_unit(&self) -> Vec<u8> {
            let mut data = self.stored.be_bytes();
            data.resize(data.len().div_ceil(FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE, 0);
            data
        }

        /// The image in the primary HDU: the in-place reader's file.
        fn plain(&self) -> Vec<u8> {
            let mut cards = self.structural(make_card("SIMPLE", "T"));
            cards.extend(self.cards.iter().cloned());
            let mut bytes = build_header_block(&cards);
            bytes.extend(self.data_unit());
            bytes
        }

        /// The image in an IMAGE extension behind an empty primary.
        fn in_extension(&self) -> Vec<u8> {
            let mut cards = self.structural(make_card("XTENSION", "'IMAGE   '"));
            cards.extend([make_card("PCOUNT", 0), make_card("GCOUNT", 1)]);
            cards.extend(self.cards.iter().cloned());
            let mut bytes = empty_primary();
            bytes.extend(build_header_block(&cards));
            bytes.extend(self.data_unit());
            bytes
        }

        /// The same stored samples and keywords, tile-compressed (`.fz`).
        fn tile_compressed(&self) -> Vec<u8> {
            let extra: Vec<_> = self
                .cards
                .iter()
                .map(|c| fitsio_pure::header::parse_card(c.as_bytes().try_into().unwrap()).unwrap())
                .collect();
            let mut bytes = empty_primary();
            bytes.extend(self.stored.compress(&self.naxes, &extra));
            bytes
        }
    }

    fn empty_primary() -> Vec<u8> {
        build_header_block(&[
            make_card("SIMPLE", "T"),
            make_card("BITPIX", 8),
            make_card("NAXIS", 0),
            make_card("EXTEND", "T"),
        ])
    }

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        use std::io::Write;
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(bytes).unwrap();
        enc.finish().unwrap()
    }

    /// Metadata as a comparable tuple (`ImageMetadata` has no `PartialEq`).
    fn meta_tuple(m: &ImageMetadata) -> (usize, usize, usize, DataType, BayerPattern, bool) {
        (m.width, m.height, m.channels, m.dtype, m.bayer_pattern, m.flip_vertical)
    }

    /// Pixels as bit patterns, so NaN compares equal to itself.
    fn pixel_bits(p: &PixelData) -> (DataType, Vec<u32>) {
        match p {
            PixelData::Uint16(v) => (DataType::Uint16, v.iter().map(|&s| s as u32).collect()),
            PixelData::Float32(v) => (DataType::Float32, v.iter().map(|s| s.to_bits()).collect()),
        }
    }

    /// Keywords minus the ones that differ by HDU kind, not by image. Numbers
    /// compare by value: fitsio-pure's writer spells `120.0` `1.2…E2`.
    fn image_keywords(h: &ImageHeader) -> Vec<(String, String)> {
        h.keywords
            .iter()
            .filter(|k| !matches!(k.key.as_str(), "SIMPLE" | "XTENSION" | "PCOUNT" | "GCOUNT"))
            .map(|k| (k.key.clone(), k.value.parse::<f64>().map_or_else(|_| k.value.clone(), |v| v.to_string())))
            .collect()
    }

    /// The gzipped, extension-HDU and tile-compressed (plain and gzipped)
    /// forms of `img` read exactly as its plain form does through the
    /// in-place reader: metadata, PixelData variant, every pixel's bits, and
    /// the header-only read with its keywords.
    fn assert_other_forms_match_plain(name: &str, img: &Synthetic) -> (ImageMetadata, PixelData) {
        let dir = std::env::temp_dir();
        let id = std::process::id();
        let plain_path = dir.join(format!("rustafits_test_{name}_{id}.fits"));
        std::fs::write(&plain_path, img.plain()).unwrap();
        let (want_meta, want_pixels) = read_fits_image(&plain_path).unwrap();
        let want_info = read_fits_header_info(&plain_path).unwrap();
        std::fs::remove_file(&plain_path).ok();

        let fz = img.tile_compressed();
        let forms = [
            ("fits.gz", gzip(&img.plain())),
            ("ext.fits", img.in_extension()),
            ("fz", fz.clone()),
            ("fz.gz", gzip(&fz)),
        ];
        for (suffix, bytes) in forms {
            let path = dir.join(format!("rustafits_test_{name}_{id}.{suffix}"));
            std::fs::write(&path, bytes).unwrap();
            let read = read_fits_image(&path);
            let info = read_fits_header_info(&path);
            std::fs::remove_file(&path).ok();
            let (meta, pixels) = read.unwrap_or_else(|e| panic!("{name}.{suffix}: {e:#}"));
            let info = info.unwrap_or_else(|e| panic!("{name}.{suffix} header: {e:#}"));
            assert_eq!(meta_tuple(&meta), meta_tuple(&want_meta), "{name}.{suffix} metadata");
            assert!(pixel_bits(&pixels) == pixel_bits(&want_pixels), "{name}.{suffix} pixels");
            assert_eq!(
                (info.width, info.height, info.channels, info.sample_type, info.bayer_pattern, info.flip_vertical),
                (
                    want_info.width,
                    want_info.height,
                    want_info.channels,
                    want_info.sample_type,
                    want_info.bayer_pattern,
                    want_info.flip_vertical
                ),
                "{name}.{suffix} header layout"
            );
            assert_eq!(image_keywords(&info), image_keywords(&want_info), "{name}.{suffix} keywords");
        }
        (want_meta, want_pixels)
    }

    #[test]
    fn compressed_and_extension_u16_match_plain() {
        let (w, h) = (37, 23);
        // Unsigned 16-bit over the whole range, stored with BZERO 32768.
        let values: Vec<u16> = (0..w * h).map(|i| (i * 7919 % 65536) as u16).collect();
        let img = Synthetic {
            naxes: vec![w, h],
            cards: vec![make_card("BZERO", 32768), make_card("BSCALE", 1), make_card("EXPTIME", "120.0 / seconds")],
            stored: Stored::I16(values.iter().map(|&v| (v as i32 - 32768) as i16).collect()),
        };
        let (meta, pixels) = assert_other_forms_match_plain("u16", &img);
        assert!(meta.flip_vertical, "no ROWORDER is bottom-up");
        assert!(matches!(pixels, PixelData::Uint16(ref v) if *v == values));
    }

    #[test]
    fn compressed_and_extension_f32_match_plain() {
        let (w, h) = (19, 11);
        let mut values: Vec<f32> = (0..w * h).map(|i| (i as f32 - 50.0) * 13.25).collect();
        values[3] = f32::NAN;
        values[40] = f32::INFINITY;
        values[41] = -0.0;
        values[42] = 1.0e30;
        let img = Synthetic {
            naxes: vec![w, h],
            cards: vec![make_card("BZERO", "100.0"), make_card("BSCALE", "2.0"), make_card("ROWORDER", "'TOP-DOWN'")],
            stored: Stored::F32(values.clone()),
        };
        let (meta, pixels) = assert_other_forms_match_plain("f32", &img);
        assert!(!meta.flip_vertical);
        let PixelData::Float32(out) = pixels else { panic!("expected Float32 pixels") };
        assert!(out[3].is_nan());
        assert_eq!(out[0], 100.0 + 2.0 * values[0]);
    }

    #[test]
    fn compressed_and_extension_osc_match_plain() {
        let (w, h) = (32, 16);
        let img = Synthetic {
            naxes: vec![w, h],
            cards: vec![
                make_card("BZERO", 32768),
                make_card("BAYERPAT", "'RGGB'"),
                make_card("XBAYROFF", 1),
                make_card("YBAYROFF", 0),
                make_card("ROWORDER", "'TOP-DOWN'"),
            ],
            stored: Stored::I16((0..w * h).map(|i| (i as i32 * 97 % 60000 - 30000) as i16).collect()),
        };
        let (meta, _) = assert_other_forms_match_plain("osc", &img);
        assert_eq!(meta.bayer_pattern, BayerPattern::Grbg, "odd XBAYROFF folds RGGB to GRBG");
    }

    #[test]
    fn compressed_and_extension_rgb_match_plain() {
        let naxes = vec![8, 6, 3];
        let img = Synthetic {
            naxes: naxes.clone(),
            cards: vec![make_card("BZERO", 32768)],
            stored: Stored::I16((0..8 * 6 * 3).map(|i| (i * 1201 - 32768) as i16).collect()),
        };
        let (meta, _) = assert_other_forms_match_plain("rgb", &img);
        assert_eq!(meta.channels, 3);
    }

    #[test]
    fn compressed_and_extension_bitpix16_bzero0_match_plain() {
        // No BZERO and mostly negative samples: the offset-data rule applies alike.
        let img = Synthetic {
            naxes: vec![100, 10],
            cards: vec![],
            stored: Stored::I16((0..1000).map(|i| -30000 + i as i16).collect()),
        };
        let (_, pixels) = assert_other_forms_match_plain("b16z0", &img);
        assert!(matches!(pixels, PixelData::Uint16(ref v) if v[0] == 2768));
    }

    #[test]
    fn tile_compressed_header_is_the_images() {
        let img = Synthetic {
            naxes: vec![16, 4],
            cards: vec![make_card("BZERO", 32768), make_card("FILTER", "'Ha'")],
            stored: Stored::I16(vec![0; 64]),
        };
        let path = std::env::temp_dir().join(format!("rustafits_test_fz_header_{}.fz", std::process::id()));
        std::fs::write(&path, img.tile_compressed()).unwrap();
        let info = read_fits_header_info(&path);
        std::fs::remove_file(&path).ok();
        let keys: Vec<String> = info.unwrap().keywords.into_iter().map(|k| k.key).collect();
        assert_eq!(keys, ["BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "BZERO", "FILTER"]);
    }

    #[test]
    fn empty_primary_without_an_image_is_an_error() {
        let path = std::env::temp_dir().join(format!("rustafits_test_no_image_{}.fits", std::process::id()));
        std::fs::write(&path, empty_primary()).unwrap();
        let read = read_fits_image(&path);
        let info = read_fits_header_info(&path);
        std::fs::remove_file(&path).ok();
        assert!(read.is_err() && info.is_err());
    }
}
