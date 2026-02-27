use std::io::{BufReader, Read};
use std::fs::File;
use std::path::Path;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;

use crate::types::{BayerPattern, DataType, ImageMetadata, PixelData};

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

fn read_fits_header(reader: &mut BufReader<File>) -> Result<FitsHeader> {
    let mut hdr = FitsHeader {
        bitpix: 0,
        naxis: 0,
        naxis1: 0,
        naxis2: 0,
        naxis3: 0,
        bzero: 0.0,
        bscale: 1.0,
        bayerpat: String::new(),
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
            } else if let Some(v) = parse_string_keyword(card, "ROWORDER") {
                hdr.roworder = v;
            }
        }
    }

    if hdr.bitpix == 0 {
        bail!("Missing BITPIX keyword in FITS header");
    }
    if hdr.naxis < 2 {
        bail!("FITS image must have at least 2 dimensions");
    }
    if hdr.naxis1 == 0 || hdr.naxis2 == 0 {
        bail!("Invalid FITS image dimensions");
    }

    Ok(hdr)
}

pub fn read_fits_image(path: &Path) -> Result<(ImageMetadata, PixelData)> {
    let file = File::open(path).context("Failed to open FITS file")?;
    let mut reader = BufReader::new(file);

    let hdr = read_fits_header(&mut reader)?;

    let channels = if hdr.naxis >= 3 && hdr.naxis3 > 0 {
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

    let flip_vertical = hdr.roworder == "TOP-DOWN";

    let num_pixels = hdr.naxis1 * hdr.naxis2 * channels;
    let bytes_per_pixel = (hdr.bitpix.unsigned_abs() as usize) / 8;
    let data_size = num_pixels * bytes_per_pixel;

    let mut raw_data = vec![0u8; data_size];
    reader
        .read_exact(&mut raw_data)
        .context("Failed to read FITS data")?;

    const CHUNK: usize = 65536;
    const PAR_THRESHOLD: usize = CHUNK * 2;

    let (dtype, pixels) = match hdr.bitpix {
        16 => {
            let mut u16_data = vec![0u16; num_pixels];
            let src = &raw_data;
            let use_par = num_pixels >= PAR_THRESHOLD;

            if hdr.bzero == 32768.0 && hdr.bscale == 1.0 {
                // Fast path: signed→unsigned via byte-swap + XOR 0x8000
                let convert = |s: &[u8], d: &mut [u16]| {
                    bswap_u16_xor(s, d);
                };
                if use_par {
                    src.par_chunks(CHUNK * 2).zip(u16_data.par_chunks_mut(CHUNK)).for_each(|(s, d)| convert(s, d));
                } else {
                    convert(src, &mut u16_data);
                }
            } else if hdr.bzero == 0.0 && hdr.bscale == 1.0 {
                let convert = |s: &[u8], d: &mut [u16]| {
                    bswap_u16(s, d);
                };
                if use_par {
                    src.par_chunks(CHUNK * 2).zip(u16_data.par_chunks_mut(CHUNK)).for_each(|(s, d)| convert(s, d));
                } else {
                    convert(src, &mut u16_data);
                }
            } else {
                let bzero = hdr.bzero;
                let bscale = hdr.bscale;
                let convert = |s: &[u8], d: &mut [u16]| {
                    for i in 0..d.len() {
                        let val = i16::from_be_bytes([s[i * 2], s[i * 2 + 1]]);
                        let scaled = bzero + bscale * val as f64;
                        d[i] = scaled.clamp(0.0, 65535.0) as u16;
                    }
                };
                if use_par {
                    src.par_chunks(CHUNK * 2).zip(u16_data.par_chunks_mut(CHUNK)).for_each(|(s, d)| convert(s, d));
                } else {
                    convert(src, &mut u16_data);
                }
            }

            (DataType::Uint16, PixelData::Uint16(u16_data))
        }
        -32 => {
            let mut f32_data = vec![0f32; num_pixels];
            let src = &raw_data;
            let bzero = hdr.bzero;
            let bscale = hdr.bscale;

            let convert = |s: &[u8], d: &mut [f32]| {
                for i in 0..d.len() {
                    let off = i * 4;
                    let val = f32::from_be_bytes([s[off], s[off + 1], s[off + 2], s[off + 3]]);
                    d[i] = (bzero + bscale * val as f64) as f32;
                }
            };
            if num_pixels >= PAR_THRESHOLD {
                src.par_chunks(CHUNK * 4).zip(f32_data.par_chunks_mut(CHUNK)).for_each(|(s, d)| convert(s, d));
            } else {
                convert(src, &mut f32_data);
            }

            (DataType::Float32, PixelData::Float32(f32_data))
        }
        8 => {
            let mut u16_data = vec![0u16; num_pixels];
            let bzero = hdr.bzero;
            let bscale = hdr.bscale;

            let convert = |s: &[u8], d: &mut [u16]| {
                for i in 0..d.len() {
                    let scaled = bzero + bscale * s[i] as f64;
                    d[i] = (scaled * 256.0) as u16;
                }
            };
            if num_pixels >= PAR_THRESHOLD {
                raw_data.par_chunks(CHUNK).zip(u16_data.par_chunks_mut(CHUNK)).for_each(|(s, d)| convert(s, d));
            } else {
                convert(&raw_data, &mut u16_data);
            }

            (DataType::Uint16, PixelData::Uint16(u16_data))
        }
        32 => {
            let mut f32_data = vec![0f32; num_pixels];
            let src = &raw_data;
            let bzero = hdr.bzero;
            let bscale = hdr.bscale;

            let convert = |s: &[u8], d: &mut [f32]| {
                for i in 0..d.len() {
                    let off = i * 4;
                    let val = i32::from_be_bytes([s[off], s[off + 1], s[off + 2], s[off + 3]]);
                    d[i] = (bzero + bscale * val as f64) as f32;
                }
            };
            if num_pixels >= PAR_THRESHOLD {
                src.par_chunks(CHUNK * 4).zip(f32_data.par_chunks_mut(CHUNK)).for_each(|(s, d)| convert(s, d));
            } else {
                convert(src, &mut f32_data);
            }

            (DataType::Float32, PixelData::Float32(f32_data))
        }
        -64 => {
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
            if num_pixels >= PAR_THRESHOLD {
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

/// SIMD-accelerated big-endian u16 byte swap.
fn bswap_u16(src: &[u8], dst: &mut [u16]) {
    let n = dst.len();
    let mut i = 0;

    #[cfg(target_arch = "aarch64")]
    {
        use std::arch::aarch64::*;
        // NEON: vrev16q_u8 swaps adjacent bytes in each 16-bit lane (16 bytes = 8 u16 at a time)
        while i + 8 <= n {
            unsafe {
                let bytes = vld1q_u8(src.as_ptr().add(i * 2));
                let swapped = vrev16q_u8(bytes);
                vst1q_u8(dst.as_mut_ptr().add(i) as *mut u8, swapped);
            }
            i += 8;
        }
    }

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("ssse3") {
            i = bswap_u16_ssse3(src, dst, i, n);
        }
    }

    // Scalar remainder
    for j in i..n {
        dst[j] = u16::from_be_bytes([src[j * 2], src[j * 2 + 1]]);
    }
}

/// SIMD-accelerated big-endian u16 byte swap with XOR 0x8000.
fn bswap_u16_xor(src: &[u8], dst: &mut [u16]) {
    let n = dst.len();
    let mut i = 0;

    #[cfg(target_arch = "aarch64")]
    {
        use std::arch::aarch64::*;
        let xor_mask = unsafe { vdupq_n_u16(0x8000) };
        while i + 8 <= n {
            unsafe {
                let bytes = vld1q_u8(src.as_ptr().add(i * 2));
                let swapped = vrev16q_u8(bytes);
                let vals = vreinterpretq_u16_u8(swapped);
                let result = veorq_u16(vals, xor_mask);
                vst1q_u16(dst.as_mut_ptr().add(i), result);
            }
            i += 8;
        }
    }

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("ssse3") {
            i = bswap_u16_xor_ssse3(src, dst, i, n);
        }
    }

    // Scalar remainder
    for j in i..n {
        let raw = u16::from_be_bytes([src[j * 2], src[j * 2 + 1]]);
        dst[j] = raw ^ 0x8000;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn bswap_u16_ssse3(src: &[u8], dst: &mut [u16], start: usize, n: usize) -> usize {
    use std::arch::x86_64::*;
    let shuffle = _mm_setr_epi8(1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14);
    let mut i = start;
    while i + 8 <= n {
        let bytes = _mm_loadu_si128(src.as_ptr().add(i * 2) as *const __m128i);
        let swapped = _mm_shuffle_epi8(bytes, shuffle);
        _mm_storeu_si128(dst.as_mut_ptr().add(i) as *mut __m128i, swapped);
        i += 8;
    }
    i
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn bswap_u16_xor_ssse3(src: &[u8], dst: &mut [u16], start: usize, n: usize) -> usize {
    use std::arch::x86_64::*;
    let shuffle = _mm_setr_epi8(1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14);
    let xor_mask = _mm_set1_epi16(0x8000u16 as i16);
    let mut i = start;
    while i + 8 <= n {
        let bytes = _mm_loadu_si128(src.as_ptr().add(i * 2) as *const __m128i);
        let swapped = _mm_shuffle_epi8(bytes, shuffle);
        let result = _mm_xor_si128(swapped, xor_mask);
        _mm_storeu_si128(dst.as_mut_ptr().add(i) as *mut __m128i, result);
        i += 8;
    }
    i
}
