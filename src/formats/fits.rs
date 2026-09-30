use std::borrow::Cow;
use std::path::Path;

use anyhow::{bail, Context, Result};
use fitsio_pure::hdu::{parse_fits, Hdu, HduInfo};
use fitsio_pure::image::{extract_bscale_bzero, image_dimensions, serialize_image};
use fitsio_pure::tiled::read_tiled_image;
use fitsio_pure::value::Value;
use rayon::prelude::*;

use crate::types::{BayerPattern, DataType, ImageMetadata, PixelData};

struct FitsHeader {
    bitpix: i64,
    naxis: usize,
    naxis1: usize,
    naxis2: usize,
    naxis3: usize,
    bzero: f64,
    bscale: f64,
    bayerpat: String,
    roworder: String,
}

fn string_keyword(hdu: &Hdu, keyword: &str) -> String {
    hdu.cards
        .iter()
        .find(|card| card.keyword_str() == keyword)
        .and_then(|card| match &card.value {
            Some(Value::String(s)) => Some(s.trim_end().to_string()),
            _ => None,
        })
        .unwrap_or_default()
}

fn read_fits_header(hdu: &Hdu) -> Result<FitsHeader> {
    let bitpix = match &hdu.info {
        HduInfo::Primary { bitpix, .. } | HduInfo::Image { bitpix, .. } => *bitpix,
        HduInfo::CompressedImage { zbitpix, .. } => *zbitpix,
        _ => bail!("FITS HDU is not an image"),
    };
    let naxes = image_dimensions(hdu)?;
    let (bscale, bzero) = extract_bscale_bzero(&hdu.cards);
    let hdr = FitsHeader {
        bitpix,
        naxis: naxes.len(),
        naxis1: naxes.first().copied().unwrap_or(0),
        naxis2: naxes.get(1).copied().unwrap_or(0),
        naxis3: naxes.get(2).copied().unwrap_or(0),
        bzero,
        bscale,
        bayerpat: string_keyword(hdu, "BAYERPAT"),
        roworder: string_keyword(hdu, "ROWORDER"),
    };

    if hdr.naxis1 == 0 || hdr.naxis2 == 0 {
        bail!("Invalid FITS image dimensions");
    }

    Ok(hdr)
}

/// The image is the primary HDU, or for files whose primary HDU is empty
/// (tile-compressed `.fz` files always are), the first image extension.
fn image_hdu(hdus: &[Hdu]) -> Result<&Hdu> {
    hdus.iter()
        .find(|hdu| image_dimensions(hdu).is_ok_and(|naxes| naxes.len() >= 2))
        .context("FITS image must have at least 2 dimensions")
}

pub fn read_fits_image(path: &Path) -> Result<(ImageMetadata, PixelData)> {
    let mut file = std::fs::read(path).context("Failed to open FITS file")?;
    if fitsio_pure::gzip::is_gzip(&file) {
        file = fitsio_pure::gzip::decompress(&file).context("Failed to decompress FITS file")?;
    }
    let fits = parse_fits(&file).context("Failed to parse FITS file")?;
    let hdu = image_hdu(&fits.hdus)?;

    let hdr = read_fits_header(hdu)?;

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

    let num_pixels = hdr.naxis1 * hdr.naxis2 * channels;
    let bytes_per_pixel = (hdr.bitpix.unsigned_abs() as usize) / 8;
    let data_size = num_pixels * bytes_per_pixel;

    // Tiles decompress to native-endian pixels; serialising them back to
    // big-endian lets both kinds of file share the conversions below.
    let raw_data: Cow<[u8]> = if matches!(hdu.info, HduInfo::CompressedImage { .. }) {
        Cow::Owned(serialize_image(&read_tiled_image(&file, hdu)?))
    } else {
        Cow::Borrowed(&file[hdu.data_start..])
    };
    let raw_data = raw_data.get(..data_size).context("Failed to read FITS data")?;

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
                convert(raw_data, &mut u16_data);
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
            i = unsafe { bswap_u16_ssse3(src, dst, i, n) };
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
            i = unsafe { bswap_u16_xor_ssse3(src, dst, i, n) };
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

#[cfg(test)]
mod tests {
    use super::*;
    use fitsio_pure::{BLOCK_SIZE as FITS_BLOCK_SIZE, CARD_SIZE as FITS_CARD_SIZE};

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
}
