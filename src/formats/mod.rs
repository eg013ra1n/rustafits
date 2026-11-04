pub mod fits;
pub mod xisf;

use std::fs::File;
use std::io::Read;
use std::path::Path;

use anyhow::{bail, Result};

use crate::types::{BayerPattern, ImageHeader, ImageMetadata, PixelData};

pub fn read_image(path: &Path) -> Result<(ImageMetadata, PixelData)> {
    if is_xisf(path) {
        xisf::read_xisf_image(path)
    } else if is_fits(path) {
        fits::read_fits_image(path)
    } else {
        bail!(
            "Unsupported file format: {}",
            path.extension()
                .and_then(|e| e.to_str())
                .unwrap_or("unknown")
        )
    }
}

/// Read a FITS/XISF header without its pixel data: layout (as `read_image`
/// would report it) and every keyword.
pub fn read_header(path: &Path) -> Result<ImageHeader> {
    if is_xisf(path) {
        xisf::read_xisf_header_info(path)
    } else if is_fits(path) {
        fits::read_fits_header_info(path)
    } else {
        bail!(
            "Unsupported file format: {}",
            path.extension()
                .and_then(|e| e.to_str())
                .unwrap_or("unknown")
        )
    }
}

/// `pattern` as seen from one site further right (`odd_x`: the 2×2 cell's
/// columns swap, RGGB → GRBG) and/or one site further down (`odd_y`: its rows
/// swap, RGGB → GBRG).
pub(crate) fn shift_bayer_pattern(pattern: BayerPattern, odd_x: bool, odd_y: bool) -> BayerPattern {
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

fn is_xisf(path: &Path) -> bool {
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        if ext.eq_ignore_ascii_case("xisf") {
            return true;
        }
    }
    // Check magic bytes
    if let Ok(mut f) = File::open(path) {
        let mut sig = [0u8; 8];
        if f.read_exact(&mut sig).is_ok() && &sig == b"XISF0100" {
            return true;
        }
    }
    false
}

fn is_fits(path: &Path) -> bool {
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        if ext.eq_ignore_ascii_case("fits") || ext.eq_ignore_ascii_case("fit") {
            return true;
        }
    }
    // Check magic bytes: "SIMPLE  ="
    if let Ok(mut f) = File::open(path) {
        let mut buf = [0u8; 9];
        if f.read_exact(&mut buf).is_ok() && &buf == b"SIMPLE  =" {
            return true;
        }
    }
    false
}
