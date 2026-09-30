pub mod fits;
pub mod xisf;

use std::fs::File;
use std::io::Read;
use std::path::Path;

use anyhow::{bail, Result};

use crate::types::{ImageMetadata, PixelData};

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

fn has_fits_extension(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|ext| ["fits", "fit", "fts", "fz"].iter().any(|f| ext.eq_ignore_ascii_case(f)))
}

fn is_fits(path: &Path) -> bool {
    if has_fits_extension(path) {
        return true;
    }
    // Gzipped FITS: `.fits.gz`, `.fit.gz`, `.fts.gz`
    if path.extension().and_then(|e| e.to_str()).is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
        && path.file_stem().is_some_and(|stem| has_fits_extension(Path::new(stem)))
    {
        return true;
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
