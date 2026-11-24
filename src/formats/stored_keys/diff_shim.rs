//! The old header rules as their callers apply them today, for `examples/header_rules_diff.rs` (starloom phase 2,
//! Task 2, ruling P2-R23: before `stored_keys` gives way to `crate::header`, every pinned input and every corpus file
//! is read through both, and the owner sees where they differ). Re-exported as `astroimage::old_header_rules` under
//! the `test-support` feature, `#[doc(hidden)]`, not semver-covered; it goes with `stored_keys`.
//!
//! Each value is read exactly as its caller reads it now: `star_measure::raw::header_facts` (FITS `BITPIX`,
//! `BZERO`, `BSCALE`, `DATAMAX`, `EGAIN` by `FitsHeader`, over the routed image header of a compressed or
//! extension-HDU file; XISF `DATAMAX`, `EGAIN` by `num` over `read_xisf_stored_keys`), the FITS reader's Bayer
//! offsets (`FitsHeader::get_i32`, absent is 0) and the XISF reader's (`offset(num(..))`). The XISF reader's callers
//! read no `BITPIX`, `BZERO` or `BSCALE`; for those three this reads `num`, the rule the XISF text would give them.

use std::path::Path;

use anyhow::Result;

use super::{num, offset, FitsHeader};
use crate::formats::fits::routed_image_header;
use crate::formats::xisf::read_xisf_stored_keys;
use crate::types::{BayerPattern, FileFormat};

/// The keys the comparison reads, in its order.
pub const KEYS: [&str; 7] = ["BITPIX", "BZERO", "BSCALE", "DATAMAX", "EGAIN", "XBAYROFF", "YBAYROFF"];

/// One file read by the old rules.
#[derive(Debug, Clone)]
pub struct OldRead {
    /// The file is XISF (`ImageConverter::read_header`'s format, as `star_measure::raw` decides it).
    pub xisf: bool,
    /// The value of each of [`KEYS`], in order, as [`int`] / [`number`] spell it; `None` is absent or unreadable.
    pub values: Vec<Option<String>>,
    /// The Bayer pattern the reader reports: `BAYERPAT` (FITS) or the XML scan's pattern (XISF), folded with
    /// [`OldRead::offsets`].
    pub bayer: BayerPattern,
    /// The `XBAYROFF` / `YBAYROFF` the reader folded (absent or unreadable is 0).
    pub offsets: (i64, i64),
}

/// An integer as the comparison prints it.
pub fn int(v: Option<i64>) -> Option<String> {
    v.map(|n| n.to_string())
}

/// A number as the comparison prints it: `f64`'s shortest round-trip form, so two strings are equal exactly when
/// the numbers are (`NaN` equals `NaN`).
pub fn number(v: Option<f64>) -> Option<String> {
    v.map(|n| format!("{n:?}"))
}

/// The reader's Bayer fold (`formats::shift_bayer_pattern`).
pub fn shift_bayer_pattern(pattern: BayerPattern, odd_x: bool, odd_y: bool) -> BayerPattern {
    crate::formats::shift_bayer_pattern(pattern, odd_x, odd_y)
}

/// Read `path` by the old rules.
pub fn read(path: &Path) -> Result<OldRead> {
    let info = crate::formats::read_header(path)?;
    let xisf = info.format == FileFormat::Xisf;
    let (values, offsets) = if xisf {
        let map = read_xisf_stored_keys(path)?;
        let n = |k: &str| number(num(&map, k));
        let (x, y) = (offset(num(&map, "XBAYROFF")), offset(num(&map, "YBAYROFF")));
        (
            vec![n("BITPIX"), n("BZERO"), n("BSCALE"), n("DATAMAX"), n("EGAIN"), int(Some(x)), int(Some(y))],
            (x, y),
        )
    } else {
        let h = match routed_image_header(path)? {
            Some(blocks) => FitsHeader::from_bytes(&blocks)?,
            None => FitsHeader::from_path(path)?,
        };
        let i = |k: &str| int(h.get_i32(k).map(i64::from));
        let n = |k: &str| number(h.get_f64(k));
        let x = i64::from(h.get_i32("XBAYROFF").unwrap_or(0));
        let y = i64::from(h.get_i32("YBAYROFF").unwrap_or(0));
        (vec![i("BITPIX"), n("BZERO"), n("BSCALE"), n("DATAMAX"), n("EGAIN"), int(Some(x)), int(Some(y))], (x, y))
    };
    Ok(OldRead { xisf, values, bayer: info.bayer_pattern, offsets })
}
