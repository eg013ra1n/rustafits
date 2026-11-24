//! A FITS / XISF header read by Athenaeum's catalog rules (starloom phase 2, ruling P2-R1).
//!
//! Athenaeum's catalog reads every header through one set of rules: the scanner's FITS card state machine, its XISF
//! `<FITSKeyword>` extraction, one key map per format, and a resolver that tries a field's keywords in order. This
//! module is a port of those rules from `athenaeum-core` at ff7e0ca3 (`fits_parser/header_keys.rs`,
//! `fits_parser/fits_header_reader.rs`, the header-text and frame-building lines of `fits_parser/mod.rs`, and the
//! default chains of `header_fallbacks`), with their unit tests. Athenaeum's user-configured header rules stay
//! Athenaeum's: only the default chains are here.
//!
//! - [`header_text`]: the text Athenaeum stores in its catalog for a file. A gzipped, tile-compressed or
//!   extension-HDU FITS file gives the text of its logical header (the image's own header: structural cards first,
//!   `ZBLANK` dropped, `ZPCOUNT` / `ZGCOUNT` / `ZTENSION` renamed, `ZHECKSUM` / `ZDATASUM` mapped, the table's own
//!   cards dropped), the header the pixel reader decodes the image by.
//! - [`header_key_map`] reads such a text; [`read_header_key_map`] builds the same map from the file's bytes, as the
//!   scanner does. [`HeaderKeyMap::number`] and [`HeaderKeyMap::integer`] read a keyword that is not one of the
//!   twelve fields (`BITPIX`, `BZERO`, `DATAMAX`, `EGAIN`, `XBAYROFF`, …) by the catalog's number and integer rules.
//! - [`resolve_default`] reads one of the twelve catalog fields ([`HeaderField`]) through its default keyword chain.
//! - [`read_frame_facts`]: the [`FrameFacts`] of a file, every field read as the catalog's scanner reads it.
//!
//! The module knows no Athenaeum card (`ATH_*`) and no card of the stacking kernel's.
//!
//! Every header value the crate reads for a decision goes through this module: the FITS and XISF readers' Bayer
//! offsets (`formats::fits`, `formats::xisf`), the header facts of `star_measure::raw`, and the `rustafits info`
//! `--facts` lines and key-map columns. What does not: the decode's own structural card parse (`ImageHeader`:
//! `BITPIX`, `BZERO`, `BSCALE`, `BAYERPAT`, `ROWORDER`, the layout), and the `ImageHeader` keyword list that `info
//! --all` / `--json` print and that the built-in and `-k` columns of `info` read.
//!
//! A typed value is read by the catalog's reader: a FITS value through `FitsHeader` on the key as written (the first
//! card winning), an XISF value through the catalog's key map, an offset as the catalog's column reads the name. The
//! Bayer offsets are read exactly as the catalog's `frames.xbayroff` / `ybayroff` columns are (FITS `get_i32` on the
//! key as written, the first card winning; XISF the scanner's `parse_offset` on the element named as written, the
//! last one winning), so the render, the Analysis and the catalog agree.
//!
//! A routed FITS file (gzip, tile-compressed, or an empty primary with the image in an extension) is read by its
//! logical header, the header the pixel reader decodes. The catalog does not do that for a `.fits` with an image
//! extension: it stores and reads the primary header, so for that class a fact here and the catalog's column can
//! differ.

#![deny(missing_docs)]

use std::path::Path;

use anyhow::Result;

mod chains;
mod facts;
mod fits_header_reader;
mod keys;
#[cfg(test)]
mod tests;
mod text;

pub use chains::{resolve_default, HeaderField, HeaderValue};
pub(crate) use facts::{fits_bayer_offsets_of_blocks, xisf_bayer_offsets};
pub(crate) use text::fits_header;
pub use facts::{read_frame_facts, FrameFacts};
pub use keys::{header_key_map, normalize_keyword, HeaderKeyMap};

/// The container format a header was read from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize), serde(rename_all = "camelCase"))]
#[non_exhaustive]
pub enum HeaderFormat {
    /// FITS: 80-column cards (plain, gzipped, tile-compressed or in an extension HDU).
    #[default]
    Fits,
    /// XISF: `<FITSKeyword>` elements of the XML header.
    Xisf,
}

/// The header text of the file at `path` and its format: the text Athenaeum's catalog stores for the file (for a
/// routed FITS file, the text of its logical header; see the module docs).
///
/// - FITS: every card through `END`, one per line (Athenaeum's `extract_fits_header`). For a routed file (gzip, an
///   empty primary with the image in an extension, a tile-compressed image), the cards of the image's logical header.
/// - XISF: the `XISF FITS Keywords:` banner and one `NAME = value` line per `<FITSKeyword/>` element (Athenaeum's
///   `extract_xisf_header`), or, when the header has none, the line naming the file that Athenaeum writes instead.
///
/// The format is XISF for a `.xisf` file or one with the `XISF0100` signature, FITS otherwise.
pub fn header_text(path: &Path) -> Result<(HeaderFormat, String)> {
    let format = text::detect_format(path);
    let text = match format {
        HeaderFormat::Fits => text::extract_fits_header(path)?,
        HeaderFormat::Xisf => text::extract_xisf_header(path)?,
    };
    Ok((format, text))
}

/// The key map of the file at `path`, built from its header bytes as Athenaeum's scanner builds it (FITS: the cards,
/// the first of a keyword winning; XISF: the `<FITSKeyword/>` elements, the last winning). It answers as
/// `header_key_map` over [`header_text`] does, except where the stored text cannot carry what the bytes say: an XISF
/// value holding a line break, an XISF name holding `" = "`, or a header with no keyword (Athenaeum's placeholder
/// text).
pub fn read_header_key_map(path: &Path) -> Result<HeaderKeyMap> {
    match text::detect_format(path) {
        HeaderFormat::Fits => Ok(HeaderKeyMap::from_fits_header(&text::fits_header(path)?).with_format(HeaderFormat::Fits)),
        HeaderFormat::Xisf => {
            let content = text::read_xisf_header(path)?;
            let pairs = text::xisf_keyword_pairs(&content, Some(path))?;
            Ok(HeaderKeyMap::from_last_wins(pairs).with_format(HeaderFormat::Xisf))
        }
    }
}
