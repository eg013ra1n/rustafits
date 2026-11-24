//! The header text of a file, as Athenaeum's catalog stores it, and the XISF `<FITSKeyword>` extraction behind it:
//! Athenaeum's `athenaeum-core/src/fits_parser/mod.rs` at ff7e0ca3 (`read_xisf_header` with its prelude check,
//! `xisf_xml_slice`, `xisf_keyword_element`, `xisf_keyword_pairs`, `xisf_header_text`, `extract_fits_header`,
//! `extract_xisf_header`), bodies unchanged apart from paths. One addition: a gzipped, tile-compressed or
//! extension-HDU FITS file is read by its logical header, the header blocks `formats::fits::routed_image_header`
//! returns (PR #5's rule; Athenaeum's catalog does not index such files).

use std::io::{BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result};

use super::fits_header_reader::FitsHeader;
use super::keys::XISF_DUMP_BANNER;
use super::HeaderFormat;

/// The format of `path` as this crate's readers decide it (`formats::read_header`): XISF by its extension or its
/// `XISF0100` signature, FITS otherwise.
pub(crate) fn detect_format(path: &Path) -> HeaderFormat {
    if path.extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case("xisf")) {
        return HeaderFormat::Xisf;
    }
    let mut sig = [0u8; 8];
    let signed = std::fs::File::open(path).and_then(|mut f| f.read_exact(&mut sig)).is_ok();
    if signed && &sig == b"XISF0100" {
        HeaderFormat::Xisf
    } else {
        HeaderFormat::Fits
    }
}

/// The FITS header a file is read by: the image's logical header when PR #5 routes the file (gzip, an empty primary
/// with the image in an extension, a tile-compressed image), else the primary header, read as Athenaeum's catalog
/// reads it (`FitsHeader::from_path`). A routing failure on a file that is not gzipped (a primary the pixel reader
/// refuses, an extension it cannot find) falls back to the primary header, which is what the catalog stores for it;
/// the failure is logged. A gzipped file has no other header to fall back to.
pub(crate) fn fits_header(path: &Path) -> Result<FitsHeader> {
    match crate::formats::fits::routed_image_header(path) {
        Ok(Some(blocks)) => FitsHeader::from_bytes(&blocks)
            .with_context(|| format!("logical header of {}", path.display())),
        Ok(None) => FitsHeader::from_path(path),
        Err(e) => {
            let mut magic = [0u8; 2];
            let gzip = std::fs::File::open(path)
                .and_then(|mut f| f.read_exact(&mut magic))
                .is_ok()
                && fitsio_pure::gzip::is_gzip(&magic);
            if gzip {
                return Err(e.context(format!("gzipped FITS header of {}", path.display())));
            }
            tracing::warn!(
                path = %path.display(),
                error = %format!("{e:#}"),
                "image header routing failed; reading the primary header as the catalog does"
            );
            FitsHeader::from_path(path)
        }
    }
}

/// Length of the XISF XML header declared by the 16-byte prelude (8-byte signature, u32 LE
/// header length, 4 reserved bytes). Shared by the disk reader and the zip byte entry point.
fn xisf_prelude_header_len(prelude: &[u8; 16]) -> Result<usize> {
    if &prelude[..8] != b"XISF0100" {
        anyhow::bail!("Not a valid XISF file: bad signature");
    }
    let header_length = u32::from_le_bytes([prelude[8], prelude[9], prelude[10], prelude[11]]) as usize;
    // Safety cap at 64MB
    if header_length > 64 * 1024 * 1024 {
        anyhow::bail!("XISF header too large: {} bytes", header_length);
    }
    Ok(header_length)
}

/// Read the XISF XML header bytes from a file.
/// Parses the 16-byte binary header (8-byte signature + 4-byte header length + 4 reserved)
/// to determine exactly how many bytes to read.
pub(crate) fn read_xisf_header(path: &Path) -> Result<Vec<u8>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("Failed to open XISF file: {}", path.display()))?;
    let mut buf_reader = BufReader::new(file);

    let mut prelude = [0u8; 16];
    buf_reader.read_exact(&mut prelude)?;
    let header_length = xisf_prelude_header_len(&prelude)?;

    let mut content = vec![0u8; header_length];
    buf_reader.read_exact(&mut content)?;
    Ok(content)
}

/// Extract full FITS header as text (all raw cards, not just a keyword whitelist). For a routed file (gzip,
/// extension HDU, tile compression), the cards of its logical header.
pub(crate) fn extract_fits_header(path: &Path) -> Result<String> {
    tracing::debug!(path = %path.display(), "reading FITS header");
    let header = fits_header(path)?;
    Ok(header.to_header_text())
}

/// Extract full XISF header as text
pub(crate) fn extract_xisf_header(path: &Path) -> Result<String> {
    tracing::debug!(path = %path.display(), "opening XISF for header");
    let content = read_xisf_header(path)?;
    xisf_header_text(&content, path)
}

/// The `<?xml`…`</xisf>` slice of an XISF header block.
pub(crate) fn xisf_xml_slice(content: &[u8]) -> Result<&[u8]> {
    let xml_start = content.windows(5)
        .position(|w| w == b"<?xml")
        .ok_or_else(|| anyhow::anyhow!("No XML header found in XISF file"))?;

    let xml_end = content.windows(7)
        .skip(xml_start)
        .position(|w| w == b"</xisf>")
        .ok_or_else(|| anyhow::anyhow!("No closing </xisf> tag found"))?;
    let xml_end = xml_start + xml_end + 7; // +7 for the length of "</xisf>"

    Ok(&content[xml_start..xml_end])
}

/// `(name, value)` of one `<FITSKeyword/>` element, the value's outer `'…'` stripped; `None`
/// when the element has no name. The one reading of the element, shared by the frame builder,
/// the header text and the stored-text key map.
pub(crate) fn xisf_keyword_element(e: &quick_xml::events::BytesStart) -> Option<(String, String)> {
    let mut name = String::new();
    let mut value = String::new();

    for attr in e.attributes().flatten() {
        match attr.key.as_ref() {
            b"name" => {
                name = String::from_utf8_lossy(&attr.value).to_string();
            }
            b"value" => {
                value = String::from_utf8_lossy(&attr.value).to_string();
                // Remove quotes if present
                if value.len() >= 2 && value.starts_with('\'') && value.ends_with('\'') {
                    value = value[1..value.len()-1].to_string();
                }
            }
            _ => {}
        }
    }

    if name.is_empty() { None } else { Some((name, value)) }
}

/// Every `FITSKeyword` element of an XISF XML block, in document order (duplicates kept).
/// `label` only names the file in logs (`None` for text read back from the catalog, which has no
/// file). A malformed tail yields the keywords read so far.
pub(crate) fn xisf_keyword_pairs(content: &[u8], label: Option<&Path>) -> Result<Vec<(String, String)>> {
    use quick_xml::events::Event;
    use quick_xml::Reader;

    let xml_content = xisf_xml_slice(content)?;
    let xml_str = String::from_utf8_lossy(xml_content);

    // Parse XML to extract FITSKeyword elements
    let mut reader = Reader::from_str(&xml_str);
    reader.config_mut().trim_text(true);
    // Permissive dangling-`&` retained deliberately (0.36-era behavior; XISF
    // headers are machine-written but a malformed one should still render —
    // cycle decision 2026-07-29).
    reader.config_mut().allow_dangling_amp = true;

    let mut pairs = Vec::new();
    let mut buf = Vec::new();
    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Empty(ref e)) if e.name().as_ref() == b"FITSKeyword" => {
                if let Some(pair) = xisf_keyword_element(e) {
                    pairs.push(pair);
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => {
                if let Some(label) = label {
                    tracing::error!(
                        path = %label.display(),
                        position = reader.buffer_position(),
                        error = %e,
                        "error parsing XISF XML header; using partial result"
                    );
                } else {
                    tracing::warn!(
                        position = reader.buffer_position(),
                        error = %e,
                        "error parsing stored XISF XML header; using partial result"
                    );
                }
                break;
            }
            _ => {}
        }
        buf.clear();
    }
    Ok(pairs)
}

/// Header text of an XISF XML block (`NAME = value` lines). `path` only labels logs and the
/// no-keywords fallback text.
pub(crate) fn xisf_header_text(content: &[u8], path: &Path) -> Result<String> {
    let pairs = xisf_keyword_pairs(content, Some(path))?;

    let mut header_text = String::from(XISF_DUMP_BANNER);
    for (name, value) in pairs {
        header_text.push_str(&format!("{} = {}\n", name, value));
    }

    // If header is still minimal, record that we tried
    if header_text.lines().count() <= 3 {
        header_text = format!("Header extracted from: {}\n(No FITS keywords found in XISF)", path.display());
    }

    Ok(header_text)
}
