//! Athenaeum's header-value rules, ported so that no Analysis number moves on any input (port plan phase 4, ruling
//! R-P4-C1, 2026-11-05).
//!
//! Before phase 4, Athenaeum's Analysis tab read the header facts its measurement needs (`DATAMAX`, `EGAIN`, and for
//! FITS `BITPIX`, `BZERO`, `BSCALE`; for XISF the Bayer offsets `XBAYROFF` / `YBAYROFF`) with its own catalog rules,
//! not with this crate's `ImageHeader`. Those rules differ from `ImageHeader`'s on unusual or malformed values (an
//! empty value before a later one, a `' / '` comment inside the value, an unclosed quote, a character reference, a
//! card with no `=` in column 9, `HIERARCH`, `CONTINUE`). No corpus file has any, but a new path copies the proven
//! path, so the rules are copied here and the two readers and `star_measure::raw` use them:
//!
//! - **FITS**: [`FitsHeader`], Athenaeum's `fits_parser/fits_header_reader.rs` at e6344ddd, file and tests unchanged
//!   (see its own doc). `star_measure::raw` reads `BITPIX`, `BZERO`, `BSCALE`, `DATAMAX` and `EGAIN` through its
//!   `get_i32` / `get_f64`. The FITS reader folds `XBAYROFF` / `YBAYROFF` into the pattern itself (its own card
//!   rules, unchanged), as it did for Athenaeum's Analysis.
//! - **XISF**: Athenaeum turns the XML into text first, then reads the text. [`xisf_keyword_text`] is the
//!   `<FITSKeyword>` loop of `fits_parser::extract_xisf_header` (`NAME = value` lines: Empty elements only, the raw
//!   attribute value, one pair of quotes stripped); [`parse_stored_xisf_keys`] is `stored_header::parse_stored_header_keys`
//!   for XISF with its `parse_xisf_xml_text` and `parse_keyword_eq_text`; [`num`] and [`offset`] are `analysis/input.rs`'s.
//!   The XISF reader (`formats::xisf`, the Bayer offsets) and `star_measure::raw` (`DATAMAX`, `EGAIN`) read their
//!   values through [`xisf_stored_keys`].
//!
//! The functions are Athenaeum's code. What is not copied: Athenaeum's file reading (the callers hand over the XML
//! string their reader already read), and the "no keywords found" fallback line of `xisf_header_text`, which wrote the
//! file's path into the text (the map read from the banner alone is empty). That is the one residual: for a header with
//! no `<FITSKeyword>`, Athenaeum read a key from a file path that spells one (a newline then `EGAIN = 2`, or a
//! `<FITSKeyword …>` element) and the copy reads nothing. So a number never differs, except through a key spelled by
//! the file's path and not by its header. Differences in the code:
//!
//! - a header text is cut to the `<?xml` … `</xisf>` section of the string as Athenaeum cuts its bytes, but a string
//!   without either marker is read whole, where Athenaeum's reader returned an error;
//! - `value="'"` (one lone quote) is not a pair of quotes: Athenaeum's loop sliced `value[1..0]` and panicked; here
//!   the value stays `'`, which the line rule then skips as an unclosed quote;
//! - the parse-abort log line carries the file's `path` field when the caller has one ([`xisf_stored_keys`] takes an
//!   optional path: `read_xisf_stored_keys` passes it, the XISF reader's Bayer-offset fold has none to hand).
//!
//! The pinned inputs are in the tests: `xisf_rule_pins_athenaeums_results` and `fits_rule_pins_athenaeums_results`
//! list what the rules return on each unusual shape, and the XISF reader's and `star_measure::raw`'s own tests pin it
//! end to end.

// Athenaeum's code as copied: the clippy style lint its attribute loop trips is allowed, never rewritten.
#![allow(clippy::manual_flatten)]

use std::collections::HashMap;
use std::path::Path;

mod fits_header_reader;

pub(crate) use fits_header_reader::FitsHeader;

/// The `NAME = value` text of an XISF header's `<FITSKeyword>` elements: Athenaeum's `fits_parser::xisf_header_text`
/// loop (banner, then one line per Empty `<FITSKeyword>` with a non-empty `name`, the raw `value` attribute with one
/// pair of quotes removed). `xml` is the header string the XISF reader read; `path` is the file it came from, if the
/// caller knows it, and goes on the log line of a parse abort as Athenaeum's does.
pub(crate) fn xisf_keyword_text(xml: &str, path: Option<&Path>) -> String {
    use quick_xml::events::Event;
    use quick_xml::Reader;

    // The `<?xml` ... `</xisf>` section, as Athenaeum cuts its bytes.
    let xml_start = xml.find("<?xml").unwrap_or(0);
    let xml_end = xml[xml_start..]
        .find("</xisf>")
        .map_or(xml.len(), |p| xml_start + p + "</xisf>".len());
    let xml_str = &xml[xml_start..xml_end];

    // Parse XML to extract FITSKeyword elements
    let mut reader = Reader::from_str(xml_str);
    reader.config_mut().trim_text(true);
    // Permissive dangling-`&` retained deliberately (0.36-era behavior; XISF
    // headers are machine-written but a malformed one should still render —
    // cycle decision 2026-07-29).
    reader.config_mut().allow_dangling_amp = true;

    let mut header_text = String::new();
    header_text.push_str("XISF FITS Keywords:\n");
    header_text.push_str("==================\n\n");

    let mut buf = Vec::new();
    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Empty(ref e)) if e.name().as_ref() == b"FITSKeyword" => {
                let mut name = String::new();
                let mut value = String::new();

                for attr in e.attributes() {
                    if let Ok(attr) = attr {
                        match attr.key.as_ref() {
                            b"name" => {
                                name = String::from_utf8_lossy(&attr.value).to_string();
                            }
                            b"value" => {
                                value = String::from_utf8_lossy(&attr.value).to_string();
                                // Remove quotes if present (a lone `'` is not a pair: Athenaeum's slice panicked)
                                if value.len() >= 2 && value.starts_with('\'') && value.ends_with('\'') {
                                    value = value[1..value.len()-1].to_string();
                                }
                            }
                            _ => {}
                        }
                    }
                }

                if !name.is_empty() {
                    header_text.push_str(&format!("{} = {}\n", name, value));
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => {
                tracing::error!(
                    path = path.map(|p| tracing::field::display(p.display())),
                    position = reader.buffer_position(),
                    error = %e,
                    "error parsing XISF XML header; using partial result"
                );
                break;
            }
            _ => {}
        }
        buf.clear();
    }

    header_text
}

/// The stored-header keys of an XISF header text: Athenaeum's `stored_header::parse_stored_header_keys` for
/// `FileFormat::XISF`. The text is read as XML first (a stored header may be the raw `<?xml` block); when that yields
/// nothing it is read as `KEY = value` lines, which is what [`xisf_keyword_text`] produces.
pub(crate) fn parse_stored_xisf_keys(header_text: &str) -> HashMap<String, String> {
    let primary = parse_xisf_xml_text(header_text);
    if !primary.is_empty() {
        return primary;
    }
    // Some capture tools persist the header as a plain "KEY = value" text dump
    // rather than 80-col FITS cards or XISF XML, so the format-specific parser
    // yields nothing. Fall back to a line parser.
    parse_keyword_eq_text(header_text)
}

/// The stored-header keys of the XISF header string `xml`, read the way Athenaeum reads an XISF file's header: the
/// XML turned into `NAME = value` text, then the text into upper-cased keys, the first non-empty value of a name
/// winning. `path` is the file `xml` came from, if known (see [`xisf_keyword_text`]).
pub(crate) fn xisf_stored_keys(xml: &str, path: Option<&Path>) -> HashMap<String, String> {
    parse_stored_xisf_keys(&xisf_keyword_text(xml, path))
}

/// Parse a plain "KEY = value" header dump (one keyword per line, banner /
/// separator / COMMENT / HISTORY / blank lines skipped). Quoted strings have
/// their quotes stripped; an unquoted value's trailing " / comment" is
/// dropped. Keys are upper-cased to match the other parsers.
fn parse_keyword_eq_text(text: &str) -> HashMap<String, String> {
    let mut result: HashMap<String, String> = HashMap::new();
    for raw in text.lines() {
        let line = raw.trim();
        let Some((k, v)) = line.split_once('=') else { continue };
        let key = k.trim().to_uppercase();
        if key.is_empty()
            || key.contains(' ')
            || key == "END"
            || key == "COMMENT"
            || key == "HISTORY"
        {
            continue;
        }
        let v = v.trim();
        let value = if let Some(s) = v.strip_prefix('\'') {
            match s.find('\'') {
                Some(end) => s[..end].trim_end().to_string(),
                None => continue,
            }
        } else {
            // Drop a FITS-style trailing comment (" / ...") without breaking
            // space-separated values like "2 41 04.824".
            v.split(" / ").next().unwrap_or("").trim().to_string()
        };
        if value.is_empty() {
            continue;
        }
        result.entry(key).or_insert(value);
    }
    result
}

/// Parse XISF XML text — `FITSKeyword name="..." value="..."` elements,
/// with quotes stripped from values to match how `parse_xisf` already
/// records them.
fn parse_xisf_xml_text(text: &str) -> HashMap<String, String> {
    use quick_xml::events::Event;
    use quick_xml::Reader;

    let mut result: HashMap<String, String> = HashMap::new();
    let mut reader = Reader::from_str(text);
    reader.config_mut().trim_text(true);
    // Permissive dangling-`&` retained deliberately (0.36-era behavior; XISF
    // headers are machine-written but a malformed one should still render —
    // cycle decision 2026-07-29).
    reader.config_mut().allow_dangling_amp = true;
    let mut buf = Vec::new();
    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Empty(ref e)) | Ok(Event::Start(ref e))
                if e.name().as_ref() == b"FITSKeyword" =>
            {
                let mut name = String::new();
                let mut value = String::new();
                for attr in e.attributes().flatten() {
                    match attr.key.as_ref() {
                        b"name" => name = String::from_utf8_lossy(&attr.value).to_string(),
                        b"value" => {
                            let mut v = String::from_utf8_lossy(&attr.value).to_string();
                            if v.starts_with('\'') && v.ends_with('\'') && v.len() >= 2 {
                                v = v[1..v.len() - 1].trim().to_string();
                            }
                            value = v;
                        }
                        _ => {}
                    }
                }
                if !name.is_empty() {
                    result.entry(name.to_uppercase()).or_insert(value);
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => {
                tracing::warn!(error = %e, "stored XISF header parse aborted");
                break;
            }
            _ => {}
        }
        buf.clear();
    }
    result
}

/// A stored-header value as a number: quotes and spaces trimmed, then `f64` (Athenaeum's `analysis/input.rs` `num`).
pub(crate) fn num(map: &HashMap<String, String>, key: &str) -> Option<f64> {
    map.get(key)
        .and_then(|v| v.trim().trim_matches('\'').trim().parse::<f64>().ok())
}

/// A header offset as an integer; absent or unreadable is 0 (Athenaeum's `analysis/input.rs` `offset`).
pub(crate) fn offset(v: Option<f64>) -> i64 {
    v.filter(|x| x.is_finite()).map_or(0, |x| x.round() as i64)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Athenaeum's own tests of the ported functions (stored_header.rs), bodies unchanged ----

    #[test]
    fn parses_xisf_xml_fitskeyword_elements() {
        let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<xisf>
  <FITSKeyword name="OBJECT" value="'M42'" comment="" />
  <FITSKeyword name="EXPTIME" value="120.0" />
  <FITSKeyword name="CCD-TEMP" value="-10.0" />
  <FITSKeyword name="XBINNING" value="1" />
  <FITSKeyword name="YBINNING" value="1" />
</xisf>"#;
        let keys = parse_xisf_xml_text(xml);
        assert_eq!(keys.get("OBJECT"), Some(&"M42".to_string()));
        assert_eq!(keys.get("EXPTIME"), Some(&"120.0".to_string()));
        assert_eq!(keys.get("CCD-TEMP"), Some(&"-10.0".to_string()));
        assert_eq!(keys.get("XBINNING"), Some(&"1".to_string()));
    }

    /// quick-xml 0.41 rejects a lone `&` (one not opening a character or
    /// entity reference) in element text by default, where 0.36/0.37 read
    /// straight through. We opt back into the permissive behavior, so a stored
    /// header carrying one must still yield its keywords rather than silently
    /// losing the whole snapshot.
    #[test]
    fn dangling_ampersand_in_xisf_text_still_yields_keywords() {
        let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<xisf>
  <Property id="note">Dark & Flat calibration</Property>
  <FITSKeyword name="OBJECT" value="'M42'" comment="" />
  <FITSKeyword name="EXPTIME" value="120.0" />
</xisf>"#;
        let keys = parse_xisf_xml_text(xml);

        // Not just "no error": the FITSKeyword elements *after* the dangling
        // `&` must still have been reached and read correctly.
        assert_eq!(
            keys.get("OBJECT"),
            Some(&"M42".to_string()),
            "keywords after a dangling `&` must still parse"
        );
        assert_eq!(keys.get("EXPTIME"), Some(&"120.0".to_string()));
    }

    /// Athenaeum's test of the plain keyword dump (a capture tool's text header), without its snapshot half
    /// (`snapshot_from_keys` is not ported) and without the tool's name in the dump.
    #[test]
    fn parses_a_plain_keyword_dump() {
        // A capture tool persists the XISF header as a plain "KEY = value" text dump
        // (not XISF XML, not 80-col FITS cards). parse_xisf_xml_text returns
        // empty for this; the dispatcher must fall back to the line parser.
        let dump = "XISF FITS Keywords:\n==================\n\nSIMPLE = T\nBITPIX = 16\nNAXIS = 2\nCOMMENT = \nFOCALLEN = 140.616\nXPIXSZ = 2.4\nOBJCTRA = 2 41 04.824\nOBJCTDEC = +62 56 44.95\nCROTA2 = 4.383973760472021\nCREATOR = a capture tool\n";
        // XML parser must yield nothing for this format.
        assert!(parse_xisf_xml_text(dump).is_empty());
        // Dispatcher must recover the keys via the fallback.
        let keys = parse_stored_xisf_keys(dump);
        assert_eq!(keys.get("FOCALLEN"), Some(&"140.616".to_string()));
        assert_eq!(keys.get("OBJCTRA"), Some(&"2 41 04.824".to_string()));
        assert_eq!(keys.get("OBJCTDEC"), Some(&"+62 56 44.95".to_string()));
    }

    // ---- the pinned inputs (R-P4-C1): Athenaeum's results on the shapes `ImageHeader`'s rules read otherwise ----

    /// An XISF header holding `children` (the `<FITSKeyword>`s), read as Athenaeum reads an XISF file's header.
    fn xisf_keys(children: &str) -> HashMap<String, String> {
        xisf_stored_keys(
            &format!(
                r#"<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0"><Image geometry="2:2:1">{children}</Image></xisf>"#
            ),
            None,
        )
    }

    fn kw(name: &str, value: &str) -> String {
        format!(r#"<FITSKeyword name="{name}" value="{value}" comment=""/>"#)
    }

    /// What Athenaeum's rule returns for `XBAYROFF` (and so for any key) on the shapes the rule `ImageHeader::get` +
    /// `num` reads differently. Every row is Athenaeum's result; the right-hand comment says what the first keyword,
    /// XML-normalized, would read instead.
    #[test]
    fn xisf_rule_pins_athenaeums_results() {
        let get = |children: String| xisf_keys(&children).get("XBAYROFF").cloned();
        // An empty value is skipped, so a later keyword of the name wins (first-keyword rule: "").
        assert_eq!(get(format!("{}{}", kw("XBAYROFF", ""), kw("XBAYROFF", "1"))), Some("1".into()));
        assert_eq!(get(format!("{}{}", kw("XBAYROFF", "''"), kw("XBAYROFF", "1"))), Some("1".into()));
        assert_eq!(
            get(format!(r#"<FITSKeyword name="XBAYROFF"/>{}"#, kw("XBAYROFF", "1"))),
            Some("1".into()),
            "no value attribute"
        );
        // A quoted value runs to the next quote; an unquoted one is cut at " / " (first-keyword rule: "1 / off").
        assert_eq!(get(kw("XBAYROFF", "1 / off")), Some("1".into()));
        assert_eq!(get(kw("XBAYROFF", "'1' / off")), Some("1".into()));
        // An unclosed quote drops the line (first-keyword rule: "'1").
        assert_eq!(get(kw("XBAYROFF", "'1")), None);
        // The name is trimmed and upper-cased by the line rule (first-keyword rule: no match).
        assert_eq!(get(kw(" XBAYROFF", "1")), Some("1".into()));
        assert_eq!(get(kw("xbayroff", "1")), Some("1".into()));
        // The raw attribute is read, not the XML-normalized one (first-keyword rule: "1").
        assert_eq!(get(kw("XBAYROFF", "&#49;")), Some("&#49;".into()));
        // The first non-empty value wins.
        assert_eq!(get(format!("{}{}", kw("XBAYROFF", "0"), kw("XBAYROFF", "1"))), Some("0".into()));
        // Empty elements only: a start / end pair is not read (first-keyword rule: read).
        assert_eq!(get(r#"<FITSKeyword name="XBAYROFF" value="1"></FITSKeyword>"#.into()), None);
        // A `<Property>` is not a keyword (`ImageHeader::get` matches `<Property id=…>` too).
        assert_eq!(get(r#"<Property id="XBAYROFF" type="Int32" value="1"/>"#.into()), None);
        // A lone quote is no pair of quotes (Athenaeum's slice panicked here); the line rule skips it.
        assert_eq!(get(kw("XBAYROFF", "'")), None);
        // One pair of quotes is stripped from the raw attribute before the line rule reads the line, so `'1'0'` is the
        // line `XBAYROFF = 1'0`: unreadable, so the offset is 0 (without the strip the line rule reads `'1'` as 1).
        assert_eq!(get(kw("XBAYROFF", "'1'0'")), Some("1'0".into()));
        let keys = xisf_keys(&kw("XBAYROFF", "'1'0'"));
        assert_eq!(num(&keys, "XBAYROFF"), None);
        assert_eq!(offset(num(&keys, "XBAYROFF")), 0);
        // The header is cut at the first `<?xml`: a keyword in front of it is not read (read whole: 1).
        let before_the_cut = format!(
            r#"{}<?xml version="1.0"?><xisf><Image geometry="2:2:1"/></xisf>"#,
            kw("XBAYROFF", "1")
        );
        assert_eq!(xisf_stored_keys(&before_the_cut, None).get("XBAYROFF"), None);
    }

    #[test]
    fn xisf_num_and_offset_read_as_athenaeum_reads_them() {
        let keys = xisf_keys(&format!(
            "{}{}{}{}{}",
            kw("A", "1.5"),
            kw("B", "'1.5'"),
            kw("C", "0.6"),
            kw("D", "-1"),
            kw("E", "nan")
        ));
        assert_eq!(num(&keys, "A"), Some(1.5));
        assert_eq!(num(&keys, "B"), Some(1.5));
        assert_eq!(num(&keys, "Z"), None);
        assert_eq!(offset(num(&keys, "C")), 1);
        assert_eq!(offset(num(&keys, "D")), -1);
        assert_eq!(offset(num(&keys, "E")), 0, "non-finite is 0");
        assert_eq!(offset(num(&keys, "Z")), 0, "absent is 0");
    }

    #[test]
    fn an_xisf_without_keywords_reads_no_keys() {
        assert!(xisf_keys("").is_empty());
        assert!(xisf_stored_keys("", None).is_empty());
        assert!(xisf_stored_keys("not xml at all", None).is_empty());
    }

    /// The `<?xml` … `</xisf>` cut: a keyword after the closing tag is not read, as in Athenaeum's reader.
    #[test]
    fn a_keyword_after_the_closing_tag_is_not_read() {
        let xml = format!(
            r#"junk<?xml version="1.0"?><xisf>{}</xisf>{}"#,
            kw("EGAIN", "1.5"),
            kw("DATAMAX", "9")
        );
        let keys = xisf_stored_keys(&xml, None);
        assert_eq!(keys.get("EGAIN").map(String::as_str), Some("1.5"));
        assert_eq!(keys.get("DATAMAX"), None);
    }

    /// The keyword pass's parse-abort ERROR line names the file when the caller has it (Athenaeum's `path` field),
    /// and has no `path` field when it has not. A `</xisf>` inside a comment ends the cut section in the comment.
    #[test]
    fn the_parse_abort_log_line_carries_the_path_when_given() {
        let xml = format!(
            r#"<?xml version="1.0"?><xisf><!-- </xisf> -->{}</xisf>"#,
            kw("EGAIN", "2")
        );
        let path = Path::new("/data/a_frame.xisf");
        let mut keys = HashMap::new();
        let with = crate::test_support::captured_logs(|| keys = xisf_stored_keys(&xml, Some(path)));
        assert!(with.contains("error parsing XISF XML header; using partial result"), "{with}");
        assert!(with.contains("path=/data/a_frame.xisf"), "{with}");
        assert_eq!(keys.get("EGAIN"), None, "the keyword is past the cut: the partial result is empty");
        let without = crate::test_support::captured_logs(|| keys = xisf_stored_keys(&xml, None));
        assert!(without.contains("error parsing XISF XML header; using partial result"), "{without}");
        assert!(!without.contains("path="), "{without}");
    }

    /// A header whose keyword text itself holds a `<FITSKeyword …/>` element reads as XML first (Athenaeum's
    /// dispatcher): the element inside the text is the key set, the dump lines are not.
    #[test]
    fn a_keyword_value_holding_an_element_is_read_as_xml_first() {
        let keys = xisf_keys(&kw("Z", "<FITSKeyword name='XBAYROFF' value='1'/>"));
        assert_eq!(keys.get("XBAYROFF").map(String::as_str), Some("1"));
        assert_eq!(keys.get("Z"), None);
    }

    /// 80-column FITS cards, `END`, padded to a 2880-byte block.
    fn fits_header(cards: &[&str]) -> FitsHeader {
        let mut block = Vec::new();
        for c in cards {
            block.extend_from_slice(format!("{c:<80}").as_bytes());
        }
        block.extend_from_slice(format!("{:<80}", "END").as_bytes());
        block.resize(block.len().div_ceil(2880) * 2880, b' ');
        FitsHeader::parse_header(&block).unwrap()
    }

    /// What Athenaeum's `FitsHeader` returns for `EGAIN` on the card shapes that `ImageHeader`'s card parser reads
    /// otherwise. The right-hand comment says what the other rule reads.
    #[test]
    fn fits_rule_pins_athenaeums_results() {
        let egain = |cards: &[&str]| fits_header(cards).get_f64("EGAIN");
        // No `= ` in columns 9-10: the card is skipped (other rule: 1.5).
        assert_eq!(egain(&["EGAIN    1.5"]), None);
        assert_eq!(egain(&["EGAIN   =1.5"]), None);
        // A value-less card is skipped, so the real card after it counts (other rule: the first card, no value).
        assert_eq!(egain(&["EGAIN", "EGAIN   = 1.5"]), Some(1.5));
        // HIERARCH cards are read under their long key (other rule: not matched).
        assert_eq!(egain(&["HIERARCH EGAIN = 1.5"]), Some(1.5));
        // Blanks in front of the key are trimmed (other rule: not matched).
        assert_eq!(egain(&[" EGAIN  = 1.5"]), Some(1.5));
        // A string ending in `&` is a continuation: the marker is dropped (other rule: "1.5&", unreadable)...
        assert_eq!(egain(&["EGAIN   = '1.5&'"]), Some(1.5));
        // ... and the CONTINUE cards after it are joined on (other rule: "1.").
        assert_eq!(egain(&["EGAIN   = '1.&'", "CONTINUE  '5'"]), Some(1.5));
        // Fortran exponents; the first card of a keyword wins; a comment after `/` is cut.
        assert_eq!(egain(&["EGAIN   = 1.5D0 / e-/ADU"]), Some(1.5));
        assert_eq!(egain(&["EGAIN   = 1.5", "EGAIN   = 2.5"]), Some(1.5));
        // A quoted number is unquoted by the card rule.
        assert_eq!(egain(&["EGAIN   = '1.5'"]), Some(1.5));
        assert_eq!(fits_header(&["BITPIX  = 16.0"]).get_i32("BITPIX"), Some(16));
    }
}
