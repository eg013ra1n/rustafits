//! One key reader per header format: Athenaeum's `athenaeum-core/src/fits_parser/header_keys.rs` at ff7e0ca3
//! (starloom phase 2, ruling P2-R1), with its unit tests.
//!
//! Every place in Athenaeum that turns a header into catalog fields (the FITS and XISF frame builders, the metadata
//! pane's snapshot, the override recomputation, Apply) reads a [`HeaderKeyMap`]. The map is built by the scanner's own
//! reading of each format, so a stored header answers exactly as the file did when it was scanned:
//!
//! - FITS: the scanner's card state machine (`FitsHeader::from_header_text`), HIERARCH, doubled quotes and CONTINUE
//!   cards included.
//! - XISF: the scanner's `FITSKeyword` extraction (`xisf_keyword_pairs`), or the `NAME = value` dump Athenaeum stores
//!   (`xisf_header_text`), which this module inverts.
//!
//! Keys are normalised with [`normalize_keyword`]; values are the raw text of the format (the readers trim).
//!
//! What differs from Athenaeum's file: `stored_header_key_map` is [`header_key_map`] here and takes a
//! [`HeaderFormat`]; the map remembers that format, because [`HeaderKeyMap::integer`] reads an integer by the rule of
//! the format's own scanner line; [`HeaderKeyMap::number`] and [`HeaderKeyMap::integer`] are new (the catalog reads
//! `BITPIX`, `XBAYROFF` and the like through `FitsHeader` / the XISF keyword map, not through a key map); the map has
//! no `Default`.

use std::collections::HashMap;

use super::fits_header_reader::{parse_fits_f64, FitsHeader};
use super::text::xisf_keyword_pairs;
use super::HeaderFormat;

/// First line of the XISF keyword dump `xisf_header_text` writes: title, rule and a blank line.
pub(crate) const XISF_DUMP_BANNER: &str = "XISF FITS Keywords:\n==================\n\n";

/// A header keyword as the rules and the readers compare it: trimmed and upper-cased, so a user
/// keyword matches a HIERARCH name whatever its case (inner spaces kept). `None` when the result
/// is empty or contains `=` (such a name cannot be a keyword).
pub fn normalize_keyword(raw: &str) -> Option<String> {
    let k = raw.trim().to_uppercase();
    if k.is_empty() || k.contains('=') {
        None
    } else {
        Some(k)
    }
}

/// Normalised keyword -> raw value of one header, and the format it was read from.
#[derive(Debug, Clone)]
pub struct HeaderKeyMap {
    map: HashMap<String, String>,
    format: HeaderFormat,
}

impl HeaderKeyMap {
    /// The raw value of a keyword (the key is normalised here, so `"ccd-temp"` finds `CCD-TEMP`).
    pub fn get(&self, key: &str) -> Option<&str> {
        match self.map.get(key) {
            Some(v) => Some(v.as_str()),
            None => normalize_keyword(key).and_then(|k| self.map.get(&k)).map(String::as_str),
        }
    }

    /// True when the header holds any keyword other than `COMMENT` / `HISTORY`.
    pub fn has_keywords(&self) -> bool {
        self.map.keys().any(|k| k != "COMMENT" && k != "HISTORY")
    }

    /// The format the map was read from.
    pub fn format(&self) -> HeaderFormat {
        self.format
    }

    /// A keyword's value as a number, by the catalog's number rule (`parse_fits_f64`, which `FitsHeader::get_f64`
    /// and the field resolver use): the value trimmed, a Fortran `D` exponent read as `E`. `NaN` and infinities are
    /// returned as written; a value that is not a number, or an absent keyword, is `None`.
    pub fn number(&self, key: &str) -> Option<f64> {
        self.get(key).and_then(|v| parse_fits_f64(v.trim()))
    }

    /// A keyword's value as an integer, by the rule of the catalog's scanner line for the map's format, the rule it
    /// reads `XBAYROFF` / `YBAYROFF` with:
    ///
    /// - FITS: `FitsHeader::get_i32` (an `i64` cut to 32 bits, else a `parse_fits_f64` number truncated to `i32`),
    ///   widened to `i64`. `'1'` is 1 (the card rule unquotes it), `0.6` is 0, `1D0` is 1.
    /// - XISF: the scanner's `parse_offset` (an `i64`, else an `f64` truncated to `i64`; no `D` exponent). `0.6` is
    ///   0, `1D0` is `None`.
    ///
    /// The key is looked up as [`HeaderKeyMap::get`] looks it up (normalised); the catalog's XISF scanner line looks
    /// the name up as written (see [`super::read_frame_facts`]).
    pub fn integer(&self, key: &str) -> Option<i64> {
        let s = self.get(key)?.trim();
        match self.format {
            HeaderFormat::Fits => s
                .parse::<i64>()
                .ok()
                .map(|n| n as i32)
                .or_else(|| parse_fits_f64(s).map(|f| f as i32))
                .map(i64::from),
            HeaderFormat::Xisf => s.parse::<i64>().ok().or_else(|| s.parse::<f64>().ok().map(|f| f as i64)),
        }
    }

    /// The map the older readers returned: UPPERCASE keys, no `COMMENT` / `HISTORY`, values
    /// trimmed and blank ones dropped.
    pub fn into_legacy_map(self) -> HashMap<String, String> {
        self.map
            .into_iter()
            .filter(|(k, _)| k != "COMMENT" && k != "HISTORY")
            .filter_map(|(k, v)| {
                let t = v.trim();
                if t.is_empty() {
                    None
                } else {
                    Some((k, t.to_string()))
                }
            })
            .collect()
    }

    /// The same map, read as `format`.
    pub(crate) fn with_format(mut self, format: HeaderFormat) -> Self {
        self.format = format;
        self
    }

    /// Card order, first occurrence wins (what a scan does).
    pub(crate) fn from_fits_header(h: &FitsHeader) -> Self {
        Self::from_first_wins(h.keywords_in_order().map(|(k, v)| (k.to_string(), v.to_string())))
    }

    /// Later pairs replace earlier ones (XISF `FITSKeyword` elements and dump lines). The map reads as FITS until
    /// [`HeaderKeyMap::with_format`] says otherwise.
    pub(crate) fn from_last_wins(pairs: impl IntoIterator<Item = (String, String)>) -> Self {
        let mut map = HashMap::new();
        for (k, v) in pairs {
            if let Some(k) = normalize_keyword(&k) {
                map.insert(k, v);
            }
        }
        Self { map, format: HeaderFormat::Fits }
    }

    /// Earlier pairs win. The map reads as FITS until [`HeaderKeyMap::with_format`] says otherwise.
    pub(crate) fn from_first_wins(pairs: impl IntoIterator<Item = (String, String)>) -> Self {
        let mut map = HashMap::new();
        for (k, v) in pairs {
            if let Some(k) = normalize_keyword(&k) {
                map.entry(k).or_insert(v);
            }
        }
        Self { map, format: HeaderFormat::Fits }
    }
}

/// The key map of a stored header text of the given format: the text [`super::header_text`] returns, which is the
/// text Athenaeum stores in its catalog (`fits_header.header`) for the file.
pub fn header_key_map(format: HeaderFormat, text: &str) -> HeaderKeyMap {
    stored_header_key_map(format, text).with_format(format)
}

/// Athenaeum's `stored_header_key_map`, body unchanged.
fn stored_header_key_map(format: HeaderFormat, text: &str) -> HeaderKeyMap {
    match format {
        HeaderFormat::Fits => {
            let keys = HeaderKeyMap::from_fits_header(&FitsHeader::from_header_text(text));
            if keys.has_keywords() {
                keys
            } else {
                // Hand-written `KEY = value` text, never a scanner output.
                HeaderKeyMap::from_first_wins(legacy_eq_text_pairs(text))
            }
        }
        HeaderFormat::Xisf => {
            if let Some(rest) = text.strip_prefix(XISF_DUMP_BANNER) {
                return HeaderKeyMap::from_last_wins(xisf_dump_pairs(rest));
            }
            if text.contains("<?xml") {
                match xisf_keyword_pairs(text.as_bytes(), None) {
                    Ok(pairs) => {
                        let keys = HeaderKeyMap::from_last_wins(pairs);
                        if keys.has_keywords() {
                            return keys;
                        }
                    }
                    Err(e) => tracing::warn!(
                        error = %e,
                        "stored XISF header xml unreadable; reading as text"
                    ),
                }
            }
            HeaderKeyMap::from_first_wins(legacy_eq_text_pairs(text))
        }
    }
}

/// Inverse of `xisf_header_text`'s body: one `NAME = value` per line, split at the first
/// `" = "`; a line ending in `" ="` has an empty value. Values are verbatim (no comment split).
fn xisf_dump_pairs(body: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for raw in body.lines() {
        let line = raw.strip_suffix('\r').unwrap_or(raw);
        if line.trim().is_empty() {
            continue;
        }
        if let Some((name, value)) = line.split_once(" = ") {
            out.push((name.to_string(), value.to_string()));
        } else if let Some(name) = line.strip_suffix(" =") {
            out.push((name.to_string(), String::new()));
        }
    }
    out
}

/// A plain `KEY = value` header dump (one keyword per line; banner, separator, COMMENT, HISTORY
/// and blank lines skipped). Quoted strings lose their quotes; an unquoted value's trailing
/// `" / comment"` is dropped. Hand-written text only: no scanner output takes this path.
fn legacy_eq_text_pairs(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
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
        out.push((key, value));
    }
    out
}

/// Athenaeum's tests of `header_keys.rs`, bodies unchanged. Two names map Athenaeum's onto this module's: the format
/// constants `FileFormat::FITS` / `FileFormat::XISF`, and `stored_header_key_map` (here the private body behind
/// [`header_key_map`]).
#[cfg(test)]
mod tests {
    use super::*;

    /// Athenaeum's `FileFormat` constants, as this module's [`HeaderFormat`].
    #[allow(non_snake_case)]
    mod FileFormat {
        pub(super) const FITS: super::HeaderFormat = super::HeaderFormat::Fits;
        pub(super) const XISF: super::HeaderFormat = super::HeaderFormat::Xisf;
    }

    #[test]
    fn normalize_keyword_trims_uppercases_and_refuses_empty_or_equals() {
        assert_eq!(normalize_keyword(" det-temp "), Some("DET-TEMP".into()));
        assert_eq!(normalize_keyword("eso det temp"), Some("ESO DET TEMP".into()));
        assert_eq!(normalize_keyword("   "), None);
        assert_eq!(normalize_keyword("A=B"), None);
    }

    #[test]
    fn xisf_dump_reader_inverts_xisf_header_text() {
        let dump = format!("{XISF_DUMP_BANNER}FILTER = Ha / 7nm\nGAIN = 100\nGAIN = 120\nOBJECT = \n");
        let m = stored_header_key_map(FileFormat::XISF, &dump);
        assert_eq!(m.get("FILTER"), Some("Ha / 7nm"));
        assert_eq!(m.get("GAIN"), Some("120"));
        assert_eq!(m.get("OBJECT"), Some(""));
    }

    #[test]
    fn xisf_dump_reader_reads_a_trimmed_empty_value() {
        let dump = format!("{XISF_DUMP_BANNER}OBJECT =\nGAIN = 5");
        let m = stored_header_key_map(FileFormat::XISF, &dump);
        assert_eq!(m.get("OBJECT"), Some(""));
        assert_eq!(m.get("GAIN"), Some("5"));
    }

    #[test]
    fn xisf_placeholder_has_no_keywords() {
        let p = "Header extracted from: /a=b/x.xisf\n(No FITS keywords found in XISF)";
        assert!(!stored_header_key_map(FileFormat::XISF, p).has_keywords());
    }

    #[test]
    fn fits_hand_dump_falls_back_to_eq_text() {
        assert_eq!(
            stored_header_key_map(FileFormat::FITS, "SIMPLE = T\nFOCALLEN = 140.6").get("FOCALLEN"),
            Some("140.6")
        );
    }

    #[test]
    fn hierarch_keys_match_case_insensitively() {
        let card = format!("{:<80}", "HIERARCH Eso Det Temp = -7.25 / sensor");
        let m = stored_header_key_map(FileFormat::FITS, &card);
        assert_eq!(m.get("ESO DET TEMP"), Some("-7.25"));
        assert_eq!(m.get("eso det temp"), Some("-7.25"));
    }

    #[test]
    fn xisf_raw_xml_reads_the_scanners_extraction() {
        let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<xisf>
  <FITSKeyword name="Filter" value="'Ha / 7nm'" comment="" />
  <FITSKeyword name="GAIN" value="100" />
  <FITSKeyword name="GAIN" value="120" />
</xisf>"#;
        let m = stored_header_key_map(FileFormat::XISF, xml);
        assert_eq!(m.get("FILTER"), Some("Ha / 7nm"));
        assert_eq!(m.get("GAIN"), Some("120"));
    }

    #[test]
    fn legacy_map_drops_comment_history_and_trims() {
        let cards = [
            format!("{:<80}", "COMMENT some text"),
            format!("{:<80}", "HISTORY something"),
            format!("{:<80}", "OBJECT  = ''"),
            format!("{:<80}", "FILTER  = '  Ha  '"),
            format!("{:<80}", "END"),
        ]
        .join("\n");
        let legacy = stored_header_key_map(FileFormat::FITS, &cards).into_legacy_map();
        assert!(!legacy.contains_key("COMMENT"));
        assert!(!legacy.contains_key("HISTORY"));
        assert!(!legacy.contains_key("OBJECT"));
        assert_eq!(legacy.get("FILTER").map(String::as_str), Some("Ha"));
    }
}

/// The two typed readers this crate adds ([`HeaderKeyMap::number`], [`HeaderKeyMap::integer`]), each against the
/// catalog reader it follows.
#[cfg(test)]
mod typed_value_tests {
    use super::*;

    fn fits_map(cards: &[&str]) -> HeaderKeyMap {
        let text: Vec<String> = cards.iter().map(|c| format!("{c:<80}")).collect();
        header_key_map(HeaderFormat::Fits, &text.join("\n"))
    }

    fn xisf_map(pairs: &[(&str, &str)]) -> HeaderKeyMap {
        let mut dump = XISF_DUMP_BANNER.to_string();
        for (k, v) in pairs {
            dump.push_str(&format!("{k} = {v}\n"));
        }
        header_key_map(HeaderFormat::Xisf, &dump)
    }

    #[test]
    fn a_fits_number_and_integer_read_as_fits_header_reads_them() {
        let cards = [
            "EGAIN   = 1.5D0 / e-/ADU",
            "EGAIN   = 2.5",
            "XBAYROFF= '1       '",
            "YBAYROFF=                  0.6",
            "BITPIX  =                 16.0",
            "BIG     =           4294967297",
        ];
        let m = fits_map(&cards);
        let block: Vec<u8> = cards
            .iter()
            .chain(&["END"])
            .flat_map(|c| format!("{c:<80}").into_bytes())
            .chain(std::iter::repeat(b' '))
            .take(2880)
            .collect();
        let h = FitsHeader::from_bytes(&block).unwrap();
        assert_eq!(m.format(), HeaderFormat::Fits);
        assert_eq!(m.number("EGAIN"), h.get_f64("EGAIN"));
        assert_eq!(m.number("EGAIN"), Some(1.5), "the first card wins, D exponent read");
        for key in ["XBAYROFF", "YBAYROFF", "BITPIX", "BIG"] {
            assert_eq!(m.integer(key), h.get_i32(key).map(i64::from), "{key}");
        }
        assert_eq!(m.integer("XBAYROFF"), Some(1));
        assert_eq!(m.integer("YBAYROFF"), Some(0), "a fraction is truncated");
        assert_eq!(m.integer("BIG"), Some(1), "cut to 32 bits, as get_i32 cuts it");
        assert_eq!(m.integer("ABSENT"), None);
    }

    #[test]
    fn an_xisf_number_and_integer_read_by_the_xisf_scanner_rules() {
        let m = xisf_map(&[
            ("EGAIN", "1.5"),
            ("EGAIN", "2.5"),
            ("DATAMAX", "6.5D4"),
            ("XBAYROFF", "0.6"),
            ("YBAYROFF", "1D0"),
            ("NAXIS1", "4294967297"),
        ]);
        assert_eq!(m.format(), HeaderFormat::Xisf);
        assert_eq!(m.number("EGAIN"), Some(2.5), "the last keyword wins");
        assert_eq!(m.number("DATAMAX"), Some(65000.0), "parse_fits_f64 reads the D exponent");
        assert_eq!(m.integer("XBAYROFF"), Some(0), "a fraction is truncated");
        assert_eq!(m.integer("YBAYROFF"), None, "no D exponent in parse_offset");
        assert_eq!(m.integer("NAXIS1"), Some(4_294_967_297), "no cut to 32 bits");
        assert_eq!(m.integer("xbayroff"), Some(0), "looked up normalised");
    }
}
