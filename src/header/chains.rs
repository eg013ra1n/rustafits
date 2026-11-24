//! The catalog's twelve rule-driven fields and their default keyword chains: Athenaeum's
//! `athenaeum-core/src/header_fallbacks/mod.rs` (`HeaderField`, `FieldKind`, `default_chain`) and `resolve.rs`
//! (`resolve_keyword`, `resolve_field`) at ff7e0ca3, with the default-chain tests of `resolve.rs`.
//!
//! Athenaeum's user can reorder a field's keywords (`HeaderFallbackConfig`); those rules stay Athenaeum's (P2-R1), so
//! this module knows the default chains only: [`resolve_default`] is Athenaeum's `resolve_field` over
//! [`HeaderField::default_chain`]. Athenaeum's `FieldValue` is [`HeaderValue`] here.

use super::fits_header_reader::parse_fits_f64;
use super::keys::{normalize_keyword, HeaderKeyMap};

/// Number or text: how a keyword's raw value is read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FieldKind {
    Number,
    Text,
}

/// The twelve catalog fields whose header keywords Athenaeum's user orders. Declaration order is the order of
/// Athenaeum's Settings section.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum HeaderField {
    /// Sensor temperature, °C (`CCD-TEMP`, then `DET-TEMP`).
    CcdTemp,
    /// Cooler set point, °C (`SET-TEMP`).
    SetTemp,
    /// Camera gain setting (`GAIN`).
    Gain,
    /// Camera offset setting (`OFFSET`).
    Offset,
    /// Exposure, seconds (`EXPTIME`, then `EXPOSURE`).
    Exptime,
    /// Focal length, mm (`FOCALLEN`).
    Focallen,
    /// Pixel width, µm (`XPIXSZ`, then `PIXSIZE1`).
    Xpixsz,
    /// Pixel height, µm (`YPIXSZ`, then `PIXSIZE2`).
    Ypixsz,
    /// Target name (`OBJECT`).
    Object,
    /// Filter name (`FILTER`).
    Filter,
    /// Camera name (`INSTRUME`).
    Instrume,
    /// Telescope name (`TELESCOP`).
    Telescop,
}

impl HeaderField {
    /// The twelve fields, in declaration order.
    pub const ALL: [HeaderField; 12] = [
        HeaderField::CcdTemp,
        HeaderField::SetTemp,
        HeaderField::Gain,
        HeaderField::Offset,
        HeaderField::Exptime,
        HeaderField::Focallen,
        HeaderField::Xpixsz,
        HeaderField::Ypixsz,
        HeaderField::Object,
        HeaderField::Filter,
        HeaderField::Instrume,
        HeaderField::Telescop,
    ];

    /// The id Athenaeum uses on the wire and in its settings document (camelCase).
    pub fn id(self) -> &'static str {
        match self {
            HeaderField::CcdTemp => "ccdTemp",
            HeaderField::SetTemp => "setTemp",
            HeaderField::Gain => "gain",
            HeaderField::Offset => "offset",
            HeaderField::Exptime => "exptime",
            HeaderField::Focallen => "focallen",
            HeaderField::Xpixsz => "xpixsz",
            HeaderField::Ypixsz => "ypixsz",
            HeaderField::Object => "object",
            HeaderField::Filter => "filter",
            HeaderField::Instrume => "instrume",
            HeaderField::Telescop => "telescop",
        }
    }

    fn kind(self) -> FieldKind {
        match self {
            HeaderField::Object
            | HeaderField::Filter
            | HeaderField::Instrume
            | HeaderField::Telescop => FieldKind::Text,
            _ => FieldKind::Number,
        }
    }

    /// The keywords tried when the user has not changed the field, in order.
    pub fn default_chain(self) -> &'static [&'static str] {
        match self {
            HeaderField::CcdTemp => &["CCD-TEMP", "DET-TEMP"],
            HeaderField::SetTemp => &["SET-TEMP"],
            HeaderField::Gain => &["GAIN"],
            HeaderField::Offset => &["OFFSET"],
            HeaderField::Exptime => &["EXPTIME", "EXPOSURE"],
            HeaderField::Focallen => &["FOCALLEN"],
            HeaderField::Xpixsz => &["XPIXSZ", "PIXSIZE1"],
            HeaderField::Ypixsz => &["YPIXSZ", "PIXSIZE2"],
            HeaderField::Object => &["OBJECT"],
            HeaderField::Filter => &["FILTER"],
            HeaderField::Instrume => &["INSTRUME"],
            HeaderField::Telescop => &["TELESCOP"],
        }
    }
}

/// A resolved header value.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum HeaderValue {
    /// A text field's value, trimmed and not blank.
    Text(String),
    /// A number field's value, finite.
    Number(f64),
}

impl HeaderValue {
    /// The number, when this is one.
    pub fn as_number(&self) -> Option<f64> {
        match self {
            HeaderValue::Number(n) => Some(*n),
            HeaderValue::Text(_) => None,
        }
    }

    /// The text, when this is one.
    pub fn as_text(&self) -> Option<&str> {
        match self {
            HeaderValue::Text(t) => Some(t),
            HeaderValue::Number(_) => None,
        }
    }
}

/// The value one keyword gives for a field, `None` when the keyword is absent or its value has
/// the wrong type. A number is read with `parse_fits_f64` (Fortran `D` exponent included) and
/// must be finite (`NaN` / `inf` count as absent); text is trimmed and a blank one is absent.
fn resolve_keyword(field: HeaderField, keyword: &str, keys: &HeaderKeyMap) -> Option<HeaderValue> {
    let key = normalize_keyword(keyword)?;
    let raw = keys.get(&key)?.trim();
    match field.kind() {
        FieldKind::Number => parse_fits_f64(raw)
            .filter(|n| n.is_finite())
            .map(HeaderValue::Number),
        FieldKind::Text => {
            if raw.is_empty() {
                None
            } else {
                Some(HeaderValue::Text(raw.to_string()))
            }
        }
    }
}

/// The first usable value of `order`, first keyword first.
fn resolve_field<S: AsRef<str>>(field: HeaderField, order: &[S], keys: &HeaderKeyMap) -> Option<HeaderValue> {
    order.iter().find_map(|k| resolve_keyword(field, k.as_ref(), keys))
}

/// A field's value under the default rules: the first keyword of [`HeaderField::default_chain`] that gives a usable
/// value. An unparsable or non-finite number, or a blank text, falls through to the next keyword.
pub fn resolve_default(field: HeaderField, map: &HeaderKeyMap) -> Option<HeaderValue> {
    resolve_field(field, field.default_chain(), map)
}

/// Athenaeum's `resolve.rs` tests that use the default chains only, bodies unchanged. Not copied: the tests of a user
/// order (`user_first_order_wins_over_the_standard_keyword`, `default_keyword_demoted_to_the_end_is_tried_last`,
/// `a_case_or_space_variant_of_a_user_keyword_matches`) and of `same_value`, which this module does not have. Three
/// names map Athenaeum's onto this module's: `FieldValue` is [`HeaderValue`], and `HeaderFallbackConfig::default()` /
/// `resolve_all` are the default rules read through [`resolve_default`].
#[cfg(test)]
mod tests {
    use super::*;
    use HeaderValue as FieldValue;

    /// Athenaeum's default rules document: every field on its default chain. Not a unit struct, so the copied
    /// `HeaderFallbackConfig::default()` reads as written.
    #[derive(Default)]
    struct HeaderFallbackConfig {
        _default_rules: (),
    }

    /// Athenaeum's `ResolvedFields`: the twelve fields of a header (the tests read some of them).
    #[derive(Debug, Default)]
    #[allow(dead_code)]
    struct ResolvedFields {
        ccd_temp: Option<f64>,
        set_temp: Option<f64>,
        gain: Option<f64>,
        offset: Option<f64>,
        exptime: Option<f64>,
        focallen: Option<f64>,
        xpixsz: Option<f64>,
        ypixsz: Option<f64>,
        object: Option<String>,
        filter: Option<String>,
        instrume: Option<String>,
        telescop: Option<String>,
    }

    /// Athenaeum's `resolve_all` under the default rules.
    fn resolve_all(_rules: &HeaderFallbackConfig, keys: &HeaderKeyMap) -> ResolvedFields {
        let n = |f| resolve_default(f, keys).and_then(|v| v.as_number());
        let t = |f| resolve_default(f, keys).and_then(|v| v.as_text().map(str::to_string));
        ResolvedFields {
            ccd_temp: n(HeaderField::CcdTemp),
            set_temp: n(HeaderField::SetTemp),
            gain: n(HeaderField::Gain),
            offset: n(HeaderField::Offset),
            exptime: n(HeaderField::Exptime),
            focallen: n(HeaderField::Focallen),
            xpixsz: n(HeaderField::Xpixsz),
            ypixsz: n(HeaderField::Ypixsz),
            object: t(HeaderField::Object),
            filter: t(HeaderField::Filter),
            instrume: t(HeaderField::Instrume),
            telescop: t(HeaderField::Telescop),
        }
    }

    fn keys(pairs: &[(&str, &str)]) -> HeaderKeyMap {
        HeaderKeyMap::from_first_wins(pairs.iter().map(|(k, v)| (k.to_string(), v.to_string())))
    }

    fn order(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    fn ccd(keys: &HeaderKeyMap) -> Option<FieldValue> {
        resolve_field(HeaderField::CcdTemp, &order(&["CCD-TEMP", "DET-TEMP"]), keys)
    }

    #[test]
    fn resolver_takes_the_first_usable_keyword_in_order() {
        let k = keys(&[("CCD-TEMP", "-10"), ("DET-TEMP", "-5.5")]);
        assert_eq!(ccd(&k), Some(FieldValue::Number(-10.0)));
    }

    #[test]
    fn resolver_skips_a_value_of_the_wrong_type() {
        let k = keys(&[("CCD-TEMP", "n/a"), ("DET-TEMP", "-5.5")]);
        assert_eq!(ccd(&k), Some(FieldValue::Number(-5.5)));
    }

    #[test]
    fn resolver_treats_nan_inf_text_and_blank_as_absent() {
        for bad in ["NaN", "inf", "-INF", "n/a", "", "   "] {
            let k = keys(&[("CCD-TEMP", bad), ("DET-TEMP", "12")]);
            assert_eq!(ccd(&k), Some(FieldValue::Number(12.0)), "value {bad:?}");
        }
        let only_bad = keys(&[("CCD-TEMP", "NaN")]);
        assert_eq!(ccd(&only_bad), None);
    }

    #[test]
    fn text_values_are_trimmed_and_blank_is_absent() {
        let k = keys(&[("FILTER", "  Ha  "), ("OBJECT", "  ")]);
        assert_eq!(
            resolve_keyword(HeaderField::Filter, "FILTER", &k),
            Some(FieldValue::Text("Ha".into()))
        );
        assert_eq!(resolve_keyword(HeaderField::Object, "OBJECT", &k), None);
    }

    #[test]
    fn fortran_d_notation_parses_on_every_path() {
        let k = keys(&[("EXPTIME", "1.2D+02")]);
        assert_eq!(
            resolve_keyword(HeaderField::Exptime, "EXPTIME", &k),
            Some(FieldValue::Number(120.0))
        );
    }

    #[test]
    fn resolve_all_reads_the_default_chains() {
        let k = keys(&[("DET-TEMP", "-5.5"), ("EXPOSURE", "120"), ("PIXSIZE1", "3.76"), ("OBJECT", "M42")]);
        let r = resolve_all(&HeaderFallbackConfig::default(), &k);
        assert_eq!(r.ccd_temp, Some(-5.5));
        assert_eq!(r.exptime, Some(120.0));
        assert_eq!(r.xpixsz, Some(3.76));
        assert_eq!(r.object.as_deref(), Some("M42"));
        assert_eq!(r.gain, None);
    }
}
