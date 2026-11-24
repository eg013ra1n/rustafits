//! A frame's facts, read by the lines of Athenaeum's scanner that build its catalog row (`athenaeum-core`'s
//! `fits_parser::build_frame_from_header` for FITS, `parse_xisf_xml` for XISF, at ff7e0ca3), under the default rules.
//!
//! Each field is read the way the scanner reads the column it fills, so a fact equals what the catalog stores for the
//! file under Athenaeum's default header rules:
//!
//! - the twelve rule-driven fields through [`resolve_default`] over the file's key map (FITS: the cards, the first of a
//!   keyword winning; XISF: the `<FITSKeyword>` elements, the last winning), and for XISF the pixel-size
//!   `<Property>` after the whole chain;
//! - FITS: every other field through `FitsHeader`'s getters (`get_str`: trimmed, blank is absent; `get_i32`), keyed
//!   as the card rule keys them (HIERARCH names as written);
//! - XISF: every other field from the scanner's own keyword map, keyed by the element's `name` as written, the last
//!   element of a name winning, the value verbatim (`IMAGETYP`, `DATE-OBS` and `BAYERPAT` untrimmed, as the scanner
//!   hands them on), with the scanner's fallbacks: `<Image geometry>` for the size, the binning `<Property>`s,
//!   `<ColorFilterArray pattern>` for a missing or blank `BAYERPAT`;
//! - `imagetyp`: Athenaeum's `ImageType::from_str` (upper-cased, the names its `ImageType` serializes to), `IMAGETYP`
//!   first, then `FRAME`;
//! - `origin`: Athenaeum's `calibration::geometry::origin_from_keys` (`XORGSUBF` and `YORGSUBF`, both or neither, an
//!   integral real counting); `calstat`: the trimmed, non-blank `CALSTAT` its calibrated-light check reads; `egain`:
//!   [`HeaderKeyMap::number`].
//!
//! Not the scanner's (it stores no such column): `planes` and `bitpix`. FITS reads them from `NAXIS` / `NAXIS3` and
//! `BITPIX`; XISF from the third `<Image geometry>` component and the `sampleFormat` of the same element (`UInt8` is 8,
//! `UInt16` 16, `UInt32` 32, `Float32` -32, `Float64` -64). A compressed or extension-HDU FITS file is read by its
//! logical header, so its shape is the image's.

use std::collections::HashMap;
use std::path::Path;

use anyhow::Result;

use super::chains::{resolve_default, HeaderField};
use super::fits_header_reader::FitsHeader;
use super::keys::HeaderKeyMap;
use super::text::{detect_format, fits_header, read_xisf_header, xisf_keyword_element, xisf_xml_slice};
use super::HeaderFormat;

/// What a frame's header says, read by the catalog's rules (see the module docs of [`super`]).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize), serde(rename_all = "camelCase"))]
pub struct FrameFacts {
    /// The header's format.
    pub format: HeaderFormat,
    /// Image width, pixels (`NAXIS1`; XISF: the keyword, else `<Image geometry>`).
    pub width: Option<i64>,
    /// Image height, pixels (`NAXIS2`; XISF: the keyword, else `<Image geometry>`).
    pub height: Option<i64>,
    /// Image planes: 1 for a 2-axis FITS image, `NAXIS3` for more axes; XISF: the geometry's channel count.
    pub planes: Option<i64>,
    /// FITS `BITPIX`; XISF: the `sampleFormat` as a `BITPIX`.
    pub bitpix: Option<i64>,
    /// The frame type, normalized as Athenaeum's `ImageType` names it (`Light`, `Dark`, `Flat`, `Bias`, `DarkFlat`,
    /// `MasterLight`, `MasterDark`, `MasterFlat`, `MasterBias`, `MasterDarkFlat`).
    pub imagetyp: Option<String>,
    /// Exposure, seconds (the `Exptime` chain).
    pub exptime_s: Option<f64>,
    /// Camera gain setting (the `Gain` chain).
    pub gain: Option<f64>,
    /// Camera offset setting (the `Offset` chain).
    pub offset: Option<f64>,
    /// `EGAIN`, e- per ADU, as [`HeaderKeyMap::number`] reads it.
    pub egain: Option<f64>,
    /// Sensor temperature, °C (the `CcdTemp` chain).
    pub ccd_temp: Option<f64>,
    /// Cooler set point, °C (the `SetTemp` chain).
    pub set_temp: Option<f64>,
    /// Filter name (the `Filter` chain).
    pub filter: Option<String>,
    /// Target name (the `Object` chain).
    pub object: Option<String>,
    /// Camera name (the `Instrume` chain).
    pub instrume: Option<String>,
    /// Telescope name (the `Telescop` chain).
    pub telescop: Option<String>,
    /// Focal length, mm (the `Focallen` chain).
    pub focallen_mm: Option<f64>,
    /// Pixel width, µm (the `Xpixsz` chain; XISF: then the sensor `<Property>`).
    pub xpixsz_um: Option<f64>,
    /// Pixel height, µm (the `Ypixsz` chain; XISF: then the sensor `<Property>`).
    pub ypixsz_um: Option<f64>,
    /// `XBINNING` (XISF: then the camera `<Property>`).
    pub xbinning: Option<i64>,
    /// `YBINNING` (XISF: then the camera `<Property>`).
    pub ybinning: Option<i64>,
    /// `DATE-OBS`, the text the scanner parses, unparsed.
    pub date_obs: Option<String>,
    /// `TIME-OBS`, the text the scanner parses with `DATE-OBS`, unparsed.
    pub time_obs: Option<String>,
    /// `BAYERPAT` (XISF: else `<ColorFilterArray pattern>`), as the catalog stores it: not folded with the offsets.
    pub bayerpat: Option<String>,
    /// `XBAYROFF`, the CFA phase offset; absent stays `None`.
    pub xbayroff: Option<i64>,
    /// `YBAYROFF`, the CFA phase offset; absent stays `None`.
    pub ybayroff: Option<i64>,
    /// `ROWORDER`, trimmed.
    pub roworder: Option<String>,
    /// `(XORGSUBF, YORGSUBF)`, the sub-frame origin, both or neither.
    pub origin: Option<(i64, i64)>,
    /// `CALSTAT`, the calibration state a calibrated frame declares.
    pub calstat: Option<String>,
}

/// The facts of the FITS or XISF file at `path`, read by the catalog's scanner lines under the default rules. Reads
/// the header only (a gzipped FITS file is inflated to find its image header).
pub fn read_frame_facts(path: &Path) -> Result<FrameFacts> {
    let facts = match detect_format(path) {
        HeaderFormat::Fits => facts_from_fits_header(&fits_header(path)?),
        HeaderFormat::Xisf => facts_from_xisf_xml(&read_xisf_header(path)?, path)?,
    };
    tracing::debug!(path = %path.display(), format = ?facts.format, "frame facts read");
    Ok(facts)
}

/// Athenaeum's `ImageType::from_str`, giving the variant's name.
fn image_type(s: &str) -> Option<&'static str> {
    let s_upper = s.to_uppercase();
    match s_upper.as_str() {
        // Light frames: IRAF, SBFITSEXT, INDI FRAME keyword
        "LIGHT" | "LIGHT FRAME" => Some("Light"),
        // Dark frames: IRAF, SBFITSEXT
        "DARK" | "DARK FRAME" => Some("Dark"),
        // Flat frames: IRAF, SBFITSEXT ("Flat Field"), INDI ("Flat Frame"), observatory variants
        "FLAT" | "FLAT FIELD" | "FLAT FRAME" | "FLATFIELD"
            | "SKY FLAT" | "DOME FLAT" | "TWILIGHT FLAT" => Some("Flat"),
        // Bias frames: IRAF, SBFITSEXT, "Offset" alias
        "BIAS" | "BIAS FRAME" | "OFFSET" | "OFFSET FRAME" => Some("Bias"),
        // Dark flats: various software
        "DARKFLAT" | "DARK FLAT" | "DARK FLAT FRAME"
            | "FLATDARK" | "FLAT DARK" | "FLAT DARK FRAME" => Some("DarkFlat"),
        // Master frames
        "MASTER LIGHT" | "MASTERLIGHT" | "MASTER LIGHT FRAME" => Some("MasterLight"),
        "MASTER DARK" | "MASTERDARK" | "MASTER DARK FRAME" => Some("MasterDark"),
        "MASTER FLAT" | "MASTERFLAT" | "MASTER FLAT FIELD" | "MASTER FLAT FRAME" => Some("MasterFlat"),
        "MASTER BIAS" | "MASTERBIAS" | "MASTER BIAS FRAME"
            | "MASTER OFFSET" | "MASTEROFFSET" => Some("MasterBias"),
        "MASTER DARK FLAT" | "MASTERDARKFLAT" | "MASTER DARKFLAT"
            | "MASTER FLAT DARK" | "MASTERFLATDARK" | "MASTER FLATDARK" => Some("MasterDarkFlat"),
        _ => None,
    }
}

/// The `IMAGETYP` → `FRAME` fallback of Athenaeum's `finalize_frame`.
fn frame_type(imagetyp: Option<&str>, frame: Option<&str>) -> Option<String> {
    imagetyp
        .and_then(image_type)
        .or_else(|| frame.and_then(image_type))
        .map(str::to_string)
}

/// An integer card value; an integral real (`1064.`) counts, anything else does not (Athenaeum's
/// `calibration::geometry::card_int`).
fn card_int(value: &str) -> Option<i32> {
    let v = value.trim();
    v.parse::<i32>().ok().or_else(|| {
        v.parse::<f64>()
            .ok()
            .filter(|f| f.is_finite() && f.fract() == 0.0 && f.abs() <= i32::MAX as f64)
            .map(|f| f as i32)
    })
}

/// `(XORGSUBF, YORGSUBF)`, both or nothing (Athenaeum's `origin_from_keys`, over the trimmed, non-blank values its
/// legacy key map holds).
fn origin(keys: &HeaderKeyMap) -> Option<(i64, i64)> {
    let read = |key: &str| keys.get(key).map(str::trim).filter(|v| !v.is_empty()).and_then(card_int);
    read("XORGSUBF").zip(read("YORGSUBF")).map(|(x, y)| (i64::from(x), i64::from(y)))
}

/// A trimmed, non-blank value (the legacy key map's `non_empty`).
fn non_blank(keys: &HeaderKeyMap, key: &str) -> Option<String> {
    keys.get(key).map(str::trim).filter(|v| !v.is_empty()).map(str::to_string)
}

fn number(keys: &HeaderKeyMap, field: HeaderField) -> Option<f64> {
    resolve_default(field, keys).and_then(|v| v.as_number())
}

fn text(keys: &HeaderKeyMap, field: HeaderField) -> Option<String> {
    resolve_default(field, keys).and_then(|v| v.as_text().map(str::to_string))
}

/// The facts of a FITS header (`build_frame_from_header`'s lines).
pub(crate) fn facts_from_fits_header(header: &FitsHeader) -> FrameFacts {
    let keys = HeaderKeyMap::from_fits_header(header).with_format(HeaderFormat::Fits);
    let int = |key: &str| header.get_i32(key).map(i64::from);
    let planes = match header.get_i32("NAXIS") {
        Some(n) if n >= 3 => int("NAXIS3"),
        Some(2) => Some(1),
        _ => None,
    };
    FrameFacts {
        format: HeaderFormat::Fits,
        width: int("NAXIS1"),
        height: int("NAXIS2"),
        planes,
        bitpix: int("BITPIX"),
        imagetyp: frame_type(header.get_str("IMAGETYP").as_deref(), header.get_str("FRAME").as_deref()),
        exptime_s: number(&keys, HeaderField::Exptime),
        gain: number(&keys, HeaderField::Gain),
        offset: number(&keys, HeaderField::Offset),
        egain: keys.number("EGAIN"),
        ccd_temp: number(&keys, HeaderField::CcdTemp),
        set_temp: number(&keys, HeaderField::SetTemp),
        filter: text(&keys, HeaderField::Filter),
        object: text(&keys, HeaderField::Object),
        instrume: text(&keys, HeaderField::Instrume),
        telescop: text(&keys, HeaderField::Telescop),
        focallen_mm: number(&keys, HeaderField::Focallen),
        xpixsz_um: number(&keys, HeaderField::Xpixsz),
        ypixsz_um: number(&keys, HeaderField::Ypixsz),
        xbinning: int("XBINNING"),
        ybinning: int("YBINNING"),
        date_obs: header.get_str("DATE-OBS"),
        time_obs: header.get_str("TIME-OBS"),
        bayerpat: header.get_str("BAYERPAT"),
        xbayroff: int("XBAYROFF"),
        ybayroff: int("YBAYROFF"),
        roworder: header.get_str("ROWORDER"),
        origin: origin(&keys),
        calstat: non_blank(&keys, "CALSTAT"),
    }
}

/// The `BITPIX` an XISF `sampleFormat` stands for.
fn sample_format_bitpix(format: &str) -> Option<i64> {
    match format {
        "UInt8" => Some(8),
        "UInt16" => Some(16),
        "UInt32" => Some(32),
        "Float32" => Some(-32),
        "Float64" => Some(-64),
        _ => None,
    }
}

/// The facts of an XISF XML header block (`parse_xisf_xml`'s loop and lines). `path` names the file in logs.
pub(crate) fn facts_from_xisf_xml(content: &[u8], path: &Path) -> Result<FrameFacts> {
    use quick_xml::events::Event;
    use quick_xml::Reader;

    // The XML section: from "<?xml" through </xisf>
    let xml_content = xisf_xml_slice(content)?;
    let xml_str = String::from_utf8_lossy(xml_content);

    let mut reader = Reader::from_str(&xml_str);
    reader.config_mut().trim_text(true);
    // Permissive dangling-`&` retained deliberately (0.36-era behavior; XISF
    // headers are machine-written but a malformed one should still render —
    // cycle decision 2026-07-29).
    reader.config_mut().allow_dangling_amp = true;

    let mut fits_keywords: HashMap<String, String> = HashMap::new();
    // Every FITSKeyword in document order: the one input of the resolver (last duplicate wins,
    // as `fits_keywords` below).
    let mut keyword_pairs: Vec<(String, String)> = Vec::new();
    let mut xisf_geometry: Option<(i32, i32)> = None;
    // Not the scanner's: the channel count and sample format of the element the geometry came from.
    let mut xisf_shape: (Option<i64>, Option<i64>) = (None, None);
    let mut xisf_properties: HashMap<String, String> = HashMap::new();
    let mut xisf_cfa_pattern: Option<String> = None;
    let mut buf = Vec::new();

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Empty(ref e)) if e.name().as_ref() == b"FITSKeyword" => {
                if let Some((name, value)) = xisf_keyword_element(e) {
                    keyword_pairs.push((name.clone(), value.clone()));
                    fits_keywords.insert(name, value);
                }
            }
            // Parse <Image geometry="w:h:c"> for NAXIS1/NAXIS2 fallback
            Ok(Event::Start(ref e)) | Ok(Event::Empty(ref e))
                if e.name().as_ref() == b"Image" =>
            {
                let mut sample_format = None;
                let mut geometry = None;
                for attr in e.attributes().flatten() {
                    if attr.key.as_ref() == b"geometry" {
                        let geom = String::from_utf8_lossy(&attr.value);
                        let parts: Vec<&str> = geom.split(':').collect();
                        if parts.len() >= 2 {
                            if let (Ok(w), Ok(h)) = (parts[0].parse::<i32>(), parts[1].parse::<i32>()) {
                                let channels = parts.get(2).and_then(|c| c.parse::<i64>().ok());
                                geometry = Some(((w, h), channels));
                            }
                        }
                    } else if attr.key.as_ref() == b"sampleFormat" {
                        sample_format = sample_format_bitpix(&String::from_utf8_lossy(&attr.value));
                    }
                }
                if let Some((wh, channels)) = geometry {
                    xisf_geometry = Some(wh);
                    xisf_shape = (channels, sample_format);
                }
            }
            // Parse <Property> elements for native XISF metadata
            Ok(Event::Empty(ref e)) if e.name().as_ref() == b"Property" => {
                let mut prop_id = String::new();
                let mut prop_value = String::new();
                for attr in e.attributes().flatten() {
                    match attr.key.as_ref() {
                        b"id" => prop_id = String::from_utf8_lossy(&attr.value).to_string(),
                        b"value" => prop_value = String::from_utf8_lossy(&attr.value).to_string(),
                        _ => {}
                    }
                }
                if !prop_id.is_empty() && !prop_value.is_empty() {
                    xisf_properties.insert(prop_id, prop_value);
                }
            }
            // Native XISF CFA declaration: <ColorFilterArray pattern="RGGB"
            // width="2" height="2" name="RGGB"/>. Taken raw here; validated and
            // adopted below, only when no BAYERPAT keyword is present.
            Ok(Event::Start(ref e)) | Ok(Event::Empty(ref e))
                if e.name().as_ref() == b"ColorFilterArray" =>
            {
                for attr in e.attributes().flatten() {
                    if attr.key.as_ref() == b"pattern" {
                        xisf_cfa_pattern = Some(String::from_utf8_lossy(&attr.value).to_string());
                    }
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => {
                tracing::error!(
                    path = %path.display(),
                    position = reader.buffer_position(),
                    error = %e,
                    "error parsing XISF XML; using partial result"
                );
                break;
            }
            _ => {}
        }
        buf.clear();
    }

    tracing::debug!(path = %path.display(), count = fits_keywords.len(), "extracted FITS keywords from XISF");

    // The twelve rule-driven fields: one resolver, the same one the FITS path calls. The
    // XISF pixel-size `<Property>` stays a last resort after the whole keyword order.
    let keys = HeaderKeyMap::from_last_wins(keyword_pairs).with_format(HeaderFormat::Xisf);

    let xbinning = fits_keywords.get("XBINNING")
        .and_then(|s| s.parse::<i32>().ok())
        .or_else(|| xisf_properties.get("Instrument:Camera:XBinning")
            .and_then(|s| s.parse::<i32>().ok()));
    let ybinning = fits_keywords.get("YBINNING")
        .and_then(|s| s.parse::<i32>().ok())
        .or_else(|| xisf_properties.get("Instrument:Camera:YBinning")
            .and_then(|s| s.parse::<i32>().ok()));
    // A blank keyword is not a declaration (the scanner's rule).
    let bayerpat = fits_keywords.get("BAYERPAT").cloned().filter(|s| !s.trim().is_empty());
    // Fall back to the native XISF <ColorFilterArray> element. An explicit BAYERPAT keyword keeps precedence. The
    // element carries no phase offsets, so xbayroff/ybayroff stay None below.
    let bayerpat = bayerpat.or_else(|| {
        let raw = xisf_cfa_pattern.as_ref()?;
        let pattern = raw.trim().to_uppercase();
        if !pattern.is_empty() && pattern.bytes().all(|b| matches!(b, b'R' | b'G' | b'B')) {
            tracing::debug!(path = %path.display(), pattern = %pattern, "adopted cfa pattern from xisf element");
            Some(pattern)
        } else {
            tracing::warn!(path = %path.display(), pattern = %raw, "unrecognized xisf cfa pattern");
            None
        }
    });
    // CFA phase offsets + row order. Absent keywords stay None; 0 is a real phase, not a default.
    let parse_offset = |s: &String| -> Option<i64> {
        let s = s.trim();
        s.parse::<i64>().ok().or_else(|| s.parse::<f64>().ok().map(|f| f as i64))
    };
    let xbayroff = fits_keywords.get("XBAYROFF").and_then(parse_offset);
    let ybayroff = fits_keywords.get("YBAYROFF").and_then(parse_offset);
    let roworder = fits_keywords.get("ROWORDER")
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty());
    let xpixsz = number(&keys, HeaderField::Xpixsz)
        .or_else(|| xisf_properties.get("Instrument:Sensor:XPixelSize")
            .and_then(|s| s.parse::<f64>().ok()));
    let ypixsz = number(&keys, HeaderField::Ypixsz)
        .or_else(|| xisf_properties.get("Instrument:Sensor:YPixelSize")
            .and_then(|s| s.parse::<f64>().ok()));

    // Image dimensions (with <Image geometry="w:h:c"> fallback)
    let naxis1 = fits_keywords.get("NAXIS1")
        .and_then(|s| s.parse::<i32>().ok())
        .or_else(|| xisf_geometry.map(|(w, _)| w));
    let naxis2 = fits_keywords.get("NAXIS2")
        .and_then(|s| s.parse::<i32>().ok())
        .or_else(|| xisf_geometry.map(|(_, h)| h));

    Ok(FrameFacts {
        format: HeaderFormat::Xisf,
        width: naxis1.map(i64::from),
        height: naxis2.map(i64::from),
        planes: xisf_shape.0,
        bitpix: xisf_shape.1,
        imagetyp: frame_type(
            fits_keywords.get("IMAGETYP").map(String::as_str),
            fits_keywords.get("FRAME").map(String::as_str),
        ),
        exptime_s: number(&keys, HeaderField::Exptime),
        gain: number(&keys, HeaderField::Gain),
        offset: number(&keys, HeaderField::Offset),
        egain: keys.number("EGAIN"),
        ccd_temp: number(&keys, HeaderField::CcdTemp),
        set_temp: number(&keys, HeaderField::SetTemp),
        filter: text(&keys, HeaderField::Filter),
        object: text(&keys, HeaderField::Object),
        instrume: text(&keys, HeaderField::Instrume),
        telescop: text(&keys, HeaderField::Telescop),
        focallen_mm: number(&keys, HeaderField::Focallen),
        xpixsz_um: xpixsz,
        ypixsz_um: ypixsz,
        xbinning: xbinning.map(i64::from),
        ybinning: ybinning.map(i64::from),
        date_obs: fits_keywords.get("DATE-OBS").cloned(),
        time_obs: fits_keywords.get("TIME-OBS").cloned(),
        bayerpat,
        xbayroff,
        ybayroff,
        roworder,
        origin: origin(&keys),
        calstat: non_blank(&keys, "CALSTAT"),
    })
}
