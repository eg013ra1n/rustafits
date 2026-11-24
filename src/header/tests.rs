//! The module end to end on files: the logical header of a compressed file, the facts of a FITS and an XISF frame,
//! the stored text against the bytes, and a gzipped file against its plain form. Athenaeum's own tests of the ported
//! rules sit next to the code they test (`keys`, `fits_header_reader`, `chains`).

use std::path::{Path, PathBuf};

use super::*;

fn fpack(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("fpack").join(name)
}

/// An 80-column `KEY = value` card.
fn card(key: &str, value: &str) -> String {
    format!("{:<80}", format!("{key:<8}= {value:>20}"))
}

/// A 16-bit FITS image of `w` x `h` samples `data` (big-endian, as stored) whose header is `cards`, then `END`.
fn write_fits(path: &Path, cards: &[String], data: &[i16]) {
    let mut bytes: Vec<u8> = cards.iter().map(|c| format!("{c:<80}")).collect::<String>().into_bytes();
    bytes.extend(format!("{:<80}", "END").bytes());
    bytes.resize(bytes.len().div_ceil(2880) * 2880, b' ');
    bytes.extend(data.iter().flat_map(|v| v.to_be_bytes()));
    bytes.resize(bytes.len().div_ceil(2880) * 2880, 0);
    std::fs::write(path, bytes).unwrap();
}

/// An XISF file whose `<Image>` holds `children`, with `w * h` UInt16 samples.
fn write_xisf(path: &Path, geometry: &str, children: &str, samples: usize) {
    let xml = |pos: usize| {
        format!(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?><xisf version=\"1.0\">\
<Image geometry=\"{geometry}\" sampleFormat=\"UInt16\" colorSpace=\"Gray\" byteOrder=\"little\" \
location=\"attachment:{pos:08}:{}\">{children}</Image></xisf>",
            samples * 2
        )
    };
    let len = xml(0).len();
    let mut bytes = b"XISF0100".to_vec();
    bytes.extend_from_slice(&(len as u32).to_le_bytes());
    bytes.extend_from_slice(&[0u8; 4]);
    bytes.extend_from_slice(xml(16 + len).as_bytes());
    bytes.extend(std::iter::repeat_n(0u8, samples * 2));
    std::fs::write(path, bytes).unwrap();
}

fn kw(name: &str, value: &str) -> String {
    format!(r#"<FITSKeyword name="{name}" value="{value}" comment=""/>"#)
}

/// Every keyword name of a header text, read by the card rule.
fn fits_keys(text: &str) -> Vec<String> {
    fits_header_reader::FitsHeader::from_header_text(text)
        .keywords_in_order()
        .map(|(k, _)| k.to_string())
        .collect()
}

#[test]
fn each_fpack_fixture_s_logical_header_carries_the_original_s_keys() {
    let (format, original_text) = header_text(&fpack("u16_stars.fits")).unwrap();
    assert_eq!(format, HeaderFormat::Fits);
    let original = header_key_map(format, &original_text);
    let original_keys = fits_keys(&original_text);
    for fz in ["fpack_rice_u16.fz", "fpack_hc_u16_s4.fz", "fpack_plio_u16.fz"] {
        let (format, text) = header_text(&fpack(fz)).unwrap();
        assert_eq!(format, HeaderFormat::Fits, "{fz}");
        let logical = header_key_map(format, &text);
        for key in &original_keys {
            if key == "CHECKSUM" || key == "DATASUM" {
                continue;
            }
            assert_eq!(logical.get(key), original.get(key), "{fz}: {key}");
        }
        for key in ["BITPIX", "NAXIS1", "NAXIS2"] {
            assert_eq!(logical.integer(key), original.integer(key), "{fz}: {key}");
        }
        for key in fits_keys(&text) {
            assert!(
                !key.starts_with('Z') && !key.starts_with("TTYPE") && !key.starts_with("TFORM"),
                "{fz}: the table's or the compression's {key} is in the logical header"
            );
        }
        let bytes = read_header_key_map(&fpack(fz)).unwrap();
        for key in fits_keys(&text) {
            assert_eq!(bytes.get(&key), logical.get(&key), "{fz}: {key}, bytes against text");
        }
    }
}

#[test]
fn frame_facts_of_a_fits_and_an_xisf_file_match_their_cards() {
    let dir = tempfile::tempdir().unwrap();
    let fits = dir.path().join("a.fits");
    let mut cards = vec![card("SIMPLE", "T"), card("BITPIX", "16"), card("NAXIS", "2"), card("NAXIS1", "4"), card("NAXIS2", "2")];
    for (k, v) in [
        ("IMAGETYP", "'LIGHT'"),
        ("EXPOSURE", "120.5"),
        ("PIXSIZE1", "3.76"),
        ("BAYERPAT", "'RGGB'"),
        ("XBAYROFF", "1"),
        ("XORGSUBF", "100"),
        ("YORGSUBF", "8"),
        ("DATE-OBS", "'2026-01-02'"),
        ("TIME-OBS", "'21:30:00.5'"),
    ] {
        cards.push(card(k, v));
    }
    write_fits(&fits, &cards, &[0; 8]);
    let xisf = dir.path().join("a.xisf");
    let children = [
        kw("IMAGETYP", "'LIGHT'"),
        kw("EXPOSURE", "120.5"),
        kw("PIXSIZE1", "3.76"),
        kw("BAYERPAT", "'RGGB'"),
        kw("XBAYROFF", "1"),
        kw("XORGSUBF", "100"),
        kw("YORGSUBF", "8"),
        kw("DATE-OBS", "'2026-01-02'"),
        kw("TIME-OBS", "'21:30:00.5'"),
    ]
    .concat();
    write_xisf(&xisf, "4:2:1", &children, 8);

    for (path, format) in [(&fits, HeaderFormat::Fits), (&xisf, HeaderFormat::Xisf)] {
        let f = read_frame_facts(path).unwrap();
        let want = FrameFacts {
            format,
            width: Some(4),
            height: Some(2),
            planes: Some(1),
            bitpix: Some(16),
            imagetyp: Some("Light".into()),
            exptime_s: Some(120.5),
            gain: None,
            offset: None,
            egain: None,
            ccd_temp: None,
            set_temp: None,
            filter: None,
            object: None,
            instrume: None,
            telescop: None,
            focallen_mm: None,
            xpixsz_um: Some(3.76),
            ypixsz_um: None,
            xbinning: None,
            ybinning: None,
            date_obs: Some("2026-01-02".into()),
            time_obs: Some("21:30:00.5".into()),
            bayerpat: Some("RGGB".into()),
            xbayroff: Some(1),
            ybayroff: None,
            roworder: None,
            origin: Some((100, 8)),
            calstat: None,
        };
        assert_eq!(f, want, "{}", path.display());
    }
}

#[test]
fn the_key_map_of_a_file_equals_the_key_map_of_its_stored_text() {
    let dir = tempfile::tempdir().unwrap();
    let fits = dir.path().join("b.fits");
    let cards = vec![
        card("SIMPLE", "T"),
        card("BITPIX", "16"),
        card("NAXIS", "2"),
        card("NAXIS1", "2"),
        card("NAXIS2", "2"),
        card("EGAIN", "1.5"),
        card("EGAIN", "2.5"),
        format!("{:<80}", "HIERARCH Eso Det Temp = -7.25 / sensor"),
        card("FILTER", "'Long filter &'"),
        format!("{:<80}", "CONTINUE  '  two'"),
        card("OBJECT", "'O''Brien'"),
        format!("{:<80}", "COMMENT a comment"),
    ];
    write_fits(&fits, &cards, &[0; 4]);
    let xisf = dir.path().join("b.xisf");
    let children = [
        kw("EGAIN", "1.5"),
        kw("EGAIN", "2.5"),
        kw("FILTER", "'Ha / 7nm'"),
        kw(" xbayroff", "1"),
        kw("OBJECT", ""),
        kw("DATAMAX", "'0.9'"),
    ]
    .concat();
    write_xisf(&xisf, "2:2:1", &children, 4);

    let keys = ["EGAIN", "ESO DET TEMP", "FILTER", "OBJECT", "COMMENT", "XBAYROFF", "DATAMAX", "BITPIX", "ABSENT"];
    for path in [&fits, &xisf] {
        let bytes = read_header_key_map(path).unwrap();
        let (format, text) = header_text(path).unwrap();
        let stored = header_key_map(format, &text);
        assert_eq!(bytes.format(), stored.format());
        for key in keys {
            assert_eq!(bytes.get(key), stored.get(key), "{}: {key}", path.display());
        }
    }
    let f = read_header_key_map(&fits).unwrap();
    assert_eq!(f.number("EGAIN"), Some(1.5), "FITS: the first card wins");
    assert_eq!(f.get("FILTER"), Some("Long filter  two"));
    let x = read_header_key_map(&xisf).unwrap();
    assert_eq!(x.number("EGAIN"), Some(2.5), "XISF: the last element wins");
    assert_eq!(x.get("FILTER"), Some("Ha / 7nm"));
    assert_eq!(x.get("OBJECT"), Some(""));
}

#[test]
fn a_gzipped_fits_reads_the_same_facts_and_pixels_as_the_plain_file() {
    use std::io::Write;
    let dir = tempfile::tempdir().unwrap();
    let plain = dir.path().join("c.fits");
    let cards = vec![
        card("SIMPLE", "T"),
        card("BITPIX", "16"),
        card("NAXIS", "2"),
        card("NAXIS1", "4"),
        card("NAXIS2", "2"),
        card("BZERO", "32768"),
        card("EGAIN", "0.8"),
        card("EXPTIME", "300"),
        card("IMAGETYP", "'Dark Frame'"),
    ];
    write_fits(&plain, &cards, &[-32768, -1000, 0, 1000, 2000, 3000, 4000, 32767]);
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(&std::fs::read(&plain).unwrap()).unwrap();
    let gzipped = dir.path().join("c.fits.gz");
    std::fs::write(&gzipped, gz.finish().unwrap()).unwrap();

    let want = read_frame_facts(&plain).unwrap();
    assert_eq!(want.imagetyp.as_deref(), Some("Dark"));
    assert_eq!(want.egain, Some(0.8));
    assert_eq!(read_frame_facts(&gzipped).unwrap(), want);
    assert_eq!(header_text(&gzipped).unwrap(), header_text(&plain).unwrap());
    let pixels = |p: &Path| match crate::ImageConverter::read_raw(p).unwrap().1 {
        crate::PixelData::Uint16(v) => v,
        crate::PixelData::Float32(_) => panic!("{}: expected u16", p.display()),
    };
    assert_eq!(pixels(&gzipped), pixels(&plain));
}

#[test]
fn an_xisf_without_keywords_gives_athenaeum_s_placeholder_text() {
    let dir = tempfile::tempdir().unwrap();
    let xisf = dir.path().join("d.xisf");
    write_xisf(&xisf, "2:2:1", "", 4);
    let (format, text) = header_text(&xisf).unwrap();
    assert_eq!(format, HeaderFormat::Xisf);
    assert_eq!(text, format!("Header extracted from: {}\n(No FITS keywords found in XISF)", xisf.display()));
    assert!(!header_key_map(format, &text).has_keywords());
    assert!(!read_header_key_map(&xisf).unwrap().has_keywords());
}

// ---- the pinned inputs of the header rules, moved here from `formats::stored_keys` when it was removed ----
//
// Each pin read through the crate's older copy of Athenaeum's Analysis rules (`stored_keys`, Athenaeum at e6344ddd)
// before starloom phase 2; it reads through `header` now, the catalog's own reader (P2-R1, P2-R23). A comment names
// the pins whose value changed (the owner saw them at checkpoint 3, 2026-11-24, and ruled 1 and 2 on them); every
// other expectation is the old one.

use super::facts::fits_bayer_offsets;
use super::text::xisf_keyword_pairs;

/// An XISF header holding `children` (the `<FITSKeyword>`s) in an `<Image>`.
fn xisf_xml(children: &str) -> String {
    format!(
        r#"<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0"><Image geometry="2:2:1">{children}</Image></xisf>"#
    )
}

/// `XBAYROFF` of an XISF header holding `children`, as the catalog's column holds it.
fn xisf_x(children: &str) -> Option<i64> {
    xisf_bayer_offsets(&xisf_xml(children)).unwrap().0
}

/// What the catalog's reader returns for `XBAYROFF` on the shapes the Analysis rules read otherwise.
#[test]
fn xisf_rule_pins_the_catalogs_results() {
    // An empty value in front is no value: the later keyword of the name wins.
    assert_eq!(xisf_x(&format!("{}{}", kw("XBAYROFF", ""), kw("XBAYROFF", "1"))), Some(1));
    assert_eq!(xisf_x(&format!("{}{}", kw("XBAYROFF", "''"), kw("XBAYROFF", "1"))), Some(1));
    assert_eq!(xisf_x(&format!(r#"<FITSKeyword name="XBAYROFF"/>{}"#, kw("XBAYROFF", "1"))), Some(1), "no value attribute");
    // The value is the raw attribute: a comment after the number makes it no number. Changed: was 1, 1 (the Analysis
    // rule cut the value at " / "): sk-x04, rm-02, xu-04; sk-x05, xu-05.
    assert_eq!(xisf_x(&kw("XBAYROFF", "1 / off")), None);
    assert_eq!(xisf_x(&kw("XBAYROFF", "'1' / off")), None);
    // An unclosed quote is no number.
    assert_eq!(xisf_x(&kw("XBAYROFF", "'1")), None);
    // The name is the one written. Changed: was 1, 1 (the Analysis rule trimmed and upper-cased the name): sk-x07,
    // sk-x08, xu-09, xv-13.
    assert_eq!(xisf_x(&kw(" XBAYROFF", "1")), None);
    assert_eq!(xisf_x(&kw("xbayroff", "1")), None);
    // The raw attribute is read, a character reference stays text.
    assert_eq!(xisf_x(&kw("XBAYROFF", "&#49;")), None);
    // The last element of a name wins. Changed: was 0 (the first non-empty value won): sk-x10, xv-15; xv-14 (1 then 0).
    assert_eq!(xisf_x(&format!("{}{}", kw("XBAYROFF", "0"), kw("XBAYROFF", "1"))), Some(1));
    assert_eq!(xisf_x(&format!("{}{}", kw("XBAYROFF", "1"), kw("XBAYROFF", "0"))), Some(0));
    // Empty elements only: a start / end pair is not read.
    assert_eq!(xisf_x(r#"<FITSKeyword name="XBAYROFF" value="1"></FITSKeyword>"#), None);
    // A `<Property>` is not a keyword.
    assert_eq!(xisf_x(r#"<Property id="XBAYROFF" type="Int32" value="1"/>"#), None);
    // A lone quote is no pair of quotes; `'1'0'` has its outer pair stripped and `1'0` is no number.
    assert_eq!(xisf_x(&kw("XBAYROFF", "'")), None);
    assert_eq!(xisf_x(&kw("XBAYROFF", "'1'0'")), None);
    // The header is cut at the first `<?xml`: a keyword in front of it is not read.
    let before_the_cut = format!(r#"{}<?xml version="1.0"?><xisf><Image geometry="2:2:1"/></xisf>"#, kw("XBAYROFF", "1"));
    assert_eq!(xisf_bayer_offsets(&before_the_cut).unwrap(), (None, None));
    // A keyword value that holds an element is a value, not a second header. Changed: was 1 (the Analysis rule read
    // the dump as XML first): sk-x23.
    assert_eq!(xisf_x(&kw("Z", "<FITSKeyword name='XBAYROFF' value='1'/>")), None);
}

#[test]
fn xisf_number_and_offset_read_as_the_catalog_reads_them() {
    let xml = xisf_xml(&format!(
        "{}{}{}{}{}{}",
        kw("A", "1.5"),
        kw("B", "'1.5'"),
        kw("XBAYROFF", "0.6"),
        kw("YBAYROFF", "-1"),
        kw("C", "1D0"),
        kw("D", "nan"),
    ));
    let map = header_key_map(HeaderFormat::Xisf, &super::text::xisf_header_text(xml.as_bytes(), Path::new("a.xisf")).unwrap());
    assert_eq!(map.number("A"), Some(1.5));
    assert_eq!(map.number("B"), Some(1.5));
    assert_eq!(map.number("Z"), None);
    assert_eq!(map.number("C"), Some(1.0), "a Fortran exponent");
    assert!(map.number("D").is_some_and(f64::is_nan), "NaN as written");
    // A fraction is cut, not rounded. Changed: 0.6 was 1: sk-x16, xv-06, xv-12.
    assert_eq!(xisf_bayer_offsets(&xml).unwrap(), (Some(0), Some(-1)));
    for (value, want) in [("1D0", None), ("nan", Some(0)), ("1.5", Some(1)), ("-0.6", Some(0)), ("x", None), ("", None)] {
        assert_eq!(xisf_x(&kw("XBAYROFF", value)), want, "XBAYROFF {value:?}");
    }
}

#[test]
fn an_xisf_without_keywords_reads_no_offsets() {
    assert_eq!(xisf_bayer_offsets(&xisf_xml("")).unwrap(), (None, None));
    assert!(xisf_bayer_offsets("").is_err());
    assert!(xisf_bayer_offsets("not xml at all").is_err());
}

/// The `<?xml` … `</xisf>` cut: a keyword after the closing tag is not read.
#[test]
fn a_keyword_after_the_closing_tag_is_not_read() {
    let xml = format!(r#"junk<?xml version="1.0"?><xisf>{}</xisf>{}"#, kw("EGAIN", "1.5"), kw("DATAMAX", "9"));
    let pairs = xisf_keyword_pairs(xml.as_bytes(), None).unwrap();
    assert_eq!(pairs, vec![("EGAIN".to_string(), "1.5".to_string())]);
}

/// The keyword pass's parse-abort line names the file when the caller has it, and has no `path` field when it has
/// not. A `</xisf>` inside a comment ends the cut section in the comment.
#[test]
fn the_parse_abort_log_line_carries_the_path_when_given() {
    let xml = format!(r#"<?xml version="1.0"?><xisf><!-- </xisf> -->{}</xisf>"#, kw("EGAIN", "2"));
    let path = Path::new("/data/a_frame.xisf");
    let mut pairs = Vec::new();
    let with = crate::test_support::captured_logs(|| pairs = xisf_keyword_pairs(xml.as_bytes(), Some(path)).unwrap());
    assert!(with.contains("error parsing XISF XML header; using partial result"), "{with}");
    assert!(with.contains("path=/data/a_frame.xisf"), "{with}");
    assert!(pairs.is_empty(), "the keyword is past the cut: the partial result is empty");
    let without = crate::test_support::captured_logs(|| pairs = xisf_keyword_pairs(xml.as_bytes(), None).unwrap());
    assert!(without.contains("error parsing stored XISF XML header; using partial result"), "{without}");
    assert!(!without.contains("path="), "{without}");
}

#[test]
fn xisf_keywords_after_a_dangling_ampersand_are_still_read() {
    let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<xisf>
  <Property id="note">Dark & Flat calibration</Property>
  <FITSKeyword name="OBJECT" value="'M42'" comment="" />
  <FITSKeyword name="EXPTIME" value="120.0" />
</xisf>"#;
    let pairs = xisf_keyword_pairs(xml.as_bytes(), None).unwrap();
    assert_eq!(pairs, vec![("OBJECT".to_string(), "M42".to_string()), ("EXPTIME".to_string(), "120.0".to_string())]);
}

/// A capture tool's plain `KEY = value` header dump reads as the stored text of an XISF header.
#[test]
fn a_plain_keyword_dump_reads_as_an_xisf_header_text() {
    let dump = "XISF FITS Keywords:\n==================\n\nSIMPLE = T\nBITPIX = 16\nNAXIS = 2\nCOMMENT = \nFOCALLEN = 140.616\nXPIXSZ = 2.4\nOBJCTRA = 2 41 04.824\nOBJCTDEC = +62 56 44.95\nCROTA2 = 4.383973760472021\nCREATOR = a capture tool\n";
    let map = header_key_map(HeaderFormat::Xisf, dump);
    assert_eq!(map.get("FOCALLEN"), Some("140.616"));
    assert_eq!(map.get("OBJCTRA"), Some("2 41 04.824"));
    assert_eq!(map.get("OBJCTDEC"), Some("+62 56 44.95"));
}

/// 80-column FITS cards, `END`, padded to a 2880-byte block.
fn fits_header_of(cards: &[&str]) -> fits_header_reader::FitsHeader {
    let mut block = Vec::new();
    for c in cards {
        block.extend_from_slice(format!("{c:<80}").as_bytes());
    }
    block.extend_from_slice(format!("{:<80}", "END").as_bytes());
    block.resize(block.len().div_ceil(2880) * 2880, b' ');
    fits_header_reader::FitsHeader::from_bytes(&block).unwrap()
}

/// What the catalog's `FitsHeader` returns for `EGAIN` on card shapes a plain card parser reads otherwise.
#[test]
fn fits_rule_pins_the_catalogs_results() {
    let egain = |cards: &[&str]| fits_header_of(cards).get_f64("EGAIN");
    assert_eq!(egain(&["EGAIN    1.5"]), None);
    assert_eq!(egain(&["EGAIN   =1.5"]), None);
    assert_eq!(egain(&["EGAIN", "EGAIN   = 1.5"]), Some(1.5));
    assert_eq!(egain(&["HIERARCH EGAIN = 1.5"]), Some(1.5));
    assert_eq!(egain(&[" EGAIN  = 1.5"]), Some(1.5));
    assert_eq!(egain(&["EGAIN   = '1.5&'"]), Some(1.5));
    assert_eq!(egain(&["EGAIN   = '1.&'", "CONTINUE  '5'"]), Some(1.5));
    assert_eq!(egain(&["EGAIN   = 1.5D0 / e-/ADU"]), Some(1.5));
    assert_eq!(egain(&["EGAIN   = 1.5", "EGAIN   = 2.5"]), Some(1.5));
    assert_eq!(egain(&["EGAIN   = '1.5'"]), Some(1.5));
    assert_eq!(fits_header_of(&["BITPIX  = 16.0"]).get_i32("BITPIX"), Some(16));
}

/// The FITS offsets are the catalog column's: `get_i32` on the key as written, the first card winning.
#[test]
fn fits_bayer_offsets_are_the_catalogs_columns() {
    let x = |cards: &[&str]| fits_bayer_offsets(&fits_header_of(cards)).0;
    assert_eq!(x(&["XBAYROFF= 1"]), Some(1));
    assert_eq!(x(&["XBAYROFF= '1'"]), Some(1));
    assert_eq!(x(&["XBAYROFF= 0", "XBAYROFF= 1"]), Some(0), "the first card wins");
    assert_eq!(x(&["XBAYROFF= 0.6"]), Some(0), "a fraction is cut");
    assert_eq!(x(&["XBAYROFF= 1D0"]), Some(1), "a Fortran exponent");
    assert_eq!(x(&["XBAYROFF= 1 / comment"]), Some(1));
    assert_eq!(x(&[]), None);
}
