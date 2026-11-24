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
