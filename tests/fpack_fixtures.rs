//! Tile-compressed images written by the reference `fpack` (see `tests/fpack/README.md`), decoded and compared
//! pixel by pixel with the image they hold, stored plain and read by the in-place reader.
use astroimage::{ImageConverter, ImageMetadata, PixelData};

const DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fpack");

fn read_u16(name: &str) -> (ImageMetadata, Vec<u16>) {
    let path = format!("{DIR}/{name}");
    let (meta, pixels) =
        ImageConverter::read_raw(&path).unwrap_or_else(|e| panic!("{name}: {e:#}"));
    let PixelData::Uint16(v) = pixels else {
        panic!("{name}: expected u16 pixels")
    };
    (meta, v)
}

/// `fz` decodes to exactly the pixels and layout of the plain file `plain`.
fn assert_decodes_to(fz: &str, plain: &str) -> Vec<u16> {
    let (want_meta, want) = read_u16(plain);
    let (meta, got) = read_u16(fz);
    assert_eq!(
        (
            meta.width,
            meta.height,
            meta.channels,
            meta.dtype,
            meta.bayer_pattern,
            meta.flip_vertical
        ),
        (
            want_meta.width,
            want_meta.height,
            want_meta.channels,
            want_meta.dtype,
            want_meta.bayer_pattern,
            want_meta.flip_vertical
        ),
        "{fz} layout"
    );
    let differing = got.iter().zip(&want).filter(|(a, b)| a != b).count();
    assert_eq!(differing, 0, "{fz}: {differing} pixels differ from {plain}");
    got
}

#[test]
fn fpack_rice_u16_is_lossless() {
    assert_decodes_to("fpack_rice_u16.fz", "u16_stars.fits");
}

/// Unsigned 16-bit PLIO: the writer offsets the samples by 32768, which the decode takes off again.
#[test]
fn fpack_plio_u16_is_lossless() {
    assert_decodes_to("fpack_plio_u16.fz", "u16_stars.fits");
}

/// Lossy HCOMPRESS next to saturated stars decodes slightly above 65535; those values clip to 65535 rather than
/// wrapping to near-black.
#[test]
fn fpack_hcompress_u16_clips_saturated_stars() {
    let got = assert_decodes_to("fpack_hc_u16_s4.fz", "fpack_hc_u16_s4.expected.fits");
    let (_, original) = read_u16("u16_stars.fits");
    let max_error = got
        .iter()
        .zip(&original)
        .map(|(&a, &b)| a.abs_diff(b))
        .max()
        .unwrap();
    assert!(
        max_error <= 1558,
        "lossy error {max_error} against the original"
    );
    assert_eq!(
        got.iter().filter(|&&v| v == u16::MAX).count(),
        213,
        "saturated pixels"
    );
}

#[test]
fn fpack_header_is_the_images() {
    for name in [
        "fpack_rice_u16.fz",
        "fpack_hc_u16_s4.fz",
        "fpack_plio_u16.fz",
    ] {
        let h = ImageConverter::read_header(format!("{DIR}/{name}")).unwrap();
        assert_eq!(
            (h.width, h.height, h.channels, h.sample_type),
            (64, 64, 1, "u16"),
            "{name}"
        );
        assert_eq!(h.get("BITPIX"), Some("16"), "{name}");
        assert_eq!(h.get("BZERO"), Some("32768"), "{name}");
        assert!(
            h.get("ZCMPTYPE").is_none() && h.get("TFORM1").is_none(),
            "{name}"
        );
    }
}
