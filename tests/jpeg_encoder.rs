//! Contract tests for `encode_jpeg`, backend-independent.
//!
//! Everything here decodes with `jpeg-decoder`, an *independent* pure-Rust
//! decoder, rather than with the encoder's own decode side — a symmetric bug in
//! the encoder would sail straight through a self-round-trip.
//!
//! The sharp edge this locks down: at 4:2:0 the encoder pads the trailing MCU
//! with synthetic samples, and different libjpeg implementations fill that
//! padding differently (see `edge_mcu_padding_never_leaks_into_the_image`).
//! Those samples must never influence a pixel that is actually inside the
//! image.

use astroimage::encode_jpeg;

/// Deterministic starfield-ish content: dark background, scattered gaussian
/// blobs, and a deliberately bright right/bottom border so edge damage shows up
/// as a large error rather than a subtle one.
fn synth(w: usize, h: usize, channels: usize) -> Vec<u8> {
    let mut v = vec![0u8; w * h * channels];
    let mut s: u32 = 0x9e37_79b9;
    let mut rng = move || {
        s ^= s << 13;
        s ^= s >> 17;
        s ^= s << 5;
        s
    };

    for px in v.chunks_exact_mut(channels) {
        let n = (rng() % 14) as u8 + 6;
        px[0] = n;
        px[1] = n;
        px[2] = n;
        if channels == 4 {
            px[3] = 255;
        }
    }

    let star_count = (w * h / 900).max(8);
    for _ in 0..star_count {
        let cx = (rng() as usize) % w;
        let cy = (rng() as usize) % h;
        let peak = 140 + (rng() % 110);
        for dy in 0..7usize {
            for dx in 0..7usize {
                let (x, y) = (cx + dx, cy + dy);
                if x >= w || y >= h {
                    continue;
                }
                let (fx, fy) = (dx as f32 - 3.0, dy as f32 - 3.0);
                let g = (-(fx * fx + fy * fy) / 3.0).exp();
                let val = (peak as f32 * g) as u32;
                let i = (y * w + x) * channels;
                for c in 0..3 {
                    v[i + c] = v[i + c].saturating_add(val.min(255) as u8);
                }
            }
        }
    }

    // Bright border: the last two columns and rows. These sit next to the
    // synthetic padding, so any bleed shows up immediately.
    for y in 0..h {
        for x in w.saturating_sub(2)..w {
            let i = (y * w + x) * channels;
            v[i] = 250;
            v[i + 1] = 230;
            v[i + 2] = 210;
        }
    }
    for y in h.saturating_sub(2)..h {
        for x in 0..w {
            let i = (y * w + x) * channels;
            v[i] = 250;
            v[i + 1] = 230;
            v[i + 2] = 210;
        }
    }
    v
}

/// Decode with `jpeg-decoder` and return interleaved RGB8.
fn decode_rgb(jpeg: &[u8], w: usize, h: usize) -> Vec<u8> {
    let mut d = jpeg_decoder::Decoder::new(std::io::Cursor::new(jpeg));
    let pixels = d.decode().expect("independent decoder rejected our JPEG");
    let info = d.info().expect("no frame info");
    assert_eq!(
        (info.width as usize, info.height as usize),
        (w, h),
        "decoded dimensions differ from the source"
    );
    assert_eq!(
        info.pixel_format,
        jpeg_decoder::PixelFormat::RGB24,
        "expected RGB24 output"
    );
    assert_eq!(pixels.len(), w * h * 3);
    pixels
}

/// Mean squared error over a caller-chosen pixel region, in RGB space.
fn region_mse(
    src: &[u8],
    channels: usize,
    dec: &[u8],
    w: usize,
    xs: impl Iterator<Item = usize> + Clone,
    ys: impl Iterator<Item = usize>,
) -> f64 {
    let mut sse = 0f64;
    let mut n = 0usize;
    for y in ys {
        for x in xs.clone() {
            for c in 0..3 {
                let a = src[(y * w + x) * channels + c] as f64;
                let b = dec[(y * w + x) * 3 + c] as f64;
                sse += (a - b) * (a - b);
                n += 1;
            }
        }
    }
    assert!(n > 0, "empty comparison region");
    sse / n as f64
}

fn psnr(mse: f64) -> f64 {
    if mse == 0.0 {
        f64::INFINITY
    } else {
        10.0 * (255.0f64 * 255.0 / mse).log10()
    }
}

/// The encoder pads the trailing MCU (16x16 at 4:2:0) with synthetic samples
/// whenever a dimension is not a multiple of 16. The format does not dictate
/// how that padding is filled, so two conforming encoders may emit different
/// (equally valid) bytes there — this test therefore asserts nothing about
/// byte-identity with any other implementation.
///
/// As it happens the pinned `libjpeg-turbo-rs` 0.8.0 *does* match C
/// libjpeg-turbo byte-for-byte at every width; that is a measured property of
/// the pin, not a contract, and it was NOT true of 0.6.x — see
/// `docs/libjpeg-turbo-rs-issue.md` (upstream issue #362) for the divergence
/// that used to exist at `width % 16` in 1..=8 and how to re-check it after a
/// version bump. Do not turn that observation into an assertion here: it would
/// fail for a legitimately different encoder.
///
/// What IS a contract, for any encoder: padding must never influence a pixel
/// inside the image. This sweeps every `width % 16` and `height % 16` residue
/// and checks the border strip specifically, at quality 95 where JPEG loss is
/// small enough that real damage cannot hide behind it.
#[test]
fn edge_mcu_padding_never_leaks_into_the_image() {
    for extra_w in 0..16usize {
        for extra_h in [0usize, 1, 7, 8, 9, 15] {
            let (w, h) = (96 + extra_w, 96 + extra_h);
            let src = synth(w, h, 3);
            let jpeg = encode_jpeg(&src, w, h, 3, 95).expect("encode failed");
            let dec = decode_rgb(&jpeg, w, h);

            // Whole image: normal JPEG loss.
            let all = psnr(region_mse(&src, 3, &dec, w, 0..w, 0..h));
            assert!(
                all > 34.0,
                "{w}x{h}: whole-image PSNR {all:.1} dB is too low"
            );

            // Right edge — the columns adjacent to the horizontal padding.
            let right = psnr(region_mse(&src, 3, &dec, w, w - 4..w, 0..h));
            assert!(
                right > 30.0,
                "{w}x{h} (w%16={}): right-edge PSNR {right:.1} dB — padding is bleeding \
                 into real pixels",
                w % 16
            );

            // Bottom edge — the rows adjacent to the vertical padding.
            let bottom = psnr(region_mse(&src, 3, &dec, w, 0..w, h - 4..h));
            assert!(
                bottom > 30.0,
                "{w}x{h} (h%16={}): bottom-edge PSNR {bottom:.1} dB — padding is bleeding \
                 into real pixels",
                h % 16
            );
        }
    }
}

/// `save_image` hands RGBA straight to the encoder instead of building an RGB
/// copy. That shortcut is only sound if the encoder discards alpha and reads
/// R,G,B from bytes 0,1,2 — assert the two paths are byte-identical, not merely
/// similar.
#[test]
fn rgba_input_encodes_identically_to_hand_stripped_rgb() {
    for (w, h) in [(96usize, 96usize), (99, 97), (104, 100), (67, 43)] {
        let rgba = synth(w, h, 4);
        let rgb: Vec<u8> = rgba
            .chunks_exact(4)
            .flat_map(|p| [p[0], p[1], p[2]])
            .collect();

        let from_rgba = encode_jpeg(&rgba, w, h, 4, 90).expect("rgba encode failed");
        let from_rgb = encode_jpeg(&rgb, w, h, 3, 90).expect("rgb encode failed");

        assert_eq!(
            from_rgba, from_rgb,
            "{w}x{h}: RGBA and hand-stripped RGB produced different JPEGs — the \
             alpha-skipping shortcut in save_image is not safe"
        );
    }
}

/// Quality must still move the size/fidelity trade-off in the expected
/// direction — a stuck or ignored quality parameter would otherwise pass every
/// other test here.
#[test]
fn quality_parameter_is_honoured() {
    let (w, h) = (192usize, 128usize);
    let src = synth(w, h, 3);

    let low = encode_jpeg(&src, w, h, 3, 40).unwrap();
    let high = encode_jpeg(&src, w, h, 3, 95).unwrap();
    assert!(
        low.len() < high.len(),
        "q40 ({} B) should be smaller than q95 ({} B)",
        low.len(),
        high.len()
    );

    let p_low = psnr(region_mse(&src, 3, &decode_rgb(&low, w, h), w, 0..w, 0..h));
    let p_high = psnr(region_mse(&src, 3, &decode_rgb(&high, w, h), w, 0..w, 0..h));
    assert!(
        p_high > p_low,
        "q95 PSNR {p_high:.1} dB should beat q40 {p_low:.1} dB"
    );
}

/// Guards the `channels` contract of `encode_jpeg` — the pipeline only ever
/// produces 3 or 4, and anything else must fail loudly rather than reinterpret
/// the buffer.
#[test]
fn rejects_unsupported_channel_counts_and_short_buffers() {
    let src = synth(32, 32, 3);

    for bad in [0usize, 1, 2, 5] {
        assert!(
            encode_jpeg(&src, 32, 32, bad, 90).is_err(),
            "channels={bad} should be rejected"
        );
    }

    let err = encode_jpeg(&src[..100], 32, 32, 3, 90).unwrap_err();
    assert!(
        err.to_string().contains("too small"),
        "expected a buffer-size error, got: {err}"
    );
}
