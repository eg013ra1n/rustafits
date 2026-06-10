//! End-to-end stretch robustness: float frames with NaN borders (the
//! registration-crop case) must render visibly, not all-black.

use astroimage::{BayerPattern, DataType, ImageConverter, ImageMetadata, PixelData};

fn meta(width: usize, height: usize, dtype: DataType) -> ImageMetadata {
    ImageMetadata {
        width,
        height,
        channels: 1,
        dtype,
        bayer_pattern: BayerPattern::None,
        flip_vertical: false,
    }
}

/// Deterministic pseudo-noise background (no rand dependency).
fn synth_frame(width: usize, height: usize, level: f32, spread: f32) -> Vec<f32> {
    (0..width * height)
        .map(|i| {
            let t = (i as f32 * 0.731).sin() * (i as f32 * 0.0137).cos();
            level + t * spread
        })
        .collect()
}

#[test]
fn nan_border_float_frame_renders_non_black() {
    let (w, h) = (128, 128);
    let mut data = synth_frame(w, h, 500.0, 30.0);

    // NaN out a wide border (~58% of pixels) — mimics a registered/cropped
    // float frame. Pre-fix, this poisoned the stretch MADN and the whole
    // preview rendered black.
    let border = 24;
    for y in 0..h {
        for x in 0..w {
            if x < border || x >= w - border || y < border || y >= h - border {
                data[y * w + x] = f32::NAN;
            }
        }
    }

    let img = ImageConverter::new()
        .process_data(meta(w, h, DataType::Float32), PixelData::Float32(data))
        .expect("processing must succeed");

    // Center pixels (real data) must be visible…
    let center_idx = ((h / 2) * w + w / 2) * img.channels as usize;
    let center_val = img.data[center_idx];
    assert!(
        center_val > 10,
        "center of frame must render visibly, got {center_val}"
    );
    // …and the NaN border must render black, not garbage.
    let corner_val = img.data[0];
    assert_eq!(corner_val, 0, "NaN border must render black");
}

#[test]
fn all_nan_float_frame_renders_black_without_error() {
    let (w, h) = (64, 64);
    let data = vec![f32::NAN; w * h];
    let img = ImageConverter::new()
        .process_data(meta(w, h, DataType::Float32), PixelData::Float32(data))
        .expect("fully-blank frame must not error");
    assert!(img.data.iter().all(|&v| v == 0), "all-NaN frame renders black");
}

#[test]
fn downscale_zero_is_clamped() {
    let (w, h) = (64, 64);
    let data = synth_frame(w, h, 500.0, 30.0);
    // Pre-fix this divided by zero inside the downscale kernel.
    let img = ImageConverter::new()
        .with_downscale(0)
        .process_data(meta(w, h, DataType::Float32), PixelData::Float32(data))
        .expect("downscale 0 must be treated as 1");
    assert_eq!(img.width, w);
    assert_eq!(img.height, h);
}
