/// Debug image rendering helpers for the debug pipeline binary.
///
/// All functions take f32 image data and produce 8-bit images via the `image` crate.

use std::path::Path;

use anyhow::Result;
use image::{GrayImage, ImageEncoder, RgbImage};

/// Auto-stretch f32 data to 8-bit grayscale PNG.
pub fn save_grayscale(data: &[f32], w: usize, h: usize, path: &Path) -> Result<()> {
    let (min, max) = min_max(data);
    let scale = if (max - min).abs() < 1e-10 {
        0.0
    } else {
        255.0 / (max - min)
    };

    let mut img = GrayImage::new(w as u32, h as u32);
    for y in 0..h {
        for x in 0..w {
            let v = ((data[y * w + x] - min) * scale).clamp(0.0, 255.0) as u8;
            img.put_pixel(x as u32, y as u32, image::Luma([v]));
        }
    }

    save_png_gray(&img, path)
}

/// Render f32 image with crosshair markers at given positions.
pub fn save_with_markers(
    data: &[f32],
    w: usize,
    h: usize,
    markers: &[(f32, f32, [u8; 3])],
    path: &Path,
) -> Result<()> {
    let img = stretch_to_rgb(data, w, h);
    let mut img = img;

    for &(mx, my, color) in markers {
        draw_crosshair(&mut img, mx, my, color, w, h);
    }

    save_png_rgb(&img, path)
}

/// Render f32 image with PSF ellipses at given positions.
pub fn save_with_ellipses(
    data: &[f32],
    w: usize,
    h: usize,
    ellipses: &[(f32, f32, f32, f32, f32, [u8; 3])],
    path: &Path,
) -> Result<()> {
    let mut img = stretch_to_rgb(data, w, h);

    for &(cx, cy, semi_a, semi_b, theta, color) in ellipses {
        draw_ellipse(&mut img, cx, cy, semi_a, semi_b, theta, color, w, h);
    }

    save_png_rgb(&img, path)
}

/// Render a small stamp upscaled 8x with centroid crosshair.
pub fn save_stamp(
    stamp: &[f32],
    w: usize,
    h: usize,
    cx: f32,
    cy: f32,
    path: &Path,
) -> Result<()> {
    let scale = 8;
    let uw = w * scale;
    let uh = h * scale;

    let (min, max) = min_max(stamp);
    let s = if (max - min).abs() < 1e-10 {
        0.0
    } else {
        255.0 / (max - min)
    };

    let mut img = RgbImage::new(uw as u32, uh as u32);
    for y in 0..h {
        for x in 0..w {
            let v = ((stamp[y * w + x] - min) * s).clamp(0.0, 255.0) as u8;
            for dy in 0..scale {
                for dx in 0..scale {
                    img.put_pixel(
                        (x * scale + dx) as u32,
                        (y * scale + dy) as u32,
                        image::Rgb([v, v, v]),
                    );
                }
            }
        }
    }

    // Draw crosshair at centroid (scaled)
    draw_crosshair(&mut img, cx * scale as f32, cy * scale as f32, [0, 255, 0], uw, uh);

    save_png_rgb(&img, path)
}

/// Render side-by-side panels (original | model | residual), each upscaled 8x.
pub fn save_fit_comparison(
    original: &[f32],
    model: &[f32],
    residual: &[f32],
    w: usize,
    h: usize,
    path: &Path,
) -> Result<()> {
    let scale = 8;
    let pw = w * scale; // panel width
    let ph = h * scale;
    let total_w = pw * 3 + 4; // 2px gap between panels

    let mut img = RgbImage::new(total_w as u32, ph as u32);

    let panels: [(&[f32], usize); 3] = [
        (original, 0),
        (model, pw + 2),
        (residual, 2 * (pw + 2)),
    ];

    for (data, x_off) in panels {
        let (min, max) = min_max(data);
        let s = if (max - min).abs() < 1e-10 {
            0.0
        } else {
            255.0 / (max - min)
        };

        for y in 0..h {
            for x in 0..w {
                let v = ((data[y * w + x] - min) * s).clamp(0.0, 255.0) as u8;
                for dy in 0..scale {
                    for dx in 0..scale {
                        let px = x_off + x * scale + dx;
                        let py = y * scale + dy;
                        if px < total_w {
                            img.put_pixel(px as u32, py as u32, image::Rgb([v, v, v]));
                        }
                    }
                }
            }
        }
    }

    save_png_rgb(&img, path)
}

// ---- Internal helpers ----

fn min_max(data: &[f32]) -> (f32, f32) {
    let mut min = f32::INFINITY;
    let mut max = f32::NEG_INFINITY;
    for &v in data {
        if v.is_finite() {
            if v < min { min = v; }
            if v > max { max = v; }
        }
    }
    if !min.is_finite() { min = 0.0; }
    if !max.is_finite() { max = 1.0; }
    (min, max)
}

fn stretch_to_rgb(data: &[f32], w: usize, h: usize) -> RgbImage {
    let (min, max) = min_max(data);
    let scale = if (max - min).abs() < 1e-10 {
        0.0
    } else {
        255.0 / (max - min)
    };

    let mut img = RgbImage::new(w as u32, h as u32);
    for y in 0..h {
        for x in 0..w {
            let v = ((data[y * w + x] - min) * scale).clamp(0.0, 255.0) as u8;
            img.put_pixel(x as u32, y as u32, image::Rgb([v, v, v]));
        }
    }
    img
}

fn draw_crosshair(img: &mut RgbImage, cx: f32, cy: f32, color: [u8; 3], w: usize, h: usize) {
    let ix = cx.round() as i32;
    let iy = cy.round() as i32;
    let arm = 5;

    for d in -arm..=arm {
        // Horizontal
        let px = ix + d;
        if px >= 0 && px < w as i32 && iy >= 0 && iy < h as i32 {
            img.put_pixel(px as u32, iy as u32, image::Rgb(color));
        }
        // Vertical
        let py = iy + d;
        if ix >= 0 && ix < w as i32 && py >= 0 && py < h as i32 {
            img.put_pixel(ix as u32, py as u32, image::Rgb(color));
        }
    }
}

fn draw_ellipse(
    img: &mut RgbImage,
    cx: f32,
    cy: f32,
    semi_a: f32,
    semi_b: f32,
    theta: f32,
    color: [u8; 3],
    w: usize,
    h: usize,
) {
    let (cos_t, sin_t) = (theta.cos(), theta.sin());
    let steps = 64;
    for i in 0..steps {
        let angle = 2.0 * std::f32::consts::PI * i as f32 / steps as f32;
        let ex = semi_a * angle.cos();
        let ey = semi_b * angle.sin();
        let px = cx + ex * cos_t - ey * sin_t;
        let py = cy + ex * sin_t + ey * cos_t;
        let ix = px.round() as i32;
        let iy = py.round() as i32;
        if ix >= 0 && ix < w as i32 && iy >= 0 && iy < h as i32 {
            img.put_pixel(ix as u32, iy as u32, image::Rgb(color));
        }
    }
}

fn save_png_gray(img: &GrayImage, path: &Path) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let buf = std::io::BufWriter::new(file);
    let encoder = image::codecs::png::PngEncoder::new(buf);
    encoder.write_image(
        img.as_raw(),
        img.width(),
        img.height(),
        image::ExtendedColorType::L8,
    )?;
    Ok(())
}

fn save_png_rgb(img: &RgbImage, path: &Path) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let buf = std::io::BufWriter::new(file);
    let encoder = image::codecs::png::PngEncoder::new(buf);
    encoder.write_image(
        img.as_raw(),
        img.width(),
        img.height(),
        image::ExtendedColorType::Rgb8,
    )?;
    Ok(())
}
