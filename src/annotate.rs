/// Star annotation overlay: draws detected-star ellipses onto converted images.

use crate::analysis::AnalysisResult;
use crate::types::ProcessedImage;

/// Color scheme for star annotations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ColorScheme {
    /// All annotations use a single color (green).
    Uniform,
    /// Color by eccentricity: green (round) → yellow → red (elongated).
    Eccentricity,
    /// Color by FWHM relative to median: green (tight) → yellow → red (bloated).
    Fwhm,
}

/// Configuration for annotation rendering.
pub struct AnnotationConfig {
    /// Color scheme for annotations.
    pub color_scheme: ColorScheme,
    /// Draw a direction tick along the elongation axis.
    pub show_direction_tick: bool,
    /// Minimum ellipse semi-axis radius in output pixels.
    pub min_radius: f32,
    /// Maximum ellipse semi-axis radius in output pixels.
    pub max_radius: f32,
    /// Line thickness: 1 = single pixel, 2 = 3px cross kernel, 3 = 5px diamond.
    pub line_width: u8,
    /// Eccentricity threshold: below this is green (good).
    pub ecc_good: f32,
    /// Eccentricity threshold: between `ecc_good` and this is yellow (warning).
    /// At or above this is red (problem).
    pub ecc_warn: f32,
    /// FWHM ratio threshold: below this is green (good). Ratio = star FWHM / median FWHM.
    pub fwhm_good: f32,
    /// FWHM ratio threshold: between `fwhm_good` and this is yellow (warning).
    /// At or above this is red (problem).
    pub fwhm_warn: f32,
}

impl Default for AnnotationConfig {
    fn default() -> Self {
        AnnotationConfig {
            color_scheme: ColorScheme::Eccentricity,
            show_direction_tick: true,
            min_radius: 6.0,
            max_radius: 60.0,
            line_width: 2,
            ecc_good: 0.5,
            ecc_warn: 0.6,
            fwhm_good: 1.3,
            fwhm_warn: 2.0,
        }
    }
}

/// Pre-computed annotation for one star, in output image coordinates.
pub struct StarAnnotation {
    /// Centroid X in output image coordinates.
    pub x: f32,
    /// Centroid Y in output image coordinates.
    pub y: f32,
    /// Semi-major axis in output pixels.
    pub semi_major: f32,
    /// Semi-minor axis in output pixels.
    pub semi_minor: f32,
    /// Rotation angle (radians), counter-clockwise from +X axis.
    pub theta: f32,
    /// Original eccentricity value.
    pub eccentricity: f32,
    /// Original geometric mean FWHM (analysis pixels).
    pub fwhm: f32,
    /// RGB color based on the chosen scheme.
    pub color: [u8; 3],
}

// ── Tier 1: Raw geometry ──

/// Compute annotation geometry for all detected stars, transformed to output coordinates.
pub fn compute_annotations(
    result: &AnalysisResult,
    output_width: usize,
    output_height: usize,
    flip_vertical: bool,
    config: &AnnotationConfig,
) -> Vec<StarAnnotation> {
    if result.stars.is_empty() || result.width == 0 || result.height == 0 {
        return Vec::new();
    }

    let scale_x = output_width as f32 / result.width as f32;
    let scale_y = output_height as f32 / result.height as f32;

    result
        .stars
        .iter()
        .map(|star| {
            let x_out = star.x * scale_x;
            let y_out = if flip_vertical {
                output_height as f32 - 1.0 - star.y * scale_y
            } else {
                star.y * scale_y
            };

            // Semi-axes: scale FWHM to output, multiply by 2.5 for visibility, then clamp.
            // Use the larger of (fwhm_x, fwhm_y) as major and smaller as minor,
            // preserving the axis ratio so elongation is clearly visible.
            let (raw_a, raw_b) = if star.fwhm_x >= star.fwhm_y {
                (star.fwhm_x * scale_x, star.fwhm_y * scale_y)
            } else {
                (star.fwhm_y * scale_y, star.fwhm_x * scale_x)
            };
            let semi_major = (raw_a * 2.5).clamp(config.min_radius, config.max_radius);
            let semi_minor = (raw_b * 2.5).clamp(config.min_radius, config.max_radius);

            let color = star_color(config, star.eccentricity, star.fwhm, result.median_fwhm);

            StarAnnotation {
                x: x_out,
                y: y_out,
                semi_major,
                semi_minor,
                theta: star.theta,
                eccentricity: star.eccentricity,
                fwhm: star.fwhm,
                color,
            }
        })
        .collect()
}

// ── Tier 2: RGBA overlay layer ──

/// Rasterize annotations into a transparent RGBA buffer (same dimensions as output image).
pub fn create_annotation_layer(
    result: &AnalysisResult,
    output_width: usize,
    output_height: usize,
    flip_vertical: bool,
    config: &AnnotationConfig,
) -> Vec<u8> {
    let mut layer = vec![0u8; output_width * output_height * 4];
    let annotations = compute_annotations(result, output_width, output_height, flip_vertical, config);
    let lw = config.line_width;

    for ann in &annotations {
        draw_ellipse_rgba(&mut layer, output_width, output_height, ann, lw);
        if config.show_direction_tick && ann.eccentricity > 0.15 {
            draw_direction_tick_rgba(&mut layer, output_width, output_height, ann, lw);
        }
    }

    layer
}

// ── Tier 3: Burn into ProcessedImage ──

/// Draw star annotations directly onto a ProcessedImage (RGB or RGBA).
pub fn annotate_image(
    image: &mut ProcessedImage,
    result: &AnalysisResult,
    config: &AnnotationConfig,
) {
    let annotations = compute_annotations(
        result,
        image.width,
        image.height,
        image.flip_vertical,
        config,
    );
    let bpp = image.channels as usize;
    let lw = config.line_width;

    for ann in &annotations {
        draw_ellipse_rgb(&mut image.data, image.width, image.height, bpp, ann, lw);
        if config.show_direction_tick && ann.eccentricity > 0.15 {
            draw_direction_tick_rgb(&mut image.data, image.width, image.height, bpp, ann, lw);
        }
    }
}

// ── Drawing primitives (private) ──

/// Bounds-checked single pixel write on an RGB/RGBA buffer.
#[inline]
fn set_pixel_one(buf: &mut [u8], width: usize, height: usize, bpp: usize, x: i32, y: i32, color: [u8; 3]) {
    if x >= 0 && y >= 0 && (x as usize) < width && (y as usize) < height {
        let idx = (y as usize * width + x as usize) * bpp;
        buf[idx] = color[0];
        buf[idx + 1] = color[1];
        buf[idx + 2] = color[2];
    }
}

/// Bounds-checked single pixel write on an RGBA buffer (sets alpha to 255).
#[inline]
fn set_pixel_one_rgba(buf: &mut [u8], width: usize, height: usize, x: i32, y: i32, color: [u8; 3]) {
    if x >= 0 && y >= 0 && (x as usize) < width && (y as usize) < height {
        let idx = (y as usize * width + x as usize) * 4;
        buf[idx] = color[0];
        buf[idx + 1] = color[1];
        buf[idx + 2] = color[2];
        buf[idx + 3] = 255;
    }
}

/// Thick pixel write: draws a kernel around (x,y).
/// lw=1: single pixel, lw=2: 3px cross (+), lw>=3: 5px diamond.
#[inline]
fn set_pixel(buf: &mut [u8], width: usize, height: usize, bpp: usize, x: i32, y: i32, color: [u8; 3], lw: u8) {
    set_pixel_one(buf, width, height, bpp, x, y, color);
    if lw >= 2 {
        set_pixel_one(buf, width, height, bpp, x - 1, y, color);
        set_pixel_one(buf, width, height, bpp, x + 1, y, color);
        set_pixel_one(buf, width, height, bpp, x, y - 1, color);
        set_pixel_one(buf, width, height, bpp, x, y + 1, color);
    }
    if lw >= 3 {
        set_pixel_one(buf, width, height, bpp, x - 2, y, color);
        set_pixel_one(buf, width, height, bpp, x + 2, y, color);
        set_pixel_one(buf, width, height, bpp, x, y - 2, color);
        set_pixel_one(buf, width, height, bpp, x, y + 2, color);
    }
}

/// Thick pixel write on RGBA buffer.
#[inline]
fn set_pixel_rgba(buf: &mut [u8], width: usize, height: usize, x: i32, y: i32, color: [u8; 3], lw: u8) {
    set_pixel_one_rgba(buf, width, height, x, y, color);
    if lw >= 2 {
        set_pixel_one_rgba(buf, width, height, x - 1, y, color);
        set_pixel_one_rgba(buf, width, height, x + 1, y, color);
        set_pixel_one_rgba(buf, width, height, x, y - 1, color);
        set_pixel_one_rgba(buf, width, height, x, y + 1, color);
    }
    if lw >= 3 {
        set_pixel_one_rgba(buf, width, height, x - 2, y, color);
        set_pixel_one_rgba(buf, width, height, x + 2, y, color);
        set_pixel_one_rgba(buf, width, height, x, y - 2, color);
        set_pixel_one_rgba(buf, width, height, x, y + 2, color);
    }
}

/// Bresenham line drawing on an RGB/RGBA buffer with thickness.
fn draw_line(buf: &mut [u8], width: usize, height: usize, bpp: usize, x0: i32, y0: i32, x1: i32, y1: i32, color: [u8; 3], lw: u8) {
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    let mut x = x0;
    let mut y = y0;

    loop {
        set_pixel(buf, width, height, bpp, x, y, color, lw);
        if x == x1 && y == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            if x == x1 { break; }
            err += dy;
            x += sx;
        }
        if e2 <= dx {
            if y == y1 { break; }
            err += dx;
            y += sy;
        }
    }
}

/// Bresenham line drawing on an RGBA buffer with thickness.
fn draw_line_rgba(buf: &mut [u8], width: usize, height: usize, x0: i32, y0: i32, x1: i32, y1: i32, color: [u8; 3], lw: u8) {
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    let mut x = x0;
    let mut y = y0;

    loop {
        set_pixel_rgba(buf, width, height, x, y, color, lw);
        if x == x1 && y == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            if x == x1 { break; }
            err += dy;
            x += sx;
        }
        if e2 <= dx {
            if y == y1 { break; }
            err += dx;
            y += sy;
        }
    }
}

/// Draw a rotated ellipse by sampling parametric points and connecting with lines.
fn draw_ellipse_rgb(buf: &mut [u8], width: usize, height: usize, bpp: usize, ann: &StarAnnotation, lw: u8) {
    let steps = 64;
    let (ct, st) = (ann.theta.cos(), ann.theta.sin());
    let mut prev_x = 0i32;
    let mut prev_y = 0i32;

    for i in 0..=steps {
        let t = (i as f32) * std::f32::consts::TAU / (steps as f32);
        let ex = ann.semi_major * t.cos();
        let ey = ann.semi_minor * t.sin();
        let rx = ex * ct - ey * st + ann.x;
        let ry = ex * st + ey * ct + ann.y;
        let px = rx.round() as i32;
        let py = ry.round() as i32;

        if i > 0 {
            draw_line(buf, width, height, bpp, prev_x, prev_y, px, py, ann.color, lw);
        }
        prev_x = px;
        prev_y = py;
    }
}

/// Draw a rotated ellipse on an RGBA layer buffer.
fn draw_ellipse_rgba(buf: &mut [u8], width: usize, height: usize, ann: &StarAnnotation, lw: u8) {
    let steps = 64;
    let (ct, st) = (ann.theta.cos(), ann.theta.sin());
    let mut prev_x = 0i32;
    let mut prev_y = 0i32;

    for i in 0..=steps {
        let t = (i as f32) * std::f32::consts::TAU / (steps as f32);
        let ex = ann.semi_major * t.cos();
        let ey = ann.semi_minor * t.sin();
        let rx = ex * ct - ey * st + ann.x;
        let ry = ex * st + ey * ct + ann.y;
        let px = rx.round() as i32;
        let py = ry.round() as i32;

        if i > 0 {
            draw_line_rgba(buf, width, height, prev_x, prev_y, px, py, ann.color, lw);
        }
        prev_x = px;
        prev_y = py;
    }
}

/// Draw a direction tick extending from the ellipse edge along theta.
fn draw_direction_tick_rgb(buf: &mut [u8], width: usize, height: usize, bpp: usize, ann: &StarAnnotation, lw: u8) {
    let tick_len = ann.semi_major * ann.eccentricity * 1.2;
    if tick_len < 2.0 {
        return;
    }
    let (ct, st) = (ann.theta.cos(), ann.theta.sin());

    // Start at the ellipse edge along major axis, extend outward
    let start_x = ann.x + ann.semi_major * ct;
    let start_y = ann.y + ann.semi_major * st;
    let end_x = start_x + tick_len * ct;
    let end_y = start_y + tick_len * st;

    draw_line(buf, width, height, bpp,
        start_x.round() as i32, start_y.round() as i32,
        end_x.round() as i32, end_y.round() as i32,
        ann.color, lw);

    // Opposite side
    let start_x2 = ann.x - ann.semi_major * ct;
    let start_y2 = ann.y - ann.semi_major * st;
    let end_x2 = start_x2 - tick_len * ct;
    let end_y2 = start_y2 - tick_len * st;

    draw_line(buf, width, height, bpp,
        start_x2.round() as i32, start_y2.round() as i32,
        end_x2.round() as i32, end_y2.round() as i32,
        ann.color, lw);
}

/// Draw a direction tick on an RGBA layer buffer.
fn draw_direction_tick_rgba(buf: &mut [u8], width: usize, height: usize, ann: &StarAnnotation, lw: u8) {
    let tick_len = ann.semi_major * ann.eccentricity * 1.2;
    if tick_len < 2.0 {
        return;
    }
    let (ct, st) = (ann.theta.cos(), ann.theta.sin());

    let start_x = ann.x + ann.semi_major * ct;
    let start_y = ann.y + ann.semi_major * st;
    let end_x = start_x + tick_len * ct;
    let end_y = start_y + tick_len * st;

    draw_line_rgba(buf, width, height,
        start_x.round() as i32, start_y.round() as i32,
        end_x.round() as i32, end_y.round() as i32,
        ann.color, lw);

    let start_x2 = ann.x - ann.semi_major * ct;
    let start_y2 = ann.y - ann.semi_major * st;
    let end_x2 = start_x2 - tick_len * ct;
    let end_y2 = start_y2 - tick_len * st;

    draw_line_rgba(buf, width, height,
        start_x2.round() as i32, start_y2.round() as i32,
        end_x2.round() as i32, end_y2.round() as i32,
        ann.color, lw);
}

/// Choose annotation color based on the color scheme and star metrics.
fn star_color(config: &AnnotationConfig, eccentricity: f32, fwhm: f32, median_fwhm: f32) -> [u8; 3] {
    match config.color_scheme {
        ColorScheme::Uniform => [0, 255, 0],
        ColorScheme::Eccentricity => {
            if eccentricity <= config.ecc_good {
                [0, 255, 0]       // Green: round, good
            } else if eccentricity <= config.ecc_warn {
                [255, 255, 0]     // Yellow: slightly elongated
            } else {
                [255, 64, 64]     // Red: problem
            }
        }
        ColorScheme::Fwhm => {
            if median_fwhm <= 0.0 {
                return [0, 255, 0];
            }
            let ratio = fwhm / median_fwhm;
            if ratio < config.fwhm_good {
                [0, 255, 0]       // Green: tight
            } else if ratio < config.fwhm_warn {
                [255, 255, 0]     // Yellow: somewhat bloated
            } else {
                [255, 64, 64]     // Red: very bloated
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::{AnalysisResult, StarMetrics};

    fn dummy_result(stars: Vec<StarMetrics>) -> AnalysisResult {
        AnalysisResult {
            width: 100,
            height: 100,
            source_channels: 1,
            background: 0.0,
            noise: 0.0,
            detection_threshold: 0.0,
            stars_detected: stars.len(),
            median_fwhm: 5.0,
            median_eccentricity: 0.2,
            median_snr: 50.0,
            median_hfr: 3.0,
            snr_db: 20.0,
            snr_weight: 100.0,
            psf_signal: 50.0,
            trail_r_squared: 0.0,
            possibly_trailed: false,
            stars,
        }
    }

    fn make_star(x: f32, y: f32, fwhm: f32, ecc: f32) -> StarMetrics {
        StarMetrics {
            x, y,
            peak: 1000.0,
            flux: 5000.0,
            fwhm_x: fwhm,
            fwhm_y: fwhm * (1.0 - ecc * ecc).sqrt(),
            fwhm,
            eccentricity: ecc,
            snr: 50.0,
            hfr: fwhm * 0.6,
            theta: 0.0,
        }
    }

    #[test]
    fn test_compute_annotations_empty() {
        let result = dummy_result(vec![]);
        let anns = compute_annotations(&result, 100, 100, false, &AnnotationConfig::default());
        assert!(anns.is_empty());
    }

    #[test]
    fn test_compute_annotations_coordinate_transform() {
        let star = make_star(50.0, 25.0, 5.0, 0.1);
        let result = dummy_result(vec![star]);

        // Same size, no flip
        let anns = compute_annotations(&result, 100, 100, false, &AnnotationConfig::default());
        assert_eq!(anns.len(), 1);
        assert!((anns[0].x - 50.0).abs() < 0.1);
        assert!((anns[0].y - 25.0).abs() < 0.1);

        // Same size, flipped
        let anns = compute_annotations(&result, 100, 100, true, &AnnotationConfig::default());
        assert!((anns[0].y - 74.0).abs() < 0.1); // 99 - 25 = 74

        // Half size (e.g. debayer)
        let anns = compute_annotations(&result, 50, 50, false, &AnnotationConfig::default());
        assert!((anns[0].x - 25.0).abs() < 0.1);
        assert!((anns[0].y - 12.5).abs() < 0.1);
    }

    #[test]
    fn test_eccentricity_colors() {
        let config = AnnotationConfig::default();
        assert_eq!(star_color(&config, 0.3, 5.0, 5.0), [0, 255, 0]);       // below 0.5 → green
        assert_eq!(star_color(&config, 0.55, 5.0, 5.0), [255, 255, 0]);    // 0.51..0.6 → yellow
        assert_eq!(star_color(&config, 0.7, 5.0, 5.0), [255, 64, 64]);     // above 0.6 → red
    }

    #[test]
    fn test_annotate_image_smoke() {
        let star = make_star(50.0, 50.0, 5.0, 0.2);
        let result = dummy_result(vec![star]);

        let mut image = ProcessedImage {
            data: vec![0u8; 100 * 100 * 3],
            width: 100,
            height: 100,
            is_color: false,
            channels: 3,
            flip_vertical: false,
        };

        annotate_image(&mut image, &result, &AnnotationConfig::default());

        // At least some pixels should have been drawn (non-zero)
        let nonzero = image.data.iter().filter(|&&b| b > 0).count();
        assert!(nonzero > 0, "Expected some drawn pixels");
    }

    #[test]
    fn test_create_annotation_layer_smoke() {
        let star = make_star(50.0, 50.0, 5.0, 0.2);
        let result = dummy_result(vec![star]);

        let layer = create_annotation_layer(&result, 100, 100, false, &AnnotationConfig::default());
        assert_eq!(layer.len(), 100 * 100 * 4);

        // Check that some alpha values are 255 (drawn pixels)
        let drawn = layer.chunks_exact(4).filter(|px| px[3] == 255).count();
        assert!(drawn > 0, "Expected some drawn pixels in layer");
    }

    #[test]
    fn test_bresenham_diagonal() {
        let mut buf = vec![0u8; 10 * 10 * 3];
        draw_line(&mut buf, 10, 10, 3, 0, 0, 9, 9, [255, 0, 0], 1);
        // Check that the diagonal has some red pixels
        let red_count = buf.chunks_exact(3).filter(|px| px[0] == 255).count();
        assert!(red_count >= 10, "Expected at least 10 red pixels on diagonal");
    }
}
