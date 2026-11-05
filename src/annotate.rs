/// Star annotation overlay: draws detected-star ellipses onto converted images.

use crate::analysis::AnalysisResult;
use crate::star_measure::analysis_mode::FrameSummary;
use crate::star_measure::PlaneMeasurement;
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
    /// Ellipse semi-axis scale in units of FWHM (semi-major = fwhm_x × scale).
    pub ellipse_scale: f32,
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
            ellipse_scale: 1.2,
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
#[derive(Debug, Clone, Copy, PartialEq)]
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

// ── The star list ──

/// One star of a [`StarOverlay`]: where it is and the shape of its ellipse, in source pixels.
///
/// Positions are in the file's row order (unflipped; [`compute_overlay`] flips when the render does). `fwhm_x` is the
/// MAJOR axis's FWHM, `fwhm_y` the minor's, and `theta` is the major axis's angle from the image x axis, radians
/// (counter-clockwise from +X): the convention of [`crate::StarMetrics`] and of `star_measure::StarFit`.
///
/// A plain, all-public struct on purpose: a caller can build one from its own measurements. Its field set is frozen for
/// 1.x (no `#[non_exhaustive]`), so adding a field is a 2.0 change.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OverlayStar {
    /// Centroid X, source pixels.
    pub x: f32,
    /// Centroid Y, source pixels, file row order.
    pub y: f32,
    /// FWHM along the major axis, source pixels.
    pub fwhm_x: f32,
    /// FWHM along the minor axis, source pixels.
    pub fwhm_y: f32,
    /// Position angle of the major axis, radians, counter-clockwise from +X.
    pub theta: f32,
    /// Eccentricity: 0 = round, approaching 1 = elongated.
    pub eccentricity: f32,
    /// Geometric mean FWHM, source pixels.
    pub fwhm: f32,
}

/// The stars to draw, in the pixel grid they were measured on.
///
/// `width` and `height` are the size of that grid; [`compute_overlay`] scales the positions and axes onto the output
/// image, which may be smaller (a Bayer mosaic is measured at its native size and rendered at half size). The
/// annotation functions take one of these, so the measurement ([`StarOverlay::from_measurement`]) and the legacy
/// [`AnalysisResult`] (`StarOverlay::from(&result)`) are drawn by the same code.
///
/// A plain, all-public struct on purpose, like [`OverlayStar`]; its field set is frozen for 1.x.
#[derive(Debug, Clone, PartialEq)]
pub struct StarOverlay {
    /// The stars to draw.
    pub stars: Vec<OverlayStar>,
    /// Width of the grid the stars were measured on, pixels.
    pub width: usize,
    /// Height of the grid the stars were measured on, pixels.
    pub height: usize,
    /// The frame's median FWHM, source pixels; 0.0 when the frame has none (every star then reads green under
    /// [`ColorScheme::Fwhm`]).
    pub median_fwhm: f32,
}

impl StarOverlay {
    /// The overlay of an Analysis-mode measurement: one [`OverlayStar`] per accepted fit in `pm.stars` (saturated and
    /// extended fits included, the stars the frame's numbers were computed on), converted from f64, and
    /// `summary.median_fwhm` as the frame median. `frame_w` and `frame_h` are the measured frame's size
    /// (`RawFrame::w`, `RawFrame::h`), also for a mosaic measured on its green sites.
    pub fn from_measurement(
        frame_w: usize,
        frame_h: usize,
        pm: &PlaneMeasurement,
        summary: &FrameSummary,
    ) -> StarOverlay {
        StarOverlay {
            stars: pm
                .stars
                .iter()
                .map(|fit| OverlayStar {
                    x: fit.x as f32,
                    y: fit.y as f32,
                    fwhm_x: fit.fwhm_x as f32,
                    fwhm_y: fit.fwhm_y as f32,
                    theta: fit.theta as f32,
                    eccentricity: fit.eccentricity() as f32,
                    fwhm: fit.fwhm() as f32,
                })
                .collect(),
            width: frame_w,
            height: frame_h,
            median_fwhm: summary.median_fwhm as f32,
        }
    }
}

impl From<&AnalysisResult> for StarOverlay {
    /// The legacy result's stars and median FWHM as they are.
    fn from(result: &AnalysisResult) -> StarOverlay {
        StarOverlay {
            stars: result
                .stars
                .iter()
                .map(|star| OverlayStar {
                    x: star.x,
                    y: star.y,
                    fwhm_x: star.fwhm_x,
                    fwhm_y: star.fwhm_y,
                    theta: star.theta,
                    eccentricity: star.eccentricity,
                    fwhm: star.fwhm,
                })
                .collect(),
            width: result.width,
            height: result.height,
            median_fwhm: result.median_fwhm,
        }
    }
}

// ── Tier 1: Raw geometry ──

/// Compute annotation geometry for all stars of an overlay, transformed to output coordinates.
///
/// Empty when the overlay has no stars or a zero size.
pub fn compute_overlay(
    overlay: &StarOverlay,
    output_width: usize,
    output_height: usize,
    flip_vertical: bool,
    config: &AnnotationConfig,
) -> Vec<StarAnnotation> {
    if overlay.stars.is_empty() || overlay.width == 0 || overlay.height == 0 {
        return Vec::new();
    }

    let scale_x = output_width as f32 / overlay.width as f32;
    let scale_y = output_height as f32 / overlay.height as f32;

    overlay
        .stars
        .iter()
        .map(|star| {
            let x_out = star.x * scale_x;
            let y_out = if flip_vertical {
                output_height as f32 - 1.0 - star.y * scale_y
            } else {
                star.y * scale_y
            };

            // Semi-axes: scale FWHM to output, multiply by ellipse_scale, then clamp.
            // The star list guarantees fwhm_x >= fwhm_y with theta along the major axis.
            let raw_a = star.fwhm_x * scale_x;
            let raw_b = star.fwhm_y * scale_y;
            let semi_major = (raw_a * config.ellipse_scale).clamp(config.min_radius, config.max_radius);
            let semi_minor = (raw_b * config.ellipse_scale).clamp(config.min_radius, config.max_radius);

            let color = star_color(config, star.eccentricity, star.fwhm, overlay.median_fwhm);

            // Y-flip reflects through the horizontal axis, negating the angle.
            let theta = if flip_vertical { -star.theta } else { star.theta };

            StarAnnotation {
                x: x_out,
                y: y_out,
                semi_major,
                semi_minor,
                theta,
                eccentricity: star.eccentricity,
                fwhm: star.fwhm,
                color,
            }
        })
        .collect()
}

/// Compute annotation geometry for all detected stars, transformed to output coordinates.
///
/// The legacy entry: [`compute_overlay`] on `StarOverlay::from(result)`.
pub fn compute_annotations(
    result: &AnalysisResult,
    output_width: usize,
    output_height: usize,
    flip_vertical: bool,
    config: &AnnotationConfig,
) -> Vec<StarAnnotation> {
    compute_overlay(&StarOverlay::from(result), output_width, output_height, flip_vertical, config)
}

// ── Tier 2: RGBA overlay layer ──

/// Rasterize an overlay's annotations into a transparent RGBA buffer (same dimensions as output image).
pub fn create_overlay_layer(
    overlay: &StarOverlay,
    output_width: usize,
    output_height: usize,
    flip_vertical: bool,
    config: &AnnotationConfig,
) -> Vec<u8> {
    let mut layer = vec![0u8; output_width * output_height * 4];
    let annotations = compute_overlay(overlay, output_width, output_height, flip_vertical, config);
    let lw = config.line_width;

    for ann in &annotations {
        draw_ellipse_rgba(&mut layer, output_width, output_height, ann, lw);
        if config.show_direction_tick && ann.eccentricity > 0.15 {
            draw_direction_tick_rgba(&mut layer, output_width, output_height, ann, lw);
        }
    }

    layer
}

/// Rasterize annotations into a transparent RGBA buffer (same dimensions as output image).
///
/// The legacy entry: [`create_overlay_layer`] on `StarOverlay::from(result)`.
pub fn create_annotation_layer(
    result: &AnalysisResult,
    output_width: usize,
    output_height: usize,
    flip_vertical: bool,
    config: &AnnotationConfig,
) -> Vec<u8> {
    create_overlay_layer(&StarOverlay::from(result), output_width, output_height, flip_vertical, config)
}

// ── Tier 3: Burn into ProcessedImage ──

/// Draw an overlay's stars directly onto a ProcessedImage (RGB or RGBA).
pub fn annotate_image_with(
    image: &mut ProcessedImage,
    overlay: &StarOverlay,
    config: &AnnotationConfig,
) {
    let annotations = compute_overlay(
        overlay,
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

/// Draw star annotations directly onto a ProcessedImage (RGB or RGBA).
///
/// The legacy entry: [`annotate_image_with`] on `StarOverlay::from(result)`.
pub fn annotate_image(
    image: &mut ProcessedImage,
    result: &AnalysisResult,
    config: &AnnotationConfig,
) {
    annotate_image_with(image, &StarOverlay::from(result), config);
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
    use crate::analysis::{AnalysisResult, StageTiming, StarMetrics};
    use crate::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
    use crate::star_measure::raw::{RawFrame, RawPlanes};
    use crate::star_measure::Domain;
    use crate::test_support::{add_noise, moffat_field, MoffatStar};

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
            snr_weight: 100.0,
            psf_signal: 50.0,
            frame_snr: 0.0,
            trail_r_squared: 0.0,
            possibly_trailed: false,
            measured_fwhm_kernel: 3.0,
            median_beta: None,
            pass1_detections: 0,
            calibrated_fwhm: 0.0,
            stars_measured: 0,
            moffat_count: 0,
            gaussian_count: 0,
            plate_scale: None,
            median_fwhm_arcsec: None,
            median_hfr_arcsec: None,
            stage_timing: StageTiming {
                background_ms: 0.0, detection_pass1_ms: 0.0, calibration_ms: 0.0,
                detection_pass2_ms: 0.0, measurement_ms: 0.0, snr_ms: 0.0,
                statistics_ms: 0.0, total_ms: 0.0,
            },
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
            beta: None,
            fit_method: crate::analysis::FitMethod::Gaussian,
            fit_residual: 0.0,
            fwhm_arcsec: None,
            hfr_arcsec: None,
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
    fn test_flip_vertical_negates_theta() {
        let theta = std::f32::consts::FRAC_PI_6; // 30°
        let star = StarMetrics {
            x: 50.0,
            y: 25.0,
            peak: 1000.0,
            flux: 5000.0,
            fwhm_x: 8.0,
            fwhm_y: 4.0,
            fwhm: 5.66,
            eccentricity: 0.87,
            snr: 50.0,
            hfr: 3.0,
            theta,
            beta: None,
            fit_method: crate::analysis::FitMethod::Gaussian,
            fit_residual: 0.0,
            fwhm_arcsec: None,
            hfr_arcsec: None,
        };
        let result = dummy_result(vec![star]);
        let config = AnnotationConfig::default();

        let anns_no_flip = compute_annotations(&result, 100, 100, false, &config);
        let anns_flipped = compute_annotations(&result, 100, 100, true, &config);

        assert!((anns_no_flip[0].theta - theta).abs() < 1e-6,
            "without flip, theta should be unchanged");
        assert!((anns_flipped[0].theta - (-theta)).abs() < 1e-6,
            "with flip, theta should be negated: got {} expected {}",
            anns_flipped[0].theta, -theta);
    }

    #[test]
    fn test_bresenham_diagonal() {
        let mut buf = vec![0u8; 10 * 10 * 3];
        draw_line(&mut buf, 10, 10, 3, 0, 0, 9, 9, [255, 0, 0], 1);
        // Check that the diagonal has some red pixels
        let red_count = buf.chunks_exact(3).filter(|px| px[0] == 255).count();
        assert!(red_count >= 10, "Expected at least 10 red pixels on diagonal");
    }

    // ── StarOverlay ──

    /// One annotation against the numbers worked out by hand: positions and semi-axes to 1e-4, the colour exactly.
    fn assert_annotation(
        got: &StarAnnotation,
        (x, y): (f32, f32),
        (semi_major, semi_minor): (f32, f32),
        theta: f32,
        (eccentricity, fwhm): (f32, f32),
        color: [u8; 3],
    ) {
        let near = |a: f32, b: f32| (a - b).abs() < 1e-4;
        assert!(near(got.x, x), "x: {} vs {x}", got.x);
        assert!(near(got.y, y), "y: {} vs {y}", got.y);
        assert!(near(got.semi_major, semi_major), "semi_major: {} vs {semi_major}", got.semi_major);
        assert!(near(got.semi_minor, semi_minor), "semi_minor: {} vs {semi_minor}", got.semi_minor);
        assert!(near(got.theta, theta), "theta: {} vs {theta}", got.theta);
        assert_eq!(got.eccentricity, eccentricity);
        assert_eq!(got.fwhm, fwhm);
        assert_eq!(got.color, color);
    }

    #[test]
    fn legacy_functions_draw_what_they_drew() {
        // A 200 x 100 source drawn at 100 x 50: both scales are 0.5. Default config: ellipse_scale 1.2, radii
        // clamped to [6, 60], ecc_good 0.5, ecc_warn 0.6.
        let round = make_star(40.0, 20.0, 8.0, 0.0); // fwhm 8 -> 8 * 0.5 * 1.2 = 4.8, up to the 6 floor
        let mut tilted = make_star(120.0, 70.0, 40.0, 0.55); // between ecc_good and ecc_warn
        tilted.theta = 0.6;
        let bloated = make_star(180.0, 90.0, 140.0, 0.8); // past ecc_warn; 84 clamps to the 60 ceiling
        let (tilted_minor, bloated_minor) = (tilted.fwhm_y, bloated.fwhm_y);
        let mut result = dummy_result(vec![round, tilted, bloated]);
        result.width = 200;
        result.height = 100;
        let config = AnnotationConfig::default();

        // No flip: y_out = y * 0.5, theta as measured.
        let plain = compute_annotations(&result, 100, 50, false, &config);
        assert_eq!(plain.len(), 3);
        assert_annotation(&plain[0], (20.0, 10.0), (6.0, 6.0), 0.0, (0.0, 8.0), [0, 255, 0]);
        assert_annotation(
            &plain[1],
            (60.0, 35.0),
            (40.0 * 0.5 * 1.2, tilted_minor * 0.5 * 1.2),
            0.6,
            (0.55, 40.0),
            [255, 255, 0],
        );
        assert_annotation(
            &plain[2],
            (90.0, 45.0),
            (60.0, bloated_minor * 0.5 * 1.2),
            0.0,
            (0.8, 140.0),
            [255, 64, 64],
        );

        // Flipped: y_out = h - 1 - y * 0.5 = 49 - y * 0.5, theta negated.
        let flipped = compute_annotations(&result, 100, 50, true, &config);
        assert_eq!(flipped.len(), 3);
        assert_annotation(&flipped[0], (20.0, 39.0), (6.0, 6.0), 0.0, (0.0, 8.0), [0, 255, 0]);
        assert_annotation(
            &flipped[1],
            (60.0, 14.0),
            (40.0 * 0.5 * 1.2, tilted_minor * 0.5 * 1.2),
            -0.6,
            (0.55, 40.0),
            [255, 255, 0],
        );
        assert_annotation(
            &flipped[2],
            (90.0, 4.0),
            (60.0, bloated_minor * 0.5 * 1.2),
            0.0,
            (0.8, 140.0),
            [255, 64, 64],
        );

        // The legacy functions and the overlay functions are the same drawing, flip on and off.
        let overlay = StarOverlay::from(&result);
        assert_eq!(overlay.stars.len(), 3);
        assert_eq!((overlay.width, overlay.height), (200, 100));
        assert_eq!(overlay.median_fwhm, 5.0);
        assert_eq!(plain, compute_overlay(&overlay, 100, 50, false, &config));
        assert_eq!(flipped, compute_overlay(&overlay, 100, 50, true, &config));
        for flip in [false, true] {
            assert_eq!(
                create_annotation_layer(&result, 100, 50, flip, &config),
                create_overlay_layer(&overlay, 100, 50, flip, &config),
            );
            let blank = || ProcessedImage {
                data: vec![0u8; 100 * 50 * 3],
                width: 100,
                height: 50,
                is_color: false,
                channels: 3,
                flip_vertical: flip,
            };
            let (mut legacy, mut shared) = (blank(), blank());
            annotate_image(&mut legacy, &result, &config);
            annotate_image_with(&mut shared, &overlay, &config);
            assert_eq!(legacy.data, shared.data);
            assert!(legacy.data.iter().any(|&b| b > 0), "nothing was drawn");
        }
    }

    #[test]
    fn overlay_flips_like_the_legacy_path() {
        let overlay = StarOverlay {
            stars: vec![OverlayStar {
                x: 10.0,
                y: 5.0,
                fwhm_x: 8.0,
                fwhm_y: 4.0,
                theta: 0.4,
                eccentricity: 0.87,
                fwhm: 5.66,
            }],
            width: 100,
            height: 50,
            median_fwhm: 5.0,
        };
        let config = AnnotationConfig::default();

        // A bottom-up file: the render flips, so row 5 lands on row 50 - 1 - 5 and the angle changes sign.
        let flipped = compute_overlay(&overlay, 100, 50, true, &config);
        assert_eq!(flipped.len(), 1);
        assert!((flipped[0].x - 10.0).abs() < 1e-4);
        assert!((flipped[0].y - 44.0).abs() < 1e-4, "y = {}", flipped[0].y);
        assert!((flipped[0].theta + 0.4).abs() < 1e-6, "theta = {}", flipped[0].theta);

        let plain = compute_overlay(&overlay, 100, 50, false, &config);
        assert!((plain[0].y - 5.0).abs() < 1e-4, "y = {}", plain[0].y);
        assert!((plain[0].theta - 0.4).abs() < 1e-6, "theta = {}", plain[0].theta);
    }

    #[test]
    fn compute_overlay_keeps_the_empty_guard() {
        let star = OverlayStar {
            x: 10.0,
            y: 5.0,
            fwhm_x: 5.0,
            fwhm_y: 5.0,
            theta: 0.0,
            eccentricity: 0.0,
            fwhm: 5.0,
        };
        let config = AnnotationConfig::default();
        let no_stars = StarOverlay { stars: vec![], width: 100, height: 50, median_fwhm: 5.0 };
        let no_width = StarOverlay { stars: vec![star], width: 0, height: 50, median_fwhm: 5.0 };
        let no_height = StarOverlay { stars: vec![star], width: 100, height: 0, median_fwhm: 5.0 };
        for overlay in [&no_stars, &no_width, &no_height] {
            assert!(compute_overlay(overlay, 100, 50, false, &config).is_empty());
            let layer = create_overlay_layer(overlay, 100, 50, false, &config);
            assert_eq!(layer.len(), 100 * 50 * 4);
            assert!(layer.iter().all(|&b| b == 0));
        }
    }

    /// `analysis_mode`'s test helper, copied (it is private to that module's tests).
    fn mono(p: Vec<f32>, w: usize, h: usize) -> RawFrame {
        RawFrame {
            planes: RawPlanes::Mono(p),
            w,
            h,
            domain: Domain::RawAdu {
                clip: 65535.0,
                full_scale: 65535.0,
                gain_e_per_adu: None,
            },
            flip_vertical: false,
            channels: 1,
        }
    }

    /// `n` Moffat (beta 4) stars on a jittered grid of a `w` x `h` plane, each 25 % wider along an axis tilted by
    /// 0.5 rad, on a sky of 1000 with noise sigma 5.
    fn tilted_star_plane(n: usize, w: usize, h: usize) -> Vec<f32> {
        let mut rng = crate::geometry::ransac::SplitMix64(23);
        let cols = ((n as f64 * w as f64 / h as f64).sqrt().ceil() as usize).max(1);
        let rows = n.div_ceil(cols);
        let (cw, ch) = (w as f64 / cols as f64, h as f64 / rows as f64);
        // Moffat beta 4: FWHM = 2 alpha sqrt(2^(1/4) - 1).
        let alpha = 3.5 / (2.0 * (2f64.powf(0.25) - 1.0).sqrt());
        let stars: Vec<MoffatStar> = (0..n)
            .map(|i| {
                let (c, r) = ((i % cols) as f64, (i / cols) as f64);
                MoffatStar {
                    x: c * cw + 12.0 + rng.next_f64() * (cw - 24.0),
                    y: r * ch + 12.0 + rng.next_f64() * (ch - 24.0),
                    amp: 500.0 + rng.next_f64() * 1500.0,
                    alpha_x: alpha * 1.25,
                    alpha_y: alpha,
                    theta: 0.5,
                }
            })
            .collect();
        let mut p = moffat_field(w, h, &stars, 4.0, 1000.0);
        add_noise(&mut p, 5.0, 24);
        p
    }

    #[test]
    fn from_measurement_keeps_every_fit() {
        let (w, h) = (640, 480);
        let frame = mono(tilted_star_plane(120, w, h), w, h);
        let opts = FrameOptions::default();
        let pm = measure_raw_frame(&frame, &opts, None);
        let summary = FrameSummary::from_measurement(&pm, &opts);
        assert!(pm.stars.len() >= 20, "the fixture should measure many stars, got {}", pm.stars.len());

        let overlay = StarOverlay::from_measurement(frame.w, frame.h, &pm, &summary);
        assert_eq!((overlay.width, overlay.height), (w, h));
        assert_eq!(overlay.median_fwhm, summary.median_fwhm as f32);
        assert!(overlay.median_fwhm > 0.0);
        assert_eq!(overlay.stars.len(), pm.stars.len(), "every accepted fit is drawn");
        for (star, fit) in overlay.stars.iter().zip(&pm.stars) {
            assert_eq!(star.x, fit.x as f32);
            assert_eq!(star.y, fit.y as f32);
            assert_eq!(star.fwhm_x, fit.fwhm_x as f32);
            assert_eq!(star.fwhm_y, fit.fwhm_y as f32);
            assert_eq!(star.theta, fit.theta as f32);
            assert_eq!(star.fwhm, fit.fwhm() as f32);
            assert_eq!(star.eccentricity, fit.eccentricity() as f32);
            assert!(star.fwhm_x >= star.fwhm_y, "fwhm_x is the major axis: {star:?}");
        }
        // The fixture is tilted, so the angle carried over is not the all-zero default.
        assert!(overlay.stars.iter().any(|s| s.theta.abs() > 0.1));
    }

    #[test]
    fn median_zero_colours_green() {
        // A frame whose only fits are saturated has no median FWHM (0.0): every star is green under the FWHM scheme,
        // however bloated, instead of dividing by zero.
        let overlay = StarOverlay {
            stars: vec![OverlayStar {
                x: 50.0,
                y: 25.0,
                fwhm_x: 30.0,
                fwhm_y: 30.0,
                theta: 0.0,
                eccentricity: 0.0,
                fwhm: 30.0,
            }],
            width: 100,
            height: 50,
            median_fwhm: 0.0,
        };
        let config = AnnotationConfig {
            color_scheme: ColorScheme::Fwhm,
            ..AnnotationConfig::default()
        };
        let anns = compute_overlay(&overlay, 100, 50, false, &config);
        assert_eq!(anns.len(), 1);
        assert_eq!(anns[0].color, [0, 255, 0]);
    }
}
