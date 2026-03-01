/// Image analysis: FWHM, eccentricity, SNR, PSF signal.

mod background;
mod detection;
mod fitting;
mod metrics;
mod snr;

use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};

use crate::formats;
use crate::processing::color::u16_to_f32;
use crate::processing::debayer;
use crate::processing::stretch::find_median;
use crate::types::{BayerPattern, PixelData};

use detection::DetectionParams;

/// Quantitative metrics for a single detected star.
pub struct StarMetrics {
    /// Subpixel centroid X.
    pub x: f32,
    /// Subpixel centroid Y.
    pub y: f32,
    /// Background-subtracted peak value (ADU).
    pub peak: f32,
    /// Total background-subtracted flux (ADU).
    pub flux: f32,
    /// FWHM along major axis (pixels).
    pub fwhm_x: f32,
    /// FWHM along minor axis (pixels).
    pub fwhm_y: f32,
    /// Geometric mean FWHM (pixels).
    pub fwhm: f32,
    /// Eccentricity: 0 = round, approaching 1 = elongated.
    pub eccentricity: f32,
    /// Per-star aperture photometry SNR.
    pub snr: f32,
    /// Half-flux radius (pixels).
    pub hfr: f32,
    /// PSF position angle in radians, counter-clockwise from +X axis.
    /// Orientation of the major axis (fwhm_x direction).
    /// 0.0 when Gaussian fit is disabled and star is nearly round.
    pub theta: f32,
}

/// Full analysis result for an image.
pub struct AnalysisResult {
    /// Image width (after debayer if applicable).
    pub width: usize,
    /// Image height (after debayer if applicable).
    pub height: usize,
    /// Number of source channels: 1 = mono, 3 = color (after debayer).
    pub source_channels: usize,
    /// Global background level (ADU).
    pub background: f32,
    /// Background noise sigma (ADU).
    pub noise: f32,
    /// Actual detection threshold used (ADU above background).
    pub detection_threshold: f32,
    /// Total stars detected before max_stars cap.
    pub stars_detected: usize,
    /// Per-star metrics, sorted by flux descending, capped at max_stars.
    pub stars: Vec<StarMetrics>,
    /// Median FWHM across all measured stars (pixels).
    pub median_fwhm: f32,
    /// Median eccentricity across all measured stars.
    pub median_eccentricity: f32,
    /// Median per-star SNR.
    pub median_snr: f32,
    /// Median half-flux radius (pixels).
    pub median_hfr: f32,
    /// Image-wide SNR in decibels: 20 × log10(mean_signal / noise). Comparable to PixInsight SNRViews.
    pub snr_db: f32,
    /// PixInsight-style SNR weight: (MeanDev / noise)².
    pub snr_weight: f32,
    /// PSF signal: median(star_peaks) / noise.
    pub psf_signal: f32,
    /// Rayleigh R² statistic for directional coherence of star position angles.
    /// 0.0 = uniform (no trail), 1.0 = all stars aligned (strong trail).
    /// Computed from detection-stage stamp moments on 2θ.
    pub trail_r_squared: f32,
    /// True if the image is likely trailed, based on the Rayleigh test.
    /// Uses configurable R² threshold (default 0.5) and eccentricity gate.
    pub possibly_trailed: bool,
}

/// Builder configuration for analysis (internal).
pub struct AnalysisConfig {
    detection_sigma: f32,
    min_star_area: usize,
    max_star_area: usize,
    saturation_fraction: f32,
    max_stars: usize,
    use_gaussian_fit: bool,
    background_mesh_size: Option<usize>,
    apply_debayer: bool,
    max_eccentricity: f32,
    trail_r_squared_threshold: f32,
}

/// Image analyzer with builder pattern.
pub struct ImageAnalyzer {
    config: AnalysisConfig,
    thread_pool: Option<Arc<rayon::ThreadPool>>,
}

impl ImageAnalyzer {
    pub fn new() -> Self {
        ImageAnalyzer {
            config: AnalysisConfig {
                detection_sigma: 5.0,
                min_star_area: 5,
                max_star_area: 2000,
                saturation_fraction: 0.95,
                max_stars: 200,
                use_gaussian_fit: true,
                background_mesh_size: None,
                apply_debayer: true,
                max_eccentricity: 1.0,
                trail_r_squared_threshold: 0.5,
            },
            thread_pool: None,
        }
    }

    /// Star detection threshold in σ above background.
    pub fn with_detection_sigma(mut self, sigma: f32) -> Self {
        self.config.detection_sigma = sigma.max(1.0);
        self
    }

    /// Reject connected components with fewer pixels than this (filters hot pixels).
    pub fn with_min_star_area(mut self, area: usize) -> Self {
        self.config.min_star_area = area.max(1);
        self
    }

    /// Reject connected components with more pixels than this (filters galaxies/nebulae).
    pub fn with_max_star_area(mut self, area: usize) -> Self {
        self.config.max_star_area = area;
        self
    }

    /// Reject stars with peak > fraction × 65535 (saturated).
    pub fn with_saturation_fraction(mut self, frac: f32) -> Self {
        self.config.saturation_fraction = frac.clamp(0.5, 1.0);
        self
    }

    /// Analyze only the brightest N stars.
    pub fn with_max_stars(mut self, n: usize) -> Self {
        self.config.max_stars = n.max(1);
        self
    }

    /// Use fast windowed-moments instead of Gaussian fit for FWHM (less accurate but faster).
    pub fn without_gaussian_fit(mut self) -> Self {
        self.config.use_gaussian_fit = false;
        self
    }

    /// Enable mesh-grid background estimation with given cell size (handles gradients).
    pub fn with_background_mesh(mut self, cell_size: usize) -> Self {
        self.config.background_mesh_size = Some(cell_size.max(16));
        self
    }

    /// Skip debayering for OSC images (less accurate but faster).
    pub fn without_debayer(mut self) -> Self {
        self.config.apply_debayer = false;
        self
    }

    /// Reject stars with eccentricity above this threshold (filters elongated objects/trails).
    /// Default: 0.5. Set to 1.0 to disable.
    pub fn with_max_eccentricity(mut self, ecc: f32) -> Self {
        self.config.max_eccentricity = ecc.clamp(0.0, 1.0);
        self
    }

    /// Set the R² threshold for trail detection.
    /// Images with Rayleigh R² above this are flagged as possibly trailed.
    /// Default: 0.5. Lower values are more aggressive (more false positives).
    pub fn with_trail_threshold(mut self, threshold: f32) -> Self {
        self.config.trail_r_squared_threshold = threshold.clamp(0.0, 1.0);
        self
    }

    /// Use a custom rayon thread pool.
    pub fn with_thread_pool(mut self, pool: Arc<rayon::ThreadPool>) -> Self {
        self.thread_pool = Some(pool);
        self
    }

    /// Analyze a FITS or XISF image file.
    pub fn analyze<P: AsRef<Path>>(&self, path: P) -> Result<AnalysisResult> {
        let path = path.as_ref();
        match &self.thread_pool {
            Some(pool) => pool.install(|| self.analyze_impl(path)),
            None => self.analyze_impl(path),
        }
    }

    /// Analyze pre-loaded f32 pixel data.
    ///
    /// `data`: planar f32 pixel data (for 3-channel: RRRGGGBBB layout).
    /// `width`: image width.
    /// `height`: image height.
    /// `channels`: 1 for mono, 3 for RGB.
    pub fn analyze_data(
        &self,
        data: &[f32],
        width: usize,
        height: usize,
        channels: usize,
    ) -> Result<AnalysisResult> {
        match &self.thread_pool {
            Some(pool) => pool.install(|| {
                self.run_analysis(data, width, height, channels, None)
            }),
            None => self.run_analysis(data, width, height, channels, None),
        }
    }

    fn analyze_impl(&self, path: &Path) -> Result<AnalysisResult> {
        let (meta, pixel_data) =
            formats::read_image(path).context("Failed to read image for analysis")?;

        // Convert to f32
        let f32_data = match pixel_data {
            PixelData::Float32(d) => d,
            PixelData::Uint16(d) => u16_to_f32(&d),
        };

        let mut data = f32_data;
        let width = meta.width;
        let height = meta.height;
        let channels = meta.channels;

        // Green-channel interpolation for OSC: native-resolution mono green image
        // (matches Siril's interpolate_nongreen — no PSF broadening from 2×2 binning)
        let green_mask = if self.config.apply_debayer
            && meta.bayer_pattern != BayerPattern::None
            && channels == 1
        {
            data = debayer::interpolate_green_f32(&data, width, height, meta.bayer_pattern);
            // width, height, channels unchanged — native resolution mono green
            Some(build_green_mask(width, height, meta.bayer_pattern))
        } else {
            None
        };

        self.run_analysis(&data, width, height, channels, green_mask.as_deref())
    }

    fn run_analysis(
        &self,
        data: &[f32],
        width: usize,
        height: usize,
        channels: usize,
        green_mask: Option<&[bool]>,
    ) -> Result<AnalysisResult> {
        // Extract luminance if multi-channel
        let lum = if channels == 3 {
            extract_luminance(data, width, height)
        } else {
            data[..width * height].to_vec()
        };

        // Background estimation
        let bg_result = if let Some(cell_size) = self.config.background_mesh_size {
            background::estimate_background_mesh(&lum, width, height, cell_size)
        } else {
            background::estimate_background(&lum, width, height)
        };

        let bg_map_ref = bg_result.background_map.as_deref();

        // Star detection
        let det_params = DetectionParams {
            detection_sigma: self.config.detection_sigma,
            min_star_area: self.config.min_star_area,
            max_star_area: self.config.max_star_area,
            saturation_limit: self.config.saturation_fraction * 65535.0,
            max_stars: self.config.max_stars,
        };

        let detected = detection::detect_stars(
            &lum,
            width,
            height,
            bg_result.background,
            bg_result.noise,
            bg_map_ref,
            &det_params,
        );

        let stars_detected = detected.len();
        let detection_threshold = self.config.detection_sigma * bg_result.noise;

        // ── Phase 1: Image-level trailing detection (Rayleigh test) ─────
        // Compute R² and advisory flag; never reject — let callers decide.
        //
        // Two paths flag trailing via Rayleigh test on 2θ:
        //
        // Path A — Strong R² (> threshold): Flag regardless of eccentricity.
        //   Real trails produce R² > 0.7 (strong directional coherence).
        //   Grid-induced coherence gives R² ≈ 0.15 (undersampled) to 0.40
        //   (oversampled with field-angle effects like coma/curvature).
        //
        // Path B — Eccentricity-gated (median ecc > 0.6): standard Rayleigh.
        //   For undersampled stars, moments are noisy and theta has grid bias.
        //   Only flag when blobs are genuinely elongated (ecc > 0.6).
        let (trail_r_squared, possibly_trailed) = if detected.len() >= 5 {
            let n = detected.len();
            let (sum_cos, sum_sin) =
                detected
                    .iter()
                    .fold((0.0f64, 0.0f64), |(sc, ss), s| {
                        let a = 2.0 * s.theta as f64;
                        (sc + a.cos(), ss + a.sin())
                    });
            let r_sq = (sum_cos * sum_cos + sum_sin * sum_sin) / (n as f64 * n as f64);
            let p = (-(n as f64) * r_sq).exp();

            let mut eccs: Vec<f32> = detected.iter().map(|s| s.eccentricity).collect();
            eccs.sort_unstable_by(|a, b| a.total_cmp(b));
            let median_ecc = eccs[eccs.len() / 2];

            let threshold = self.config.trail_r_squared_threshold as f64;
            let trailed = (r_sq > threshold && p < 0.01)
                || (median_ecc > 0.6 && p < 0.05);

            (r_sq as f32, trailed)
        } else {
            (0.0, false)
        };

        let zero_result = || {
            Ok(AnalysisResult {
                width,
                height,
                source_channels: channels,
                background: bg_result.background,
                noise: bg_result.noise,
                detection_threshold,
                stars_detected,
                stars: Vec::new(),
                median_fwhm: 0.0,
                median_eccentricity: 0.0,
                median_snr: 0.0,
                median_hfr: 0.0,
                snr_db: snr::compute_snr_db(&lum, bg_result.noise),
                snr_weight: snr::compute_snr_weight(&lum, bg_result.background, bg_result.noise),
                psf_signal: 0.0,
                trail_r_squared,
                possibly_trailed,
            })
        };

        if stars_detected == 0 {
            return zero_result();
        }

        // Measure PSF metrics
        let mut measured = metrics::measure_stars(
            &lum,
            width,
            height,
            &detected,
            bg_result.background,
            bg_map_ref,
            self.config.use_gaussian_fit,
            green_mask,
        );

        if measured.is_empty() {
            return zero_result();
        }

        // ── Phase 2: Per-star eccentricity filter (after measurement) ───
        measured.retain(|s| s.eccentricity <= self.config.max_eccentricity);

        if measured.is_empty() {
            return zero_result();
        }

        // Compute median FWHM for aperture sizing
        let mut fwhm_vals: Vec<f32> = measured.iter().map(|s| s.fwhm).collect();
        let median_fwhm = find_median(&mut fwhm_vals);

        // Per-star SNR
        snr::compute_star_snr(&lum, width, height, &mut measured, median_fwhm);

        // Summary statistics
        let mut ecc_vals: Vec<f32> = measured.iter().map(|s| s.eccentricity).collect();
        let mut snr_vals: Vec<f32> = measured.iter().map(|s| s.snr).collect();
        let mut hfr_vals: Vec<f32> = measured.iter().map(|s| s.hfr).collect();

        let median_eccentricity = find_median(&mut ecc_vals);
        let median_snr = find_median(&mut snr_vals);
        let median_hfr = find_median(&mut hfr_vals);

        let snr_db = snr::compute_snr_db(&lum, bg_result.noise);
        let snr_weight = snr::compute_snr_weight(&lum, bg_result.background, bg_result.noise);
        let psf_signal = snr::compute_psf_signal(&measured, bg_result.noise);

        // Convert to public StarMetrics
        let stars: Vec<StarMetrics> = measured
            .into_iter()
            .map(|m| StarMetrics {
                x: m.x,
                y: m.y,
                peak: m.peak,
                flux: m.flux,
                fwhm_x: m.fwhm_x,
                fwhm_y: m.fwhm_y,
                fwhm: m.fwhm,
                eccentricity: m.eccentricity,
                snr: m.snr,
                hfr: m.hfr,
                theta: m.theta,
            })
            .collect();

        Ok(AnalysisResult {
            width,
            height,
            source_channels: channels,
            background: bg_result.background,
            noise: bg_result.noise,
            detection_threshold,
            stars_detected,
            stars,
            median_fwhm,
            median_eccentricity,
            median_snr,
            median_hfr,
            snr_db,
            snr_weight,
            psf_signal,
            trail_r_squared,
            possibly_trailed,
        })
    }
}

impl Default for ImageAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

/// Build a boolean mask marking green CFA pixel positions.
///
/// Returns a `Vec<bool>` of length `width * height` where `true` marks pixels
/// that are at green positions in the Bayer pattern. For GBRG/GRBG green is
/// at (row + col) even; for RGGB/BGGR green is at (row + col) odd.
fn build_green_mask(width: usize, height: usize, pattern: BayerPattern) -> Vec<bool> {
    let green_at_even = matches!(pattern, BayerPattern::Gbrg | BayerPattern::Grbg);
    (0..height)
        .flat_map(|y| {
            (0..width).map(move |x| {
                let parity = (x + y) & 1;
                if green_at_even { parity == 0 } else { parity == 1 }
            })
        })
        .collect()
}

/// Extract luminance from planar RGB data: L = 0.2126R + 0.7152G + 0.0722B
fn extract_luminance(data: &[f32], width: usize, height: usize) -> Vec<f32> {
    use rayon::prelude::*;

    let plane_size = width * height;
    let r = &data[..plane_size];
    let g = &data[plane_size..2 * plane_size];
    let b = &data[2 * plane_size..3 * plane_size];

    let mut lum = vec![0.0_f32; plane_size];
    const CHUNK: usize = 8192;
    lum.par_chunks_mut(CHUNK)
        .enumerate()
        .for_each(|(ci, chunk)| {
            let off = ci * CHUNK;
            for (i, dst) in chunk.iter_mut().enumerate() {
                let idx = off + i;
                *dst = 0.2126 * r[idx] + 0.7152 * g[idx] + 0.0722 * b[idx];
            }
        });
    lum
}
