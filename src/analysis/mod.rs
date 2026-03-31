/// Image analysis: FWHM, eccentricity, SNR, PSF signal.

#[cfg(feature = "debug-pipeline")]
pub mod background;
#[cfg(not(feature = "debug-pipeline"))]
mod background;

#[cfg(feature = "debug-pipeline")]
pub mod convolution;
#[cfg(not(feature = "debug-pipeline"))]
mod convolution;

#[cfg(feature = "debug-pipeline")]
pub mod detection;
#[cfg(not(feature = "debug-pipeline"))]
mod detection;

#[cfg(feature = "debug-pipeline")]
pub mod fitting;
#[cfg(not(feature = "debug-pipeline"))]
mod fitting;

#[cfg(feature = "debug-pipeline")]
pub mod metrics;
#[cfg(not(feature = "debug-pipeline"))]
mod metrics;

#[cfg(feature = "debug-pipeline")]
pub mod snr;
#[cfg(not(feature = "debug-pipeline"))]
mod snr;

#[cfg(feature = "debug-pipeline")]
pub mod render;

use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};

use crate::formats;
use crate::processing::color::u16_to_f32;
use crate::processing::debayer;
use crate::processing::stretch::find_median;
use crate::types::{BayerPattern, ImageMetadata, PixelData};

use detection::DetectionParams;

/// Method used to measure this star's PSF.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FitMethod {
    /// Free-beta Moffat (8 params) — highest accuracy.
    FreeMoffat,
    /// Fixed-beta Moffat (7 params) — field median beta.
    FixedMoffat,
    /// Gaussian fallback (7 params).
    Gaussian,
    /// Windowed moments — lowest accuracy, flagged unreliable.
    Moments,
}

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
    /// Moffat β parameter (None if Gaussian/moments fit was used).
    pub beta: Option<f32>,
    /// Which PSF fitting method produced this measurement.
    pub fit_method: FitMethod,
    /// Normalized fit residual (quality weight: w = 1/(1+r)).
    /// Lower = better fit. 1.0 for moments fallback.
    pub fit_residual: f32,
    /// FWHM in arcseconds (None if optics not provided via with_optics).
    pub fwhm_arcsec: Option<f32>,
    /// HFR in arcseconds (None if optics not provided via with_optics).
    pub hfr_arcsec: Option<f32>,
}

/// Per-stage timing in milliseconds for the analysis pipeline.
pub struct StageTiming {
    pub background_ms: f64,
    pub detection_pass1_ms: f64,
    pub calibration_ms: f64,
    pub detection_pass2_ms: f64,
    pub measurement_ms: f64,
    pub snr_ms: f64,
    pub statistics_ms: f64,
    pub total_ms: f64,
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
    /// Total stars detected (raw detection count, before measure cap).
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
    /// SNR weight for frame ranking: median(star_flux)² / (noise² × background).
    /// Star-based metric immune to background gradients. Higher = better frame.
    pub snr_weight: f32,
    /// PSF signal: median(star_peaks) / noise.
    pub psf_signal: f32,
    /// Per-frame SNR: background / noise (linear ratio).
    /// Use for stacking prediction: stacked_snr = sqrt(sum(frame_snr_i²)).
    pub frame_snr: f32,
    /// Rayleigh R̄² (squared mean resultant length) for directional coherence
    /// of star position angles. Uses 2θ doubling for axial orientation data.
    /// 0.0 = uniform (no trail), 1.0 = all stars aligned (strong trail).
    /// A threshold of 0.5 corresponds to R̄ ≈ 0.71 (strong coherence).
    pub trail_r_squared: f32,
    /// True if the image is likely trailed, based on the Rayleigh test on
    /// PSF-fit stars. Fires when R̄² > threshold with significant p-value,
    /// or when R̄² > 0.15 with high median eccentricity. Suppressed if the
    /// elongation pattern is optical (radial coma or ecc-distance tilt gradient).
    /// Requires ≥20 measured stars and FWHM ≥ 2.0 px.
    pub possibly_trailed: bool,
    /// Measured FWHM from adaptive two-pass detection (pixels).
    /// This is the FWHM used for the final matched filter kernel.
    /// If the first-pass FWHM was within 30% of 3.0, this equals 3.0.
    pub measured_fwhm_kernel: f32,
    /// Median Moffat β across all stars (None if Moffat fitting not used).
    /// Typical range: 2.0-5.0 for real optics. Lower = broader wings.
    pub median_beta: Option<f32>,
    /// Pass 1 detection count (before calibration/re-detection).
    pub pass1_detections: usize,
    /// Calibrated field FWHM from Moffat pass (pixels, before capping).
    pub calibrated_fwhm: f32,
    /// Number of stars that survived PSF fitting (before late truncation).
    pub stars_measured: usize,
    /// Number of Moffat fits (FreeMoffat + FixedMoffat) among measured stars.
    pub moffat_count: usize,
    /// Number of Gaussian fits among measured stars.
    pub gaussian_count: usize,
    /// Plate scale in arcsec/pixel (None if optics not provided).
    pub plate_scale: Option<f32>,
    /// Median FWHM in arcseconds (None if optics not provided).
    pub median_fwhm_arcsec: Option<f32>,
    /// Median HFR in arcseconds (None if optics not provided).
    pub median_hfr_arcsec: Option<f32>,
    /// Per-stage timing breakdown for the analysis pipeline.
    pub stage_timing: StageTiming,
}

/// Builder configuration for analysis (internal).
pub struct AnalysisConfig {
    detection_sigma: f32,
    min_star_area: usize,
    max_star_area: usize,
    saturation_fraction: f32,
    max_stars: usize,
    apply_debayer: bool,
    trail_r_squared_threshold: f32,
    /// MRS wavelet noise layers (default 4).
    noise_layers: usize,
    /// Max stars to PSF-fit for statistics. 0 = measure all.
    measure_cap: usize,
    /// LM max iterations for pass-2 measurement fits.
    fit_max_iter: usize,
    /// LM convergence tolerance for pass-2 measurement fits.
    fit_tolerance: f64,
    /// Consecutive LM step rejects before early bailout.
    fit_max_rejects: usize,
    /// Telescope focal length in millimeters (for arcsecond measurements).
    focal_length_mm: Option<f64>,
    /// Camera pixel size in micrometers (for arcsecond measurements).
    pixel_size_um: Option<f64>,
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
                apply_debayer: true,
                trail_r_squared_threshold: 0.5,
                noise_layers: 0,
                measure_cap: 500,
                fit_max_iter: 25,
                fit_tolerance: 1e-4,
                fit_max_rejects: 5,
                focal_length_mm: None,
                pixel_size_um: None,
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

    /// Keep only the brightest N stars in the returned result.
    pub fn with_max_stars(mut self, n: usize) -> Self {
        self.config.max_stars = n.max(1);
        self
    }

    /// Skip debayering for OSC images (less accurate but faster).
    pub fn without_debayer(mut self) -> Self {
        self.config.apply_debayer = false;
        self
    }

    /// Set the R² threshold for trail detection.
    /// Images with Rayleigh R² above this are flagged as possibly trailed.
    /// Default: 0.5. Lower values are more aggressive (more false positives).
    pub fn with_trail_threshold(mut self, threshold: f32) -> Self {
        self.config.trail_r_squared_threshold = threshold.clamp(0.0, 1.0);
        self
    }

    /// Set optics parameters for arcsecond-based measurements.
    /// `focal_length_mm`: telescope focal length in millimeters.
    /// `pixel_size_um`: camera pixel size in micrometers.
    /// When both are set, FWHM and HFR are reported in arcseconds alongside pixels.
    pub fn with_optics(mut self, focal_length_mm: f64, pixel_size_um: f64) -> Self {
        self.config.focal_length_mm = Some(focal_length_mm);
        self.config.pixel_size_um = Some(pixel_size_um);
        self
    }

    /// Set MRS wavelet noise layers for noise estimation.
    /// Default: 0 (fast MAD noise from mesh-grid cell sigmas).
    /// Set to 1-6 for MRS wavelet noise (more robust against nebulosity/gradients,
    /// ~200ms slower per frame). 4 is the recommended MRS setting.
    pub fn with_mrs_layers(mut self, layers: usize) -> Self {
        self.config.noise_layers = layers;
        self
    }

    /// Max stars to PSF-fit for statistics. Default 2000.
    /// Stars are sorted by flux (brightest first) before capping.
    /// Set to 0 to measure all detected stars (catalog export mode).
    pub fn with_measure_cap(mut self, n: usize) -> Self {
        self.config.measure_cap = n;
        self
    }

    /// LM max iterations for pass-2 measurement fits. Default 25.
    /// Calibration pass always uses 50 iterations.
    pub fn with_fit_max_iter(mut self, n: usize) -> Self {
        self.config.fit_max_iter = n.max(1);
        self
    }

    /// LM convergence tolerance for pass-2 measurement fits. Default 1e-4.
    /// Calibration pass always uses 1e-6.
    pub fn with_fit_tolerance(mut self, tol: f64) -> Self {
        self.config.fit_tolerance = tol;
        self
    }

    /// Consecutive LM step rejects before early bailout. Default 5.
    pub fn with_fit_max_rejects(mut self, n: usize) -> Self {
        self.config.fit_max_rejects = n.max(1);
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
                self.run_analysis(data, width, height, channels)
            }),
            None => self.run_analysis(data, width, height, channels),
        }
    }

    /// Analyze pre-read raw pixel data (skips file I/O).
    ///
    /// Accepts `ImageMetadata` and borrows `PixelData`, handling u16→f32
    /// conversion and green-channel interpolation for OSC images internally.
    pub fn analyze_raw(
        &self,
        meta: &ImageMetadata,
        pixels: &PixelData,
    ) -> Result<AnalysisResult> {
        match &self.thread_pool {
            Some(pool) => pool.install(|| self.analyze_raw_impl(meta, pixels)),
            None => self.analyze_raw_impl(meta, pixels),
        }
    }

    /// Analyze multiple images in parallel.
    ///
    /// `concurrency` controls how many frames are analyzed simultaneously.
    /// `progress` is called after each frame completes with (completed, total, path).
    /// Returns results in approximate completion order.
    pub fn analyze_batch<P, F>(
        &self,
        paths: &[P],
        concurrency: usize,
        progress: F,
    ) -> Vec<(std::path::PathBuf, Result<AnalysisResult>)>
    where
        P: AsRef<std::path::Path> + Sync,
        F: Fn(usize, usize, &std::path::Path) + Send + Sync,
    {
        use std::sync::atomic::{AtomicUsize, Ordering};
        use rayon::prelude::*;

        let total = paths.len();
        let completed = AtomicUsize::new(0);
        let concurrency = concurrency.max(1);

        let do_batch = || {
            let mut results = Vec::with_capacity(total);
            for chunk in paths.chunks(concurrency) {
                let chunk_results: Vec<_> = chunk
                    .into_par_iter()
                    .map(|p| {
                        let path = p.as_ref();
                        let result = self.analyze_impl(path);
                        let n = completed.fetch_add(1, Ordering::Relaxed) + 1;
                        progress(n, total, path);
                        (path.to_path_buf(), result)
                    })
                    .collect();
                results.extend(chunk_results);
            }
            results
        };

        match &self.thread_pool {
            Some(pool) => pool.install(do_batch),
            None => do_batch(),
        }
    }

    fn analyze_raw_impl(
        &self,
        meta: &ImageMetadata,
        pixels: &PixelData,
    ) -> Result<AnalysisResult> {
        let f32_data = match pixels {
            PixelData::Float32(d) => std::borrow::Cow::Borrowed(d.as_slice()),
            PixelData::Uint16(d) => std::borrow::Cow::Owned(u16_to_f32(d)),
        };

        let mut data = f32_data.into_owned();
        let width = meta.width;
        let height = meta.height;
        let channels = meta.channels;

        // OSC green interpolation: replace R/B pixels with weighted average of
        // neighboring green values.  PSF fitting uses all pixels — no green mask.
        if self.config.apply_debayer
            && meta.bayer_pattern != BayerPattern::None
            && channels == 1
        {
            data = debayer::interpolate_green_f32(&data, width, height, meta.bayer_pattern);
        }

        self.run_analysis(&data, width, height, channels)
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

        // OSC green interpolation for Bayer data.
        if self.config.apply_debayer
            && meta.bayer_pattern != BayerPattern::None
            && channels == 1
        {
            data = debayer::interpolate_green_f32(&data, width, height, meta.bayer_pattern);
        }

        self.run_analysis(&data, width, height, channels)
    }

    fn run_analysis(
        &self,
        data: &[f32],
        width: usize,
        height: usize,
        channels: usize,
    ) -> Result<AnalysisResult> {
        let pipeline_start = std::time::Instant::now();

        // Extract luminance if multi-channel
        let lum = if channels == 3 {
            extract_luminance(data, width, height)
        } else {
            data[..width * height].to_vec()
        };

        let det_params = DetectionParams {
            detection_sigma: self.config.detection_sigma,
            min_star_area: self.config.min_star_area,
            max_star_area: self.config.max_star_area,
            saturation_limit: self.config.saturation_fraction * 65535.0,
        };

        // ── Stage 1: Background & Noise ──────────────────────────────────
        let t = std::time::Instant::now();
        let cell_size = background::auto_cell_size(width, height);
        let mut bg_result = background::estimate_background_mesh(&lum, width, height, cell_size);
        if self.config.noise_layers > 0 {
            // MRS wavelet noise: accurate but ~500ms. Layers 1-6.
            bg_result.noise = background::estimate_noise_mrs(
                &lum, width, height, self.config.noise_layers.max(1),
            ).max(0.001);
        }
        // noise_layers == 0: keep MAD noise from mesh-grid (already in bg_result.noise)
        let background_ms = t.elapsed().as_secs_f64() * 1000.0;

        // ── Stage 2, Pass 1: Discovery ───────────────────────────────────
        let t = std::time::Instant::now();
        let initial_fwhm = 3.0_f32;
        let pass1_stars = {
            let bg_map = bg_result.background_map.as_deref();
            let noise_map = bg_result.noise_map.as_deref();
            detection::detect_stars(
                &lum, width, height,
                bg_result.background, bg_result.noise,
                bg_map, noise_map, &det_params, initial_fwhm,
                None,
            )
        };
        let detection_pass1_ms = t.elapsed().as_secs_f64() * 1000.0;

        // Select calibration stars: brightest, not saturated, not too elongated
        let t = std::time::Instant::now();
        let calibration_stars: Vec<&detection::DetectedStar> = pass1_stars
            .iter()
            .filter(|s| s.eccentricity < 0.5 && s.area >= 5)
            .take(100)
            .collect();

        // Free-beta Moffat on calibration stars to discover field PSF model.
        // calibration_ok tracks whether we got reliable results — if not,
        // moments screening is disabled for the measurement pass.
        let field_beta: Option<f64>;
        let field_fwhm: f32;
        if calibration_stars.len() >= 3 {
            let cal_owned: Vec<detection::DetectedStar> = calibration_stars
                .iter()
                .map(|s| detection::DetectedStar {
                    x: s.x, y: s.y, peak: s.peak, flux: s.flux,
                    area: s.area, theta: s.theta, eccentricity: s.eccentricity,
                })
                .collect();
            let cal_measured = metrics::measure_stars(
                &lum, width, height, &cal_owned,
                bg_result.background,
                bg_result.background_map.as_deref(),
                None, // fit all pixels in green-interpolated image (no green mask)
                None, // free-beta Moffat
                50, 1e-6, 5, // calibration always uses full precision
                None,   // no screening for calibration
            );

            let mut beta_vals: Vec<f32> = cal_measured.iter().filter_map(|s| s.beta).collect();
            let mut fwhm_vals: Vec<f32> = cal_measured.iter().map(|s| s.fwhm).collect();

            if beta_vals.len() >= 3 {
                field_beta = Some(sigma_clipped_median(&beta_vals) as f64);
            } else if !beta_vals.is_empty() {
                field_beta = Some(find_median(&mut beta_vals) as f64);
            } else {
                field_beta = None;
            }

            if fwhm_vals.len() >= 3 {
                field_fwhm = sigma_clipped_median(&fwhm_vals);
            } else if !fwhm_vals.is_empty() {
                field_fwhm = find_median(&mut fwhm_vals);
            } else {
                field_fwhm = estimate_fwhm_from_stars(
                    &lum, width, height, &pass1_stars,
                    bg_result.background, bg_result.background_map.as_deref(),
                );
            }
        } else {
            // Too few calibration stars (e.g., trailed image where most stars
            // have ecc > 0.5).  Fall back to halfmax estimate.
            field_beta = None;
            field_fwhm = estimate_fwhm_from_stars(
                &lum, width, height, &pass1_stars,
                bg_result.background, bg_result.background_map.as_deref(),
            );
        }

        let calibration_ms = t.elapsed().as_secs_f64() * 1000.0;

        // Source-mask background re-estimation skipped for speed.
        // The sigma-clipped cell stats already reject >30% star-contaminated
        // cells and 3-round sigma clipping handles residual star flux within
        // cells. The source mask adds ~400ms for marginal improvement.

        // ── Stage 2, Pass 2: Full detection with refined kernel ──────────
        let t = std::time::Instant::now();
        // Clamp minimum FWHM to 2.0px — no real optics produce sub-2px stars,
        // and tiny kernels have poor noise rejection (OSC green-channel fits
        // can underestimate FWHM due to Bayer grid undersampling).
        let clamped_fwhm = field_fwhm.max(2.0);
        let final_fwhm = if clamped_fwhm > 1.0
            && ((clamped_fwhm - initial_fwhm) / initial_fwhm).abs() > 0.30
        {
            clamped_fwhm.min(initial_fwhm * 6.0)
        } else {
            initial_fwhm
        };

        let detected = {
            let bg_map = bg_result.background_map.as_deref();
            let noise_map = bg_result.noise_map.as_deref();
            detection::detect_stars(
                &lum, width, height,
                bg_result.background, bg_result.noise,
                bg_map, noise_map, &det_params, final_fwhm,
                Some(clamped_fwhm),
            )
        };
        let detection_pass2_ms = t.elapsed().as_secs_f64() * 1000.0;

        let bg_map_ref = bg_result.background_map.as_deref();
        let detection_threshold = self.config.detection_sigma * bg_result.noise;

        // Trail detection defaults — computed after PSF measurement (see below).
        let mut trail_r_squared = 0.0_f32;
        let mut possibly_trailed = false;

        let frame_snr = if bg_result.noise > 0.0 { bg_result.background / bg_result.noise } else { 0.0 };

        let pass1_detections = pass1_stars.len();

        let make_zero_result = |stars_detected: usize| {
            Ok(AnalysisResult {
                width, height, source_channels: channels,
                background: bg_result.background, noise: bg_result.noise,
                detection_threshold, stars_detected,
                stars: Vec::new(),
                median_fwhm: 0.0, median_eccentricity: 0.0,
                median_snr: 0.0, median_hfr: 0.0,
                snr_weight: 0.0, psf_signal: 0.0, frame_snr,
                trail_r_squared, possibly_trailed,
                measured_fwhm_kernel: final_fwhm,
                median_beta: field_beta.map(|b| b as f32),
                pass1_detections: 0, calibrated_fwhm: 0.0,
                stars_measured: 0, moffat_count: 0, gaussian_count: 0,
                plate_scale: None, median_fwhm_arcsec: None, median_hfr_arcsec: None,
                stage_timing: StageTiming {
                    background_ms: 0.0, detection_pass1_ms: 0.0, calibration_ms: 0.0,
                    detection_pass2_ms: 0.0, measurement_ms: 0.0, snr_ms: 0.0,
                    statistics_ms: 0.0, total_ms: pipeline_start.elapsed().as_secs_f64() * 1000.0,
                },
            })
        };

        if detected.is_empty() {
            return make_zero_result(0);
        }

        // ── Stage 3: PSF Measurement (with measure cap) ─────────────────
        let t = std::time::Instant::now();
        let stars_detected = detected.len();

        // Apply measure cap with spatial grid balancing.
        // Divide image into 4×4 grid, round-robin select from each cell
        // to ensure spatial coverage across the field.
        let effective_cap = if self.config.measure_cap == 0 {
            detected.len()
        } else {
            self.config.measure_cap
        };

        let to_measure: Vec<detection::DetectedStar> = if detected.len() <= effective_cap {
            detected.clone()
        } else {
            debug_assert!(
                detected.windows(2).all(|w| w[0].flux >= w[1].flux),
                "detected stars must be sorted by flux descending"
            );
            const GRID_N: usize = 4;
            let cell_w = width as f32 / GRID_N as f32;
            let cell_h = height as f32 / GRID_N as f32;
            let mut buckets: Vec<Vec<&detection::DetectedStar>> =
                vec![Vec::new(); GRID_N * GRID_N];

            for star in &detected {
                let gx = ((star.x / cell_w) as usize).min(GRID_N - 1);
                let gy = ((star.y / cell_h) as usize).min(GRID_N - 1);
                buckets[gy * GRID_N + gx].push(star);
            }

            let mut selected: Vec<detection::DetectedStar> = Vec::with_capacity(effective_cap);
            let mut idx = vec![0usize; GRID_N * GRID_N];
            loop {
                let mut added_any = false;
                for cell in 0..(GRID_N * GRID_N) {
                    if selected.len() >= effective_cap { break; }
                    if idx[cell] < buckets[cell].len() {
                        selected.push(buckets[cell][idx[cell]].clone());
                        idx[cell] += 1;
                        added_any = true;
                    }
                }
                if !added_any || selected.len() >= effective_cap { break; }
            }
            selected
        };

        let mut measured = metrics::measure_stars(
            &lum, width, height, &to_measure,
            bg_result.background, bg_map_ref,
            None, field_beta, // fit all pixels in green-interpolated image
            self.config.fit_max_iter,
            self.config.fit_tolerance,
            self.config.fit_max_rejects,
            if field_fwhm > 1.0 { Some(field_fwhm) } else { None },  // adaptive screening when FWHM is reliable
        );
        let measurement_ms = t.elapsed().as_secs_f64() * 1000.0;
        let stars_measured = measured.len();
        let moffat_count = measured.iter()
            .filter(|s| matches!(s.fit_method, FitMethod::FreeMoffat | FitMethod::FixedMoffat))
            .count();
        let gaussian_count = measured.iter()
            .filter(|s| matches!(s.fit_method, FitMethod::Gaussian))
            .count();

        if measured.is_empty() {
            return make_zero_result(stars_detected);
        }

        // ── Trail detection (Rayleigh test on PSF-fit stars) ─────────────
        // Uses PSF-fit theta and eccentricity (more accurate than detection moments).
        // Requires FWHM >= 2.0 px (pixel grid quantization) and ≥20 measured stars.
        if measured.len() >= 20 && field_fwhm >= 2.0 {
            let n = measured.len();
            let (sum_cos, sum_sin) =
                measured.iter().fold((0.0f64, 0.0f64), |(sc, ss), s| {
                    let a = 2.0 * s.theta as f64;
                    (sc + a.cos(), ss + a.sin())
                });
            let r_sq = (sum_cos * sum_cos + sum_sin * sum_sin) / (n as f64 * n as f64);
            let p = (-(n as f64) * r_sq).exp();
            let mut eccs: Vec<f32> = measured.iter().map(|s| s.eccentricity).collect();
            eccs.sort_unstable_by(|a, b| a.total_cmp(b));
            let median_ecc = if eccs.len() % 2 == 1 {
                eccs[eccs.len() / 2]
            } else {
                (eccs[eccs.len() / 2 - 1] + eccs[eccs.len() / 2]) * 0.5
            };
            let threshold = self.config.trail_r_squared_threshold as f64;
            trail_r_squared = r_sq as f32;
            possibly_trailed = (r_sq > threshold && p < 0.01)
                || (r_sq > 0.15 && median_ecc > 0.7 && p < 0.05);

            // ── Optical aberration suppression ───────────────────────────
            // If Rayleigh fires, check if the elongation pattern is optical
            // (coma/tilt) rather than parallel trailing.
            if possibly_trailed {
                // 1. Radial angle coherence (catches coma/field curvature):
                //    If elongation angles point radially from image center,
                //    it's optics, not trailing.
                let cx_img = width as f64 / 2.0;
                let cy_img = height as f64 / 2.0;
                let (rad_cos, rad_sin) = measured.iter()
                    .fold((0.0f64, 0.0f64), |(sc, ss), s| {
                        let phi = (s.y as f64 - cy_img).atan2(s.x as f64 - cx_img);
                        let delta = 2.0 * (s.theta as f64 - phi);
                        (sc + delta.cos(), ss + delta.sin())
                    });
                let nf = n as f64;
                let radial_r_sq = (rad_cos * rad_cos + rad_sin * rad_sin) / (nf * nf);
                let radial_p = (-nf * radial_r_sq).exp();
                if radial_r_sq > 0.15 && radial_p < 0.05 {
                    possibly_trailed = false;
                }
            }
            if possibly_trailed {
                // 2. Eccentricity-distance correlation (catches tilt):
                //    If eccentricity increases with distance from center,
                //    it's a field-dependent optical effect, not trailing.
                let cx_img = width as f32 / 2.0;
                let cy_img = height as f32 / 2.0;
                let nf = n as f64;
                let (sum_d, sum_e) = measured.iter()
                    .fold((0.0f64, 0.0f64), |(sd, se), s| {
                        let d = ((s.x - cx_img).powi(2) + (s.y - cy_img).powi(2)).sqrt();
                        (sd + d as f64, se + s.eccentricity as f64)
                    });
                let mean_d = sum_d / nf;
                let mean_e = sum_e / nf;
                let (mut cov, mut var_d, mut var_e) = (0.0f64, 0.0f64, 0.0f64);
                for s in measured.iter() {
                    let d = ((s.x - cx_img).powi(2) + (s.y - cy_img).powi(2)).sqrt() as f64;
                    let dd = d - mean_d;
                    let de = s.eccentricity as f64 - mean_e;
                    cov += dd * de;
                    var_d += dd * dd;
                    var_e += de * de;
                }
                let denom = (var_d * var_e).sqrt();
                if denom > 0.0 && cov / denom > 0.25 {
                    possibly_trailed = false;
                }
            }
        }

        // ── Stage 4: Metrics ─────────────────────────────────────────────
        let t = std::time::Instant::now();

        // FWHM & HFR: ecc ≤ 0.8 filter — elongated profiles inflate
        // geometric-mean FWHM. On trailed frames bypass it.
        const FWHM_ECC_MAX: f32 = 0.8;
        let fwhm_filtered: Vec<&metrics::MeasuredStar> = if possibly_trailed {
            measured.iter().collect()
        } else {
            let round: Vec<&metrics::MeasuredStar> = measured.iter()
                .filter(|s| s.eccentricity <= FWHM_ECC_MAX)
                .collect();
            if round.len() >= 3 { round } else { measured.iter().collect() }
        };
        let (fwhm_vals, hfr_vals, shape_weights) = (
            fwhm_filtered.iter().map(|s| s.fwhm).collect::<Vec<f32>>(),
            fwhm_filtered.iter().map(|s| s.hfr).collect::<Vec<f32>>(),
            fwhm_filtered.iter().map(|s| 1.0 / (1.0 + s.fit_residual)).collect::<Vec<f32>>(),
        );
        let median_fwhm = sigma_clipped_weighted_median(&fwhm_vals, &shape_weights);

        // Eccentricity: on normal frames, ecc ≤ 0.8 cutoff removes noise from
        // faint detections. On trailed frames, elongation IS the signal — bypass
        // the cutoff so the reported ecc reflects actual frame quality.
        let ecc_use_all = possibly_trailed;
        let ecc_filtered: Vec<&metrics::MeasuredStar> = if ecc_use_all {
            measured.iter().collect()
        } else {
            let filtered: Vec<&metrics::MeasuredStar> = measured.iter()
                .filter(|s| s.eccentricity <= FWHM_ECC_MAX)
                .collect();
            if filtered.len() >= 3 { filtered } else { measured.iter().collect() }
        };
        let ecc_vals: Vec<f32> = ecc_filtered.iter().map(|s| s.eccentricity).collect();
        let ecc_weights: Vec<f32> = ecc_filtered.iter()
            .map(|s| 1.0 / (1.0 + s.fit_residual))
            .collect();

        let statistics_ms_before_snr = t.elapsed().as_secs_f64() * 1000.0;

        let t = std::time::Instant::now();
        snr::compute_star_snr(&lum, width, height, &mut measured, median_fwhm);
        let snr_ms = t.elapsed().as_secs_f64() * 1000.0;

        let t = std::time::Instant::now();
        let mut snr_vals: Vec<f32> = measured.iter().map(|s| s.snr).collect();

        let median_eccentricity = sigma_clipped_weighted_median(&ecc_vals, &ecc_weights);
        let median_snr = find_median(&mut snr_vals);
        let median_hfr = sigma_clipped_weighted_median(&hfr_vals, &shape_weights);
        let psf_signal = snr::compute_psf_signal(&measured, bg_result.noise);
        let snr_weight = snr::compute_snr_weight(&measured, bg_result.background, bg_result.noise);

        // Median beta: use field_beta from calibration, or compute from all stars
        let median_beta = if let Some(fb) = field_beta {
            Some(fb as f32)
        } else {
            let mut beta_vals: Vec<f32> = measured.iter().filter_map(|s| s.beta).collect();
            if beta_vals.is_empty() { None } else { Some(find_median(&mut beta_vals)) }
        };

        let plate_scale = match (self.config.focal_length_mm, self.config.pixel_size_um) {
            (Some(fl), Some(ps)) if fl > 0.0 && ps > 0.0 => {
                Some((ps / fl * 206.265) as f32)
            }
            _ => None,
        };

        let median_fwhm_arcsec = plate_scale.map(|s| median_fwhm * s);
        let median_hfr_arcsec = plate_scale.map(|s| median_hfr * s);

        // Late cap: truncate to max_stars AFTER all statistics are computed
        measured.truncate(self.config.max_stars);

        let stars: Vec<StarMetrics> = measured
            .into_iter()
            .map(|m| StarMetrics {
                x: m.x, y: m.y, peak: m.peak, flux: m.flux,
                fwhm_x: m.fwhm_x, fwhm_y: m.fwhm_y, fwhm: m.fwhm,
                eccentricity: m.eccentricity, snr: m.snr, hfr: m.hfr,
                theta: m.theta, beta: m.beta, fit_method: m.fit_method,
                fit_residual: m.fit_residual,
                fwhm_arcsec: plate_scale.map(|s| m.fwhm * s),
                hfr_arcsec: plate_scale.map(|s| m.hfr * s),
            })
            .collect();
        let statistics_ms = statistics_ms_before_snr + t.elapsed().as_secs_f64() * 1000.0;
        let total_ms = pipeline_start.elapsed().as_secs_f64() * 1000.0;

        Ok(AnalysisResult {
            width, height, source_channels: channels,
            background: bg_result.background, noise: bg_result.noise,
            detection_threshold, stars_detected, stars,
            median_fwhm, median_eccentricity, median_snr, median_hfr,
            snr_weight, psf_signal, frame_snr,
            trail_r_squared, possibly_trailed,
            measured_fwhm_kernel: final_fwhm,
            median_beta,
            pass1_detections, calibrated_fwhm: field_fwhm,
            stars_measured, moffat_count, gaussian_count,
            plate_scale, median_fwhm_arcsec, median_hfr_arcsec,
            stage_timing: StageTiming {
                background_ms, detection_pass1_ms, calibration_ms,
                detection_pass2_ms, measurement_ms, snr_ms,
                statistics_ms, total_ms,
            },
        })
    }
}

impl Default for ImageAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

/// Sigma-clipped median: 2-iteration, 3σ MAD-based clipping.
///
/// Standard in SExtractor/DAOPHOT for robust statistics:
///   MAD = median(|x_i − median|)
///   σ_MAD = 1.4826 × MAD
///   reject: |x − median| > 3 × σ_MAD
///
/// Returns plain median if fewer than 3 values remain after clipping.
pub fn sigma_clipped_median(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut v: Vec<f32> = values.to_vec();
    for _ in 0..2 {
        if v.len() < 3 {
            break;
        }
        let med = find_median(&mut v);
        let mut abs_devs: Vec<f32> = v.iter().map(|&x| (x - med).abs()).collect();
        let mad = find_median(&mut abs_devs);
        let sigma_mad = 1.4826 * mad;
        if sigma_mad < 1e-10 {
            break; // all values identical
        }
        let clip = 3.0 * sigma_mad;
        v.retain(|&x| (x - med).abs() <= clip);
    }
    if v.is_empty() {
        return find_median(&mut values.to_vec());
    }
    find_median(&mut v)
}

/// Weighted median: walk sorted (value, weight) pairs until cumulative weight >= total/2.
///
/// Returns 0.0 if inputs are empty or total weight is zero.
pub fn weighted_median(values: &[f32], weights: &[f32]) -> f32 {
    if values.is_empty() || values.len() != weights.len() {
        return 0.0;
    }
    let mut pairs: Vec<(f32, f32)> = values.iter().copied()
        .zip(weights.iter().copied())
        .filter(|(_, w)| *w > 0.0)
        .collect();
    if pairs.is_empty() {
        return 0.0;
    }
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: f32 = pairs.iter().map(|(_, w)| w).sum();
    if total <= 0.0 {
        return 0.0;
    }
    let half = total * 0.5;
    let mut cumulative = 0.0_f32;
    for &(val, w) in &pairs {
        cumulative += w;
        if cumulative >= half {
            return val;
        }
    }
    pairs.last().unwrap().0
}

/// Sigma-clipped weighted median: 2-iteration 3σ MAD clipping, then weighted median.
///
/// Combines outlier rejection (via MAD) with continuous quality weighting.
/// Falls back to plain weighted median if fewer than 3 values survive clipping.
pub fn sigma_clipped_weighted_median(values: &[f32], weights: &[f32]) -> f32 {
    if values.is_empty() || values.len() != weights.len() {
        return 0.0;
    }
    let mut v: Vec<f32> = values.to_vec();
    let mut w: Vec<f32> = weights.to_vec();
    for _ in 0..2 {
        if v.len() < 3 {
            break;
        }
        let med = weighted_median(&v, &w);
        let abs_devs: Vec<f32> = v.iter().map(|&x| (x - med).abs()).collect();
        // Unweighted MAD for clipping threshold (weights affect median, not clip boundary)
        let mut sorted_devs = abs_devs.clone();
        sorted_devs.sort_by(|a, b| a.total_cmp(b));
        let mad = sorted_devs[sorted_devs.len() / 2];
        let sigma_mad = 1.4826 * mad;
        if sigma_mad < 1e-10 {
            break;
        }
        let clip = 3.0 * sigma_mad;
        let mut new_v = Vec::with_capacity(v.len());
        let mut new_w = Vec::with_capacity(w.len());
        for (val, wt) in v.iter().zip(w.iter()) {
            if (*val - med).abs() <= clip {
                new_v.push(*val);
                new_w.push(*wt);
            }
        }
        v = new_v;
        w = new_w;
    }
    if v.is_empty() {
        return weighted_median(values, weights);
    }
    weighted_median(&v, &w)
}

/// Estimate FWHM from the brightest detected stars by extracting stamps
/// and using `estimate_sigma_halfmax`. Returns median FWHM, or 0.0 if
/// fewer than 3 stars yield valid measurements.
pub fn estimate_fwhm_from_stars(
    lum: &[f32],
    width: usize,
    height: usize,
    stars: &[detection::DetectedStar],
    background: f32,
    bg_map: Option<&[f32]>,
) -> f32 {
    // Scan top 50 brightest (already sorted by flux descending), select up to 20
    // with low eccentricity (≤ 0.7) to avoid elongated non-stellar objects
    // poisoning the kernel estimate.
    let scan_n = stars.len().min(50);
    if scan_n < 3 {
        return 0.0;
    }

    let round_stars: Vec<&detection::DetectedStar> = stars[..scan_n]
        .iter()
        .filter(|s| s.eccentricity <= 0.7)
        .take(20)
        .collect();
    if round_stars.len() < 3 {
        return 0.0;
    }

    let mut fwhm_vals = Vec::with_capacity(round_stars.len());
    for star in &round_stars {
        let stamp_radius = 30_usize; // enough for FWHM up to ~25px (defocused)
        let cx = star.x.round() as i32;
        let cy = star.y.round() as i32;
        let sr = stamp_radius as i32;
        if cx - sr <= 0 || cy - sr <= 0
            || cx + sr >= width as i32 - 1
            || cy + sr >= height as i32 - 1
        {
            continue;
        }
        let x0 = (cx - sr) as usize;
        let y0 = (cy - sr) as usize;
        let x1 = (cx + sr) as usize;
        let y1 = (cy + sr) as usize;
        let stamp_w = x1 - x0 + 1;
        let mut stamp = Vec::with_capacity(stamp_w * (y1 - y0 + 1));
        for sy in y0..=y1 {
            for sx in x0..=x1 {
                let bg = bg_map.map_or(background, |m| m[sy * width + sx]);
                stamp.push(lum[sy * width + sx] - bg);
            }
        }
        let rel_cx = star.x - x0 as f32;
        let rel_cy = star.y - y0 as f32;
        let sigma = metrics::estimate_sigma_halfmax(&stamp, stamp_w, rel_cx, rel_cy);
        let fwhm = sigma * 2.3548;
        if fwhm > 1.0 && fwhm < 20.0 {
            fwhm_vals.push(fwhm);
        }
    }

    if fwhm_vals.len() < 3 {
        return 0.0;
    }
    find_median(&mut fwhm_vals)
}

/// Build a boolean mask marking green CFA pixel positions.
///
/// Returns a `Vec<bool>` of length `width * height` where `true` marks pixels
/// Extract luminance from planar RGB data: L = 0.2126R + 0.7152G + 0.0722B
pub fn extract_luminance(data: &[f32], width: usize, height: usize) -> Vec<f32> {
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

/// Prepare luminance data from raw metadata + pixels.
///
/// Handles u16→f32 conversion, green-channel interpolation for OSC,
/// and luminance extraction for multi-channel images.
/// Returns `(luminance, width, height, channels, green_mask)`.
#[cfg(feature = "debug-pipeline")]
pub fn prepare_luminance(
    meta: &crate::types::ImageMetadata,
    pixels: &crate::types::PixelData,
    apply_debayer: bool,
) -> (Vec<f32>, usize, usize, usize, Option<Vec<bool>>) {
    use crate::processing::color::u16_to_f32;
    use crate::processing::debayer;

    let f32_data = match pixels {
        PixelData::Float32(d) => std::borrow::Cow::Borrowed(d.as_slice()),
        PixelData::Uint16(d) => std::borrow::Cow::Owned(u16_to_f32(d)),
    };

    let mut data = f32_data.into_owned();
    let width = meta.width;
    let height = meta.height;
    let channels = meta.channels;

    // OSC green interpolation (matching Siril's interpolate_nongreen).
    // No green mask — PSF fitting uses all pixels in the interpolated image.
    if apply_debayer
        && meta.bayer_pattern != BayerPattern::None
        && channels == 1
    {
        data = debayer::interpolate_green_f32(&data, width, height, meta.bayer_pattern);
    }
    let green_mask: Option<Vec<bool>> = None;

    let lum = if channels == 3 {
        extract_luminance(&data, width, height)
    } else {
        data[..width * height].to_vec()
    };

    (lum, width, height, channels, green_mask)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weighted_median_equal_weights() {
        // Equal weights → same as unweighted median
        let vals = [1.0_f32, 3.0, 5.0, 7.0, 9.0];
        let wts = [1.0_f32; 5];
        let wm = weighted_median(&vals, &wts);
        assert!((wm - 5.0).abs() < 0.01, "Equal-weight median should be 5, got {}", wm);
    }

    #[test]
    fn test_weighted_median_skewed_weights() {
        // Heavy weight on low value should pull median down
        let vals = [1.0_f32, 10.0];
        let wts = [9.0_f32, 1.0]; // 90% weight on 1.0
        let wm = weighted_median(&vals, &wts);
        assert!((wm - 1.0).abs() < 0.01, "Skewed-weight median should be 1, got {}", wm);
    }

    #[test]
    fn test_weighted_median_empty() {
        let wm = weighted_median(&[], &[]);
        assert_eq!(wm, 0.0);
    }

    #[test]
    fn test_weighted_median_single() {
        let wm = weighted_median(&[42.0], &[1.0]);
        assert!((wm - 42.0).abs() < 0.01);
    }

    #[test]
    fn test_sigma_clipped_weighted_median_basic() {
        // With an outlier, sigma clipping should reject it
        let vals = [3.0_f32, 3.1, 3.0, 3.2, 3.0, 100.0]; // 100.0 is outlier
        let wts = [1.0_f32; 6];
        let scwm = sigma_clipped_weighted_median(&vals, &wts);
        assert!(scwm < 4.0, "Outlier should be clipped, got {}", scwm);
    }
}
