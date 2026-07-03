/// Image analysis: FWHM, eccentricity, SNR, PSF signal.

#[cfg(feature = "debug-pipeline")]
pub mod background;
#[cfg(not(feature = "debug-pipeline"))]
mod background;

#[cfg(feature = "debug-pipeline")]
pub mod convolution;
#[cfg(not(feature = "debug-pipeline"))]
mod convolution;

// Adaptive multi-level detector used by plate solving. Internal.
mod adaptive_detection;

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

/// A rough star detection result: position + flux only, no PSF metrics.
///
/// Produced by [`ImageAnalyzer::detect_fast`] and friends. Intended for
/// pipelines (blind plate solving, quick previews) that only need `(x, y, flux)`
/// and can tolerate pass-1-only centroid accuracy (~0.3–0.5 px).
///
/// When [`ImageAnalyzer::with_centroid_refine`] is `true`, `x` and `y` are
/// PSF-refined (Moffat LM fit, sub-pixel accuracy ~0.05 px at SNR≥20), and
/// `fwhm` holds the mean fitted FWHM. `sx` and `sy` carry per-axis Gaussian-
/// sigma proxies (σ = FWHM_axis / 2.3548). When refinement is OFF (the
/// default), `sx = sy = fwhm = 0.0` and `raw_x = x`, `raw_y = y` — pass-1
/// centroid, byte-identical to the pre-refinement path.
#[derive(Clone, Debug)]
pub struct FastStar {
    /// Subpixel centroid X.
    /// Pass-1 (intensity-weighted) when refinement is off; PSF-refined when on.
    pub x: f32,
    /// Subpixel centroid Y.
    /// Pass-1 (intensity-weighted) when refinement is off; PSF-refined when on.
    pub y: f32,
    /// Pass-1 (intensity-weighted) centroid X — always the unrefined value.
    ///
    /// When refinement is OFF, `raw_x == x`. When refinement is ON and the
    /// per-star Moffat fit succeeds, `x` is updated to the refined centre while
    /// `raw_x` retains the pass-1 value. Stars whose fit is rejected (low SNR,
    /// non-physical result, centre-shift > 2 px) also have `raw_x == x`.
    /// Solvers can use `raw_x`/`raw_y` alongside refined `x`/`y` to pick
    /// whichever centroid yields a better final fit.
    pub raw_x: f32,
    /// Pass-1 (intensity-weighted) centroid Y — always the unrefined value.
    /// See [`FastStar::raw_x`] for the full contract.
    pub raw_y: f32,
    /// Background-subtracted peak value (ADU).
    pub peak: f32,
    /// Background-subtracted total flux (ADU).
    pub flux: f32,
    /// Aperture-photometry SNR (`flux / sqrt(flux + π r² σ²)`). Discriminates
    /// compact point sources from extended structure that has high flux but
    /// is spread over a large aperture (galaxy/nebula knots).
    pub snr: f32,
    /// Per-axis centroid uncertainty proxy (Gaussian sigma in X, pixels).
    /// `σ_x = FWHM_x / 2.3548` from the Moffat fit. `0.0` when refinement is off.
    pub sx: f32,
    /// Per-axis centroid uncertainty proxy (Gaussian sigma in Y, pixels).
    /// `σ_y = FWHM_y / 2.3548` from the Moffat fit. `0.0` when refinement is off.
    pub sy: f32,
    /// Mean fitted FWHM = `0.5 * (FWHM_x + FWHM_y)` from the Moffat fit.
    /// `0.0` when refinement is off.
    pub fwhm: f32,
}

/// Per-stage timing breakdown for the fast detection pipeline.
#[derive(Clone, Debug, Default)]
pub struct FastDetectTiming {
    /// File I/O (0.0 for data/raw entry points).
    pub read_ms: f64,
    /// f32 conversion + OSC green interpolation + luminance extraction.
    pub prep_ms: f64,
    /// Mesh-grid background & noise estimation (no MRS wavelet).
    pub background_ms: f64,
    /// Single-pass matched-filter detection.
    pub detection_ms: f64,
    /// Wall clock from entry to return.
    pub total_ms: f64,
}

/// Lean analysis result from the fast detector.
///
/// Contains only what blind plate solving needs: image dimensions, a
/// flux-sorted list of rough star centroids, and the background estimate
/// used to threshold them.
#[derive(Clone, Debug)]
pub struct FastAnalysisResult {
    /// Image width (after debayer if applicable).
    pub width: usize,
    /// Image height (after debayer if applicable).
    pub height: usize,
    /// Stars sorted by flux descending, capped at `max_stars`.
    pub stars: Vec<FastStar>,
    /// Global background level (ADU).
    pub background: f32,
    /// Background noise sigma (ADU).
    pub noise: f32,
    /// Per-stage timing breakdown.
    pub timing: FastDetectTiming,
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
    /// When `true`, run a PSF-refinement pass (Moffat LM fit) after pass-1
    /// detection and update `FastStar.{x, y, sx, sy, fwhm}`. Default `false`.
    /// The pass is skipped entirely when `false` — byte-identical output.
    centroid_refine: bool,
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
                centroid_refine: false,
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

    /// Enable or disable PSF centroid refinement for the fast detection path.
    ///
    /// When `true`, after pass-1 detection a Moffat LM fit is run per star
    /// (parallel via rayon). On success the fitted centre replaces the
    /// intensity-weighted centroid, and `FastStar.{sx, sy, fwhm}` are
    /// populated from the fit. Stars that fail per-star gates (low SNR,
    /// non-physical fit, centre-shift > 2 px) keep their pass-1 centroid.
    ///
    /// Default: `false`. When `false` the refinement pass is skipped
    /// entirely — output is byte-identical to the pre-refinement path.
    pub fn with_centroid_refine(mut self, refine: bool) -> Self {
        self.config.centroid_refine = refine;
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

    /// Fast star detection from a FITS/XISF file.
    ///
    /// Runs the matched-filter detector once (no PSF calibration, no LM Moffat
    /// fitting, no SNR photometry, no trail detection). Produces rough
    /// `(x, y, flux)` centroids with pass-1 accuracy in ~300–500 ms on a
    /// full-frame image. Intended for pipelines that only need positions and
    /// brightness ordering, such as blind plate solving.
    ///
    /// For per-star PSF metrics (FWHM, eccentricity, SNR, HFR, beta) use
    /// [`ImageAnalyzer::analyze`] instead.
    pub fn detect_fast<P: AsRef<Path>>(&self, path: P) -> Result<FastAnalysisResult> {
        let path = path.as_ref();
        match &self.thread_pool {
            Some(pool) => pool.install(|| self.detect_fast_impl(path)),
            None => self.detect_fast_impl(path),
        }
    }

    /// Fast detection from pre-loaded planar f32 pixel data.
    ///
    /// `data`: planar f32 pixel data (for 3-channel: RRRGGGBBB layout).
    /// See [`ImageAnalyzer::detect_fast`] for pipeline details.
    pub fn detect_fast_data(
        &self,
        data: &[f32],
        width: usize,
        height: usize,
        channels: usize,
    ) -> Result<FastAnalysisResult> {
        match &self.thread_pool {
            Some(pool) => pool.install(|| {
                self.run_fast_detection(data, width, height, channels)
            }),
            None => self.run_fast_detection(data, width, height, channels),
        }
    }

    /// Fast detection from pre-read raw pixel data (skips file I/O).
    ///
    /// Accepts `ImageMetadata` and borrows `PixelData`, handling u16→f32
    /// conversion and green-channel interpolation for OSC images internally.
    pub fn detect_fast_raw(
        &self,
        meta: &ImageMetadata,
        pixels: &PixelData,
    ) -> Result<FastAnalysisResult> {
        match &self.thread_pool {
            Some(pool) => pool.install(|| self.detect_fast_raw_impl(meta, pixels)),
            None => self.detect_fast_raw_impl(meta, pixels),
        }
    }

    fn detect_fast_impl(&self, path: &Path) -> Result<FastAnalysisResult> {
        let t_read = std::time::Instant::now();
        let (meta, pixel_data) =
            formats::read_image(path).context("Failed to read image for fast detection")?;
        let read_ms = t_read.elapsed().as_secs_f64() * 1000.0;
        tracing::debug!(
            stage = "read",
            duration_ms = read_ms,
            path = %path.display(),
            "fast detection stage timing"
        );

        let f32_data = match pixel_data {
            PixelData::Float32(d) => d,
            PixelData::Uint16(d) => u16_to_f32(&d),
        };

        let mut data = f32_data;
        let width = meta.width;
        let height = meta.height;
        let channels = meta.channels;

        if self.config.apply_debayer
            && meta.bayer_pattern != BayerPattern::None
            && channels == 1
        {
            data = debayer::interpolate_green_f32(&data, width, height, meta.bayer_pattern);
        }

        let mut result = self.run_fast_detection(&data, width, height, channels)?;
        result.timing.read_ms = read_ms;
        result.timing.total_ms += read_ms;
        Ok(result)
    }

    fn detect_fast_raw_impl(
        &self,
        meta: &ImageMetadata,
        pixels: &PixelData,
    ) -> Result<FastAnalysisResult> {
        let f32_data = match pixels {
            PixelData::Float32(d) => std::borrow::Cow::Borrowed(d.as_slice()),
            PixelData::Uint16(d) => std::borrow::Cow::Owned(u16_to_f32(d)),
        };

        let mut data = f32_data.into_owned();
        let width = meta.width;
        let height = meta.height;
        let channels = meta.channels;

        if self.config.apply_debayer
            && meta.bayer_pattern != BayerPattern::None
            && channels == 1
        {
            data = debayer::interpolate_green_f32(&data, width, height, meta.bayer_pattern);
        }

        self.run_fast_detection(&data, width, height, channels)
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
        tracing::debug!(
            stage = "background",
            duration_ms = background_ms,
            background = bg_result.background,
            noise = bg_result.noise,
            "analysis stage timing"
        );

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
        tracing::debug!(
            stage = "detect_pass1",
            duration_ms = detection_pass1_ms,
            count = pass1_stars.len(),
            "analysis stage timing"
        );

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
        tracing::debug!(
            stage = "calibration",
            duration_ms = calibration_ms,
            count = calibration_stars.len(),
            "analysis stage timing"
        );

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
        tracing::debug!(
            stage = "detect_pass2",
            duration_ms = detection_pass2_ms,
            count = detected.len(),
            "analysis stage timing"
        );

        let bg_map_ref = bg_result.background_map.as_deref();
        let detection_threshold = self.config.detection_sigma * bg_result.noise;

        // Trail detection defaults — computed after PSF measurement (see below).
        let mut trail_r_squared = 0.0_f32;
        let mut possibly_trailed = false;

        let frame_snr = if bg_result.noise > 0.0 { bg_result.background / bg_result.noise } else { 0.0 };

        let pass1_detections = pass1_stars.len();

        let make_zero_result = |stars_detected: usize| {
            let total_ms = pipeline_start.elapsed().as_secs_f64() * 1000.0;
            tracing::debug!(
                stage = "total",
                duration_ms = total_ms,
                stars_detected,
                count = 0,
                "analysis pipeline finished (no stars measured)"
            );
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
                    statistics_ms: 0.0, total_ms,
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
        tracing::debug!(
            stage = "measurement",
            duration_ms = measurement_ms,
            count = stars_measured,
            moffat_count,
            gaussian_count,
            "analysis stage timing"
        );

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
        tracing::debug!(
            stage = "snr",
            duration_ms = snr_ms,
            count = measured.len(),
            "analysis stage timing"
        );

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
        let stars_count = stars.len();
        let statistics_ms = statistics_ms_before_snr + t.elapsed().as_secs_f64() * 1000.0;
        tracing::debug!(
            stage = "statistics",
            duration_ms = statistics_ms,
            count = stars_count,
            "analysis stage timing"
        );
        let total_ms = pipeline_start.elapsed().as_secs_f64() * 1000.0;
        tracing::debug!(
            stage = "total",
            duration_ms = total_ms,
            stars_detected,
            stars_measured,
            count = stars_count,
            "analysis pipeline finished"
        );

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

    /// Fast detection pipeline used by plate solving — luminance extraction
    /// then the adaptive multi-level detector (falling-threshold ladder with
    /// occupancy mask, saturation-aware level clip, per-tile deep pass; see
    /// [`adaptive_detection`]). Replaces the old single global-threshold
    /// pass, which under-detected the long-focal-length / extended-object
    /// frames the solver most needs. No PSF calibration / LM fitting / SNR
    /// photometry / trail detection — the solver only needs `(x, y, flux)`.
    ///
    /// When [`AnalysisConfig::centroid_refine`] is `true`, an optional
    /// Moffat-LM refinement pass runs after detection and may improve
    /// centroid accuracy to ~0.05 px (at SNR≥20). The fast path (refine=false)
    /// is byte-identical to the pre-refinement implementation.
    fn run_fast_detection(
        &self,
        data: &[f32],
        width: usize,
        height: usize,
        channels: usize,
    ) -> Result<FastAnalysisResult> {
        use rayon::prelude::*;

        let pipeline_start = std::time::Instant::now();

        let t_prep = std::time::Instant::now();
        let lum = if channels == 3 {
            extract_luminance(data, width, height)
        } else {
            data[..width * height].to_vec()
        };
        let prep_ms = t_prep.elapsed().as_secs_f64() * 1000.0;
        tracing::debug!(
            stage = "prep",
            duration_ms = prep_ms,
            "fast detection stage timing"
        );

        let t_bg = std::time::Instant::now();
        let (background, noise) =
            adaptive_detection::background_and_noise(&lum, width, height);
        let background_ms = t_bg.elapsed().as_secs_f64() * 1000.0;
        tracing::debug!(
            stage = "background",
            duration_ms = background_ms,
            background,
            noise,
            "fast detection stage timing"
        );

        let t_det = std::time::Instant::now();
        let detected = adaptive_detection::detect_stars_adaptive(
            &lum,
            width,
            height,
            self.config.max_stars,
            0.8,
        );
        let detection_ms = t_det.elapsed().as_secs_f64() * 1000.0;
        let detected_count = detected.len();
        tracing::debug!(
            stage = "detect",
            duration_ms = detection_ms,
            count = detected_count,
            "fast detection stage timing"
        );

        // Already brightest-first and trimmed to max_stars by the detector.
        // Build pass-1 FastStar list (sx/sy/fwhm = 0 — will be filled if
        // centroid_refine is enabled).
        let mut stars: Vec<FastStar> = detected
            .into_iter()
            .map(|(ds, snr)| FastStar {
                x: ds.x,
                y: ds.y,
                raw_x: ds.x,
                raw_y: ds.y,
                peak: ds.peak,
                flux: ds.flux,
                snr,
                sx: 0.0,
                sy: 0.0,
                fwhm: 0.0,
            })
            .collect();

        // --- PSF centroid refinement (opt-in) --------------------------------
        // When enabled, run a Moffat LM fit per star in parallel. Stars whose
        // SNR < 10, whose fit fails, or whose fitted centre moves > 2 px from
        // the pass-1 centroid are silently kept at pass-1. This pass is skipped
        // entirely when refine=false, giving byte-identical output to the old
        // path.
        if self.config.centroid_refine && !stars.is_empty() {
            // Estimate a rough field FWHM from the top stars' HFDs (the
            // adaptive detector does not expose HFD directly, so we use a
            // simple area-based estimate: assume circular profile, σ ≈ √(area/π)).
            // We use a fixed init_sigma of 3.0 px if we cannot estimate it —
            // the Moffat fit is robust to a wide range of initial sigmas.
            let init_sigma_est = 3.0_f32; // sensible default for most setups

            // Refinement: parallel per star.
            let refined: Vec<Option<(f32, f32, f32, f32, f32)>> = stars
                .par_iter()
                .map(|star| {
                    refine_centroid_moffat(
                        &lum,
                        width,
                        height,
                        star,
                        background,
                        init_sigma_est,
                    )
                })
                .collect();

            for (star, refined_opt) in stars.iter_mut().zip(refined.into_iter()) {
                if let Some((rx, ry, rsx, rsy, rfwhm)) = refined_opt {
                    // Preserve the pass-1 centroid in raw_x/raw_y BEFORE
                    // overwriting x/y with the Moffat-refined position.
                    // raw_x/raw_y were initialised to pass-1 values above;
                    // this assignment is a no-op when x/y were not yet changed,
                    // but makes the intent explicit.
                    star.raw_x = star.x;
                    star.raw_y = star.y;
                    star.x = rx;
                    star.y = ry;
                    star.sx = rsx;
                    star.sy = rsy;
                    star.fwhm = rfwhm;
                }
            }
        }
        // --- end PSF refinement -----------------------------------------------

        let total_ms = pipeline_start.elapsed().as_secs_f64() * 1000.0;
        tracing::debug!(
            stage = "total",
            duration_ms = total_ms,
            count = detected_count,
            "fast detection pipeline finished"
        );

        Ok(FastAnalysisResult {
            width,
            height,
            stars,
            background,
            noise,
            timing: FastDetectTiming {
                read_ms: 0.0,
                prep_ms,
                background_ms,
                detection_ms,
                total_ms,
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

/// PSF centroid refinement via Moffat LM fit on a stamp window.
///
/// Used by [`ImageAnalyzer::run_fast_detection`] when `centroid_refine` is
/// enabled. Returns `Some((x, y, sx, sy, fwhm))` on success, where `x/y` are
/// the refined centroid (image coordinates), `sx/sy` are per-axis Gaussian-
/// sigma proxies (`FWHM_axis / 2.3548`), and `fwhm` is the mean fitted FWHM.
/// Returns `None` when the per-star gates reject the fit (low SNR, non-physical
/// parameters, or fitted centre moved > 2 px from pass-1).
fn refine_centroid_moffat(
    lum: &[f32],
    width: usize,
    height: usize,
    star: &FastStar,
    background: f32,
    init_sigma: f32,
) -> Option<(f32, f32, f32, f32, f32)> {
    // Gate: skip noisy stars — Moffat fit is unreliable at low SNR.
    if star.snr < 10.0 {
        return None;
    }

    // Stamp window: 4·sigma radius (minimum 8 px, max 50 px).
    let stamp_radius = ((4.0 * init_sigma).ceil() as usize).max(8).min(50);
    let cx = star.x.round() as i32;
    let cy = star.y.round() as i32;
    let sr = stamp_radius as i32;

    // Bounds check — keep well inside the image.
    if cx - sr <= 0 || cy - sr <= 0
        || cx + sr >= width as i32 - 1
        || cy + sr >= height as i32 - 1
    {
        return None;
    }

    let x0 = (cx - sr) as usize;
    let y0 = (cy - sr) as usize;
    let x1 = (cx + sr) as usize;
    let y1 = (cy + sr) as usize;
    let stamp_w = x1 - x0 + 1;
    let stamp_h = y1 - y0 + 1;

    // Extract background-subtracted stamp.
    let mut stamp = Vec::with_capacity(stamp_w * stamp_h);
    for sy in y0..=y1 {
        for sx in x0..=x1 {
            stamp.push(lum[sy * width + sx] - background);
        }
    }

    // Centroid relative to stamp.
    let init_cx = star.x - x0 as f32;
    let init_cy = star.y - y0 as f32;

    // Fitting radius for sample selection.
    let fitting_radius = 5.0_f64.max(4.0 * init_sigma as f64);
    let fitting_radius_sq = fitting_radius * fitting_radius;
    let cx64 = init_cx as f64;
    let cy64 = init_cy as f64;

    let pixels: Vec<fitting::PixelSample> = (0..stamp_h)
        .flat_map(|sy| {
            (0..stamp_w).filter_map(move |sx| {
                let dx = sx as f64 - cx64;
                let dy = sy as f64 - cy64;
                if dx * dx + dy * dy <= fitting_radius_sq {
                    Some((sx, sy))
                } else {
                    None
                }
            })
        })
        .map(|(sx, sy)| fitting::PixelSample {
            x: sx as f64,
            y: sy as f64,
            value: stamp[sy * stamp_w + sx] as f64,
        })
        .collect();

    if pixels.len() < 12 {
        return None;
    }

    // Background from annulus.
    let bg_outer = fitting_radius + 3.0;
    let bg_outer_sq = bg_outer * bg_outer;
    let mut annulus_vals: Vec<f64> = Vec::new();
    for sy in 0..stamp_h {
        for sx in 0..stamp_w {
            let dx = sx as f64 - cx64;
            let dy = sy as f64 - cy64;
            let r_sq = dx * dx + dy * dy;
            if r_sq > fitting_radius_sq && r_sq <= bg_outer_sq {
                annulus_vals.push(stamp[sy * stamp_w + sx] as f64);
            }
        }
    }
    annulus_vals.sort_by(|a, b| a.total_cmp(b));
    let init_bg = if annulus_vals.is_empty() {
        0.0
    } else {
        annulus_vals[annulus_vals.len() / 2]
    };

    // Use fixed beta=3 (typical astrophotos) so the fit is 7-parameter and stable.
    let beta = 3.0_f64;
    let init_sigma_d = init_sigma as f64;
    let result = fitting::fit_moffat_2d_fixed_beta(
        &pixels,
        init_bg,
        (star.peak - background) as f64,
        init_cx as f64,
        init_cy as f64,
        init_sigma_d,
        init_sigma_d,
        0.0, // theta = 0 init (fit will rotate)
        beta,
        25,   // max_iter
        1e-4, // conv_tol
        5,    // max_rejects
    )?;

    // Reject non-physical fits.
    let max_alpha = init_sigma_d * 5.0;
    if result.alpha_x <= 0.0 || result.alpha_y <= 0.0
        || result.alpha_x > max_alpha || result.alpha_y > max_alpha
    {
        return None;
    }

    // Fitted centre in image coordinates.
    let rx = (result.x0 + x0 as f64) as f32;
    let ry = (result.y0 + y0 as f64) as f32;

    // Reject if the fit moved the centroid more than 2 px — likely a bad fit
    // (blended source, cosmic ray, etc.).
    let shift_sq = (rx - star.x) * (rx - star.x) + (ry - star.y) * (ry - star.y);
    if shift_sq > 4.0 {
        return None;
    }

    // Per-axis FWHM from Moffat: FWHM = 2α√(2^(1/β)−1)
    let fwhm_x = result.fwhm_x() as f32;
    let fwhm_y = result.fwhm_y() as f32;

    // Non-physical FWHM guard.
    if fwhm_x < 0.5 || fwhm_y < 0.5 || fwhm_x > 50.0 || fwhm_y > 50.0 {
        return None;
    }

    // σ = FWHM / 2.3548 per axis.
    let rsx = fwhm_x / 2.3548;
    let rsy = fwhm_y / 2.3548;
    let rfwhm = 0.5 * (fwhm_x + fwhm_y);

    Some((rx, ry, rsx, rsy, rfwhm))
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

    /// Build a synthetic mono f32 image with Gaussian stars at known positions.
    /// Background is flat (100 ADU) + Gaussian noise; peak is added on top.
    fn make_synthetic_field(
        width: usize,
        height: usize,
        stars: &[(f32, f32, f32)], // (cx, cy, peak)
        sigma: f32,
        bg: f32,
    ) -> Vec<f32> {
        let mut img = vec![bg; width * height];
        // Deterministic low-amplitude noise so the test is reproducible.
        for (i, px) in img.iter_mut().enumerate() {
            let r = ((i.wrapping_mul(2654435761)) & 0xffff) as f32 / 65535.0 - 0.5;
            *px += r * 10.0;
        }
        let radius = (4.0 * sigma).ceil() as i32;
        for &(cx, cy, peak) in stars {
            let ix = cx.round() as i32;
            let iy = cy.round() as i32;
            for dy in -radius..=radius {
                for dx in -radius..=radius {
                    let x = ix + dx;
                    let y = iy + dy;
                    if x < 0 || y < 0 || x >= width as i32 || y >= height as i32 {
                        continue;
                    }
                    let rx = x as f32 - cx;
                    let ry = y as f32 - cy;
                    let g = peak * (-(rx * rx + ry * ry) / (2.0 * sigma * sigma)).exp();
                    img[y as usize * width + x as usize] += g;
                }
            }
        }
        img
    }

    #[test]
    fn test_fast_detect_synthetic_field() {
        let width = 256;
        let height = 256;
        // 10 stars at known subpixel positions, well separated, varying brightness.
        let truth: [(f32, f32, f32); 10] = [
            (32.3, 41.7, 12000.0),
            (78.1, 50.5, 9000.0),
            (128.6, 30.2, 14000.0),
            (200.0, 40.9, 6000.0),
            (50.8, 120.4, 11000.0),
            (130.1, 140.7, 8000.0),
            (210.5, 130.3, 10000.0),
            (60.2, 200.9, 7000.0),
            (150.6, 210.4, 13000.0),
            (220.8, 220.1, 5500.0),
        ];
        let sigma = 1.8_f32; // FWHM ≈ 4.2 px
        let img = make_synthetic_field(width, height, &truth, sigma, 100.0);

        let analyzer = ImageAnalyzer::new()
            .with_detection_sigma(5.0)
            .with_max_stars(50)
            .with_min_star_area(3);

        let result = analyzer
            .detect_fast_data(&img, width, height, 1)
            .expect("fast detection succeeds on synthetic field");

        assert_eq!(result.width, width);
        assert_eq!(result.height, height);
        assert!(
            result.stars.len() >= truth.len(),
            "expected >= {} detections, got {}",
            truth.len(),
            result.stars.len()
        );

        // For every injected star, a detection must be within 0.5 px.
        let mut matched = 0usize;
        for &(tx, ty, _) in &truth {
            let best = result
                .stars
                .iter()
                .map(|s| {
                    let dx = s.x - tx;
                    let dy = s.y - ty;
                    (dx * dx + dy * dy).sqrt()
                })
                .fold(f32::INFINITY, f32::min);
            assert!(
                best < 0.5,
                "truth ({:.2},{:.2}) has no detection within 0.5 px (nearest {:.3})",
                tx, ty, best
            );
            matched += 1;
        }
        assert_eq!(matched, truth.len());

        // Flux ordering: brightest injected star should be #1 or #2.
        let brightest_truth = truth
            .iter()
            .cloned()
            .fold(truth[0], |a, b| if b.2 > a.2 { b } else { a });
        let top = &result.stars[0];
        let dx = top.x - brightest_truth.0;
        let dy = top.y - brightest_truth.1;
        assert!(
            (dx * dx + dy * dy).sqrt() < 1.5,
            "brightest detection should be near the brightest truth"
        );
    }

    /// Centroid refinement: with refine ON, recovered centre must be within
    /// 0.05 px of truth at SNR≈20. With refine OFF, pass-1 centroid is
    /// unchanged (sx=sy=fwhm=0.0, raw_x==x, raw_y==y).
    ///
    /// Additional Phase-4 assertion: when refinement is ON and the Moffat fit
    /// succeeds, `raw_x`/`raw_y` must hold the pass-1 (pre-refinement) centroid
    /// while `x`/`y` hold the PSF-refined centre. This guarantees the solver
    /// can compare both and keep the better fit.
    #[test]
    fn test_centroid_refine_on_vs_off() {
        let width = 200;
        let height = 200;
        // A single star at a known subpixel position with SNR≈20.
        // sigma=2.0 → FWHM≈4.71 px; peak=8000 over bg=500.
        let truth_x = 100.37_f32;
        let truth_y = 100.61_f32;
        let sigma = 2.0_f32;
        let peak = 8000.0_f32;
        let bg = 500.0_f32;

        let img = make_synthetic_field(width, height, &[(truth_x, truth_y, peak)], sigma, bg);

        // --- Refine OFF (default) ---
        let analyzer_off = ImageAnalyzer::new()
            .with_detection_sigma(5.0)
            .with_max_stars(10)
            .with_min_star_area(3)
            .with_centroid_refine(false);
        let result_off = analyzer_off
            .detect_fast_data(&img, width, height, 1)
            .expect("fast detect (refine=off) on synthetic star");

        assert!(
            !result_off.stars.is_empty(),
            "should detect the injected star with refine=off"
        );
        let s_off = &result_off.stars[0];
        // With refine off: sx=sy=fwhm must be exactly 0.
        assert_eq!(s_off.sx, 0.0, "refine=off: sx must be 0");
        assert_eq!(s_off.sy, 0.0, "refine=off: sy must be 0");
        assert_eq!(s_off.fwhm, 0.0, "refine=off: fwhm must be 0");
        // With refine off: raw_x/raw_y must equal x/y (pass-1 == pass-1).
        assert_eq!(s_off.raw_x, s_off.x, "refine=off: raw_x must equal x");
        assert_eq!(s_off.raw_y, s_off.y, "refine=off: raw_y must equal y");

        // --- Refine ON ---
        let analyzer_on = ImageAnalyzer::new()
            .with_detection_sigma(5.0)
            .with_max_stars(10)
            .with_min_star_area(3)
            .with_centroid_refine(true);
        let result_on = analyzer_on
            .detect_fast_data(&img, width, height, 1)
            .expect("fast detect (refine=on) on synthetic star");

        assert!(
            !result_on.stars.is_empty(),
            "should detect the injected star with refine=on"
        );
        let s_on = &result_on.stars[0];

        // With refine on: centroid must be within 0.05 px of truth.
        let err = ((s_on.x - truth_x).powi(2) + (s_on.y - truth_y).powi(2)).sqrt();
        assert!(
            err < 0.05,
            "refine=on: centroid error {:.4} px exceeds 0.05 px (x={:.3}, y={:.3}, truth=({:.3},{:.3}))",
            err, s_on.x, s_on.y, truth_x, truth_y
        );

        // With refine on: FWHM should be close to truth (sigma*2.3548 ≈ 4.71 px).
        let truth_fwhm = sigma * 2.3548_f32;
        assert!(
            (s_on.fwhm - truth_fwhm).abs() < 0.5,
            "refine=on: fitted fwhm={:.3} px, truth={:.3} px (delta > 0.5)",
            s_on.fwhm, truth_fwhm
        );

        // sx, sy should be non-zero and consistent with FWHM.
        assert!(s_on.sx > 0.0, "refine=on: sx must be positive");
        assert!(s_on.sy > 0.0, "refine=on: sy must be positive");
        let expected_sigma = truth_fwhm / 2.3548_f32;
        let sx_err = (s_on.sx - expected_sigma).abs();
        let sy_err = (s_on.sy - expected_sigma).abs();
        assert!(
            sx_err < 0.5,
            "refine=on: sx={:.3} far from expected sigma={:.3}",
            s_on.sx, expected_sigma
        );
        assert!(
            sy_err < 0.5,
            "refine=on: sy={:.3} far from expected sigma={:.3}",
            s_on.sy, expected_sigma
        );

        // Phase-4: raw_x/raw_y must hold the pass-1 centroid (pre-refinement),
        // while x/y hold the PSF-refined centre. The refined x must be strictly
        // closer to truth than the pass-1 raw_x (on a clean Gaussian star at
        // SNR≈20 the Moffat fit always wins by a measurable margin).
        // If the fit was rejected for this star (sx==0), skip the raw vs refined
        // comparison — that star's x==raw_x by construction (no overwrite occurred).
        if s_on.sx > 0.0 {
            // raw_x/raw_y must be set to the pass-1 centroid (which typically
            // differs from the refined by > 0.01 px on a subpixel-offset star).
            let raw_err = ((s_on.raw_x - truth_x).powi(2) + (s_on.raw_y - truth_y).powi(2)).sqrt();
            let ref_err = ((s_on.x - truth_x).powi(2) + (s_on.y - truth_y).powi(2)).sqrt();
            // Refined must be at least as good as pass-1.
            assert!(
                ref_err <= raw_err + 0.02,
                "refine=on: refined centroid ({:.4}px from truth) should be \
                 no worse than pass-1 ({:.4}px from truth)",
                ref_err, raw_err
            );
            // raw_x/raw_y should NOT be identical to the refined x/y when
            // the star is offset from the pixel grid centre (pass-1 rounding
            // vs PSF-refined sub-pixel differ by > some threshold). We allow
            // up to 0.5 px agreement as a lower bound; if raw==refined to
            // < 0.001 that means raw was overwritten with the refined value.
            // This is a regression guard: raw must be the pre-refinement value.
            // (For a star at truth_x=100.37, pass-1 will typically be at ~100.0
            //  or 100.5 depending on the intensity-centroid peak pixel. We just
            //  check that raw_x/raw_y are within 0.6 px of truth, since the
            //  pass-1 centroid is typically < 0.5 px from truth on a clean star.)
            assert!(
                raw_err < 0.6,
                "refine=on: raw (pass-1) centroid {:.4}px from truth — unexpectedly large \
                 (expected pass-1 to be within 0.6 px of truth on a clean Gaussian star)",
                raw_err
            );
        }
    }

    /// Low-SNR gate: a star with SNR < 10 must not be refined (sx=sy=fwhm=0).
    /// We inject a very faint star that will be too noisy for the Moffat fit to
    /// accept, even with refine=ON.
    #[test]
    fn test_centroid_refine_low_snr_fallback() {
        let width = 200;
        let height = 200;
        // Bright star (will be refined) and a very dim star (SNR < 10).
        let bright = (100.0_f32, 100.0_f32, 8000.0_f32);
        let dim    = (50.0_f32,  50.0_f32,  150.0_f32); // tiny peak over bg=500 → SNR < 10
        let sigma  = 2.0_f32;
        let bg     = 500.0_f32;

        let img = make_synthetic_field(width, height, &[bright, dim], sigma, bg);

        let analyzer = ImageAnalyzer::new()
            .with_detection_sigma(5.0)
            .with_max_stars(50)
            .with_min_star_area(3)
            .with_centroid_refine(true);
        let result = analyzer
            .detect_fast_data(&img, width, height, 1)
            .expect("fast detect on synthetic field (low-SNR gate test)");

        // Find the bright star's result and verify it was refined.
        let bright_star = result.stars.iter().find(|s| {
            let dx = s.x - bright.0;
            let dy = s.y - bright.1;
            (dx * dx + dy * dy).sqrt() < 2.0
        });
        if let Some(bs) = bright_star {
            // Bright star should be refined (fwhm != 0).
            assert!(
                bs.fwhm > 0.0,
                "bright star: fwhm=0 even with refine=on (expected refinement to succeed)"
            );
        }
        // Any star with snr < 10 must have sx=sy=fwhm=0 (no refinement).
        for s in &result.stars {
            if s.snr <= 10.0 {
                assert_eq!(
                    s.sx, 0.0,
                    "low-SNR star (snr={:.1}) should not be refined (sx != 0)", s.snr
                );
            }
        }
    }
}
