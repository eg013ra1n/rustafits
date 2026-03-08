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
    /// Total stars with valid PSF measurements (statistics computed from all).
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
    /// Image-wide SNR in decibels: 20 × log10(mean_signal / noise).
    pub snr_db: f32,
    /// SNR weight for frame ranking: (MeanDev / noise)².
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
    /// Measured FWHM from adaptive two-pass detection (pixels).
    /// This is the FWHM used for the final matched filter kernel.
    /// If the first-pass FWHM was within 30% of 3.0, this equals 3.0.
    pub measured_fwhm_kernel: f32,
    /// Median Moffat β across all stars (None if Moffat fitting not used).
    /// Typical range: 2.0-5.0 for real optics. Lower = broader wings.
    pub median_beta: Option<f32>,
}

/// Builder configuration for analysis (internal).
pub struct AnalysisConfig {
    detection_sigma: f32,
    min_star_area: usize,
    max_star_area: usize,
    saturation_fraction: f32,
    max_measure: Option<usize>,
    max_stars: usize,
    use_gaussian_fit: bool,
    background_mesh_size: Option<usize>,
    apply_debayer: bool,
    trail_r_squared_threshold: f32,
    use_moffat_fit: bool,
    iterative_background: usize,
    /// MRS wavelet noise layers: 0 = legacy MAD (default), 1+ = MRS wavelet.
    noise_layers: usize,
    /// Fixed Moffat beta: None = free (default), Some(b) = fixed.
    moffat_beta: Option<f32>,
    /// Maximum distortion (eccentricity) filter: None = no filter (default).
    max_distortion: Option<f32>,
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
                max_measure: None,
                max_stars: 200,
                use_gaussian_fit: true,
                background_mesh_size: None,
                apply_debayer: true,
                trail_r_squared_threshold: 0.5,
                use_moffat_fit: true,
                iterative_background: 1,
                noise_layers: 0,
                moffat_beta: None,
                max_distortion: None,
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

    /// Limit the number of stars measured for PSF fitting.
    /// Brightest N detected stars are measured; rest are skipped.
    /// Statistics are computed from all measured stars.
    /// Default: None (measure all). E.g. 5000 for faster dense fields.
    pub fn with_max_measure(mut self, n: usize) -> Self {
        self.config.max_measure = Some(n.max(100));
        self
    }

    /// Keep only the brightest N stars in the returned result.
    pub fn with_max_stars(mut self, n: usize) -> Self {
        let n = n.max(1);
        // Clamp to max_measure if set
        self.config.max_stars = match self.config.max_measure {
            Some(mm) => n.min(mm),
            None => n,
        };
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

    /// Set the R² threshold for trail detection.
    /// Images with Rayleigh R² above this are flagged as possibly trailed.
    /// Default: 0.5. Lower values are more aggressive (more false positives).
    pub fn with_trail_threshold(mut self, threshold: f32) -> Self {
        self.config.trail_r_squared_threshold = threshold.clamp(0.0, 1.0);
        self
    }

    /// Enable iterative source-masked background estimation.
    /// After initial detection, source pixels are masked and background is re-estimated.
    /// Default: 1 (no iteration). Set to 2-3 for images with many bright sources.
    /// Requires `with_background_mesh` to be set (no-op with global background).
    pub fn with_iterative_background(mut self, iterations: usize) -> Self {
        self.config.iterative_background = iterations.max(1);
        self
    }

    /// Enable MRS (Multiresolution Support) wavelet noise estimation.
    /// Uses à trous B3-spline wavelet to isolate noise from nebulosity/gradients.
    /// `layers`: number of wavelet layers (1 is sufficient).
    /// Default: 0 (sigma-clipped MAD).
    pub fn with_mrs_noise(mut self, layers: usize) -> Self {
        self.config.noise_layers = layers;
        self
    }

    /// Fix Moffat beta to a constant value during PSF fitting.
    /// Removes the beta/axis-ratio tradeoff, improving FWHM stability.
    /// Typical value: 4.0 (well-corrected optics).
    /// Default: None (free beta).
    pub fn with_moffat_beta(mut self, beta: f32) -> Self {
        self.config.moffat_beta = Some(beta);
        self
    }

    /// Reject detected stars with eccentricity above this threshold before measurement.
    /// Pre-filters elongated candidates (optical ghosts, diffraction spikes, edge effects).
    /// Typical value: 0.6. Default: None (no filter).
    pub fn with_max_distortion(mut self, ecc: f32) -> Self {
        self.config.max_distortion = Some(ecc.clamp(0.0, 1.0));
        self
    }

    /// Use Moffat PSF fitting instead of Gaussian for more accurate wing modeling.
    /// Falls back to Gaussian on non-convergence. Reports per-star β values.
    /// This is the default since v0.7.0.
    pub fn with_moffat_fit(mut self) -> Self {
        self.config.use_moffat_fit = true;
        self
    }

    /// Disable Moffat PSF fitting; use Gaussian fit instead.
    pub fn without_moffat_fit(mut self) -> Self {
        self.config.use_moffat_fit = false;
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

        let green_mask = if self.config.apply_debayer
            && meta.bayer_pattern != BayerPattern::None
            && channels == 1
        {
            data = debayer::interpolate_green_f32(&data, width, height, meta.bayer_pattern);
            Some(build_green_mask(width, height, meta.bayer_pattern))
        } else {
            None
        };

        self.run_analysis(&data, width, height, channels, green_mask.as_deref())
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

        // Star detection parameters
        let det_params = DetectionParams {
            detection_sigma: self.config.detection_sigma,
            min_star_area: self.config.min_star_area,
            max_star_area: self.config.max_star_area,
            saturation_limit: self.config.saturation_fraction * 65535.0,
        };

        // Iterative background estimation + detection loop
        let n_iterations = if self.config.background_mesh_size.is_some() {
            self.config.iterative_background
        } else {
            1 // no iteration without mesh background
        };

        let mut bg_result = if let Some(cell_size) = self.config.background_mesh_size {
            background::estimate_background_mesh(&lum, width, height, cell_size)
        } else {
            background::estimate_background(&lum, width, height)
        };

        // MRS wavelet noise override
        if self.config.noise_layers > 0 {
            bg_result.noise = background::estimate_noise_mrs(
                &lum, width, height, self.config.noise_layers,
            ).max(0.001);
        }

        let mut detected;
        let final_fwhm;

        // First detection pass (always runs)
        {
            let bg_map_ref = bg_result.background_map.as_deref();
            let noise_map_ref = bg_result.noise_map.as_deref();
            let initial_fwhm = 3.0_f32;
            let pass1 = detection::detect_stars(
                &lum, width, height,
                bg_result.background, bg_result.noise,
                bg_map_ref, noise_map_ref, &det_params, initial_fwhm,
            );

            let measured_kernel_fwhm = estimate_fwhm_from_stars(
                &lum, width, height, &pass1, bg_result.background, bg_map_ref,
            );

            // Cap second-pass kernel to 2× initial to prevent runaway cascade
            let capped_kernel_fwhm = measured_kernel_fwhm.min(initial_fwhm * 2.0);

            if capped_kernel_fwhm > 0.0
                && ((capped_kernel_fwhm - initial_fwhm) / initial_fwhm).abs() > 0.30
            {
                detected = detection::detect_stars(
                    &lum, width, height,
                    bg_result.background, bg_result.noise,
                    bg_map_ref, noise_map_ref, &det_params, capped_kernel_fwhm,
                );
                final_fwhm = capped_kernel_fwhm;
            } else {
                detected = pass1;
                final_fwhm = initial_fwhm;
            }
        }

        // Iterative source-masked background re-estimation (iterations 2..n)
        if let Some(cell_size) = self.config.background_mesh_size {
            for _ in 1..n_iterations {
                if detected.is_empty() {
                    break;
                }

                // Build source mask: circle of r = 2.5 × FWHM around each star
                let mask_radius = (2.5 * final_fwhm).ceil() as i32;
                let mask_r_sq = (mask_radius * mask_radius) as f32;
                let mut source_mask = vec![false; width * height];
                for star in &detected {
                    let cx = star.x.round() as i32;
                    let cy = star.y.round() as i32;
                    for dy in -mask_radius..=mask_radius {
                        let py = cy + dy;
                        if py < 0 || py >= height as i32 { continue; }
                        for dx in -mask_radius..=mask_radius {
                            let px = cx + dx;
                            if px < 0 || px >= width as i32 { continue; }
                            if (dx * dx + dy * dy) as f32 <= mask_r_sq {
                                source_mask[py as usize * width + px as usize] = true;
                            }
                        }
                    }
                }

                // Re-estimate background excluding masked pixels
                bg_result = background::estimate_background_mesh_masked(
                    &lum, width, height, cell_size, &source_mask,
                );

                // Re-apply MRS noise override
                if self.config.noise_layers > 0 {
                    bg_result.noise = background::estimate_noise_mrs(
                        &lum, width, height, self.config.noise_layers,
                    ).max(0.001);
                }

                // Re-detect with updated background and noise map
                let bg_map_ref = bg_result.background_map.as_deref();
                let noise_map_ref = bg_result.noise_map.as_deref();
                detected = detection::detect_stars(
                    &lum, width, height,
                    bg_result.background, bg_result.noise,
                    bg_map_ref, noise_map_ref, &det_params, final_fwhm,
                );
            }
        }

        let bg_map_ref = bg_result.background_map.as_deref();

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

        let snr_db = snr::compute_snr_db(&lum, bg_result.noise);
        let snr_weight = snr::compute_snr_weight(&lum, bg_result.background, bg_result.noise);

        // Helper for zero-star results
        let make_zero_result = |stars_detected: usize| {
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
                snr_db,
                snr_weight,
                psf_signal: 0.0,
                trail_r_squared,
                possibly_trailed,
                measured_fwhm_kernel: final_fwhm,
                median_beta: None,
            })
        };

        if detected.is_empty() {
            return make_zero_result(0);
        }

        // Optional early cap: limit stars measured for speed
        // (detected are already sorted by flux descending)
        if let Some(max_m) = self.config.max_measure {
            detected.truncate(max_m);
        }

        // Eccentricity pre-filter: reject elongated candidates before measurement
        if let Some(max_ecc) = self.config.max_distortion {
            detected.retain(|s| s.eccentricity <= max_ecc);
        }

        // Measure PSF metrics on all (possibly capped) detected stars
        let mut measured = metrics::measure_stars(
            &lum,
            width,
            height,
            &detected,
            bg_result.background,
            bg_map_ref,
            self.config.use_gaussian_fit,
            self.config.use_moffat_fit,
            green_mask,
            self.config.moffat_beta.map(|b| b as f64),
        );

        if measured.is_empty() {
            return make_zero_result(0);
        }

        // True count of stars with valid PSF measurements
        let stars_detected = measured.len();

        // Compute statistics from ALL measured stars (no filtering, no cap)
        // Use sigma-clipped median (2-iter, 3σ MAD) for FWHM, ecc, HFR to reject outliers
        let fwhm_vals: Vec<f32> = measured.iter().map(|s| s.fwhm).collect();
        let median_fwhm = sigma_clipped_median(&fwhm_vals);

        // Per-star SNR on ALL measured stars
        snr::compute_star_snr(&lum, width, height, &mut measured, median_fwhm);

        // Summary statistics from full population
        let ecc_vals: Vec<f32> = measured.iter().map(|s| s.eccentricity).collect();
        let mut snr_vals: Vec<f32> = measured.iter().map(|s| s.snr).collect();
        let hfr_vals: Vec<f32> = measured.iter().map(|s| s.hfr).collect();

        let median_eccentricity = sigma_clipped_median(&ecc_vals);
        let median_snr = find_median(&mut snr_vals); // SNR left unclipped
        let median_hfr = sigma_clipped_median(&hfr_vals);

        let psf_signal = snr::compute_psf_signal(&measured, bg_result.noise);

        // Compute median beta if Moffat fitting was used
        let median_beta = if self.config.use_moffat_fit {
            let mut beta_vals: Vec<f32> = measured.iter().filter_map(|s| s.beta).collect();
            if beta_vals.is_empty() {
                None
            } else {
                Some(find_median(&mut beta_vals))
            }
        } else {
            None
        };

        // Late cap: truncate to max_stars AFTER all statistics are computed
        measured.truncate(self.config.max_stars);

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
                beta: m.beta,
                fit_method: m.fit_method,
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
            measured_fwhm_kernel: final_fwhm,
            median_beta,
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
        let stamp_radius = 12_usize; // enough for FWHM up to ~10px
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

    let green_mask = if apply_debayer
        && meta.bayer_pattern != BayerPattern::None
        && channels == 1
    {
        data = debayer::interpolate_green_f32(&data, width, height, meta.bayer_pattern);
        Some(build_green_mask(width, height, meta.bayer_pattern))
    } else {
        None
    };

    let lum = if channels == 3 {
        extract_luminance(&data, width, height)
    } else {
        data[..width * height].to_vec()
    };

    (lum, width, height, channels, green_mask)
}
