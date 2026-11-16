//! `StarFit::test_at`, the fixture Athenaeum's tests and this module's own build fits with (moved from Athenaeum's
//! `test_support` at 0b6fb183). Compiled for this crate's tests and under the `test-support` feature; not
//! semver-covered.

/// A plausible accepted fit at `(x, y)` with `amplitude` (Plan A Task 4): the fixture the `StarFit` literal sites can
/// call instead of spelling every field (ruling F36). Every other field is a round, well-sampled star — FWHM 3 px,
/// FWTM 6 px, β 4, a small residual, a positive signal over its FWTM-ellipse area — with no flags; `peak5` is
/// unmeasured (NaN), as on a fit decoded from a `fits` artifact.
impl crate::star_measure::StarFit {
    #[doc(hidden)]
    pub fn test_at(x: f64, y: f64, amplitude: f64) -> Self {
        let fwtm = 6.0;
        crate::star_measure::StarFit {
            x,
            y,
            background: 100.0,
            amplitude,
            fwhm_x: 3.0,
            fwhm_y: 3.0,
            fwtm_x: fwtm,
            fwtm_y: fwtm,
            theta: 0.0,
            beta: 4.0,
            residual: 0.02,
            signal: 8.0 * amplitude,
            area: std::f64::consts::PI * 0.25 * fwtm * fwtm,
            peak5: f64::NAN,
            sky_ring: f64::NAN,
            sky_ring_sigma: f64::NAN,
            sky_ring_n: 0,
            ap_sum: f64::NAN,
            snr: f64::NAN,
            hfr: f64::NAN,
            flags: crate::star_measure::StarFlags::default(),
        }
    }
}

/// The list entry's rows for a caller's own tests (`ListFit` and `ListMeasurement` are `#[non_exhaustive]`, so a caller
/// outside this crate cannot spell them): build one, then set the public fields a case needs.
impl crate::star_measure::analysis_mode::ListFit {
    /// A kept fit at `(x, y)` with `flux`: accepted, kept, unflagged, centred on its position, SNR 50, model flux = `flux`
    /// (`c = 1`).
    #[doc(hidden)]
    pub fn test_kept(x: f64, y: f64, flux: f64) -> Self {
        crate::star_measure::analysis_mode::ListFit {
            x,
            y,
            model_flux: flux,
            flux,
            snr: 50.0,
            flags: crate::star_measure::StarFlags::default(),
            accepted: true,
            kept: true,
            centre_offset_px: 0.0,
        }
    }

    /// A position that was refused: NaN numbers, both flags false.
    #[doc(hidden)]
    pub fn test_refused() -> Self {
        crate::star_measure::list::refused_fit()
    }
}

impl crate::star_measure::analysis_mode::ListMeasurement {
    /// A list of `fits` with aperture correction `aperture_correction` over `aperture_stars` stars: β 4, ring width 3 px,
    /// no standard error, no refusals.
    #[doc(hidden)]
    pub fn test_new(
        fits: Vec<crate::star_measure::analysis_mode::ListFit>,
        aperture_correction: Option<f64>,
        aperture_stars: usize,
    ) -> Self {
        crate::star_measure::analysis_mode::ListMeasurement {
            fits,
            aperture_correction,
            aperture_stars,
            aperture_se: None,
            beta: 4.0,
            ring_fwhm: Some(3.0),
            refused: crate::star_measure::RefusalCounts::default(),
            refusal: None,
        }
    }
}
