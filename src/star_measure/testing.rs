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
