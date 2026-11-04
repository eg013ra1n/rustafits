//! Robust location and scale statistics (Athenaeum's math reference §2.4–2.5, spec §5.1): median / MAD / average
//! deviation / biweight midvariance evaluated two-sided about the median, the 1/16 stratified pixel sample every caller
//! shares, and the noise-scale factors of the classic SNR weight. Pure functions over `f32` samples.
//!
//! Moved verbatim from Athenaeum's `integration::stats` (at 0b6fb183) for [`crate::star_measure`]; Athenaeum's
//! integration engine keeps its per-frame `(scale, offset)` normalization pairs and re-exports these. Their numbers are
//! Athenaeum's stacking numbers: no arithmetic here changes without its `MEASURE_VERSION` policy.

// Athenaeum's code as copied (port plan G-3): the lint is allowed, not "fixed".
#![allow(clippy::neg_cmp_op_on_partial_ord)]

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// σ-consistency factor for the median absolute deviation.
pub const MAD_TO_SIGMA: f32 = 1.4826;
/// σ-consistency factor for the average absolute deviation (√(π/2)).
pub const AVGDEV_TO_SIGMA: f32 = 1.2533;
/// σ-consistency factor applied to √BWMV.
pub const BWMV_TO_SIGMA: f32 = 0.991;
/// Samples outside `(CLIP_LO, CLIP_HI)` are excluded from scale estimates:
/// zeros mark missing coverage, the top value marks saturation.
pub const CLIP_LO: f32 = 1.0 / 65535.0;
pub const CLIP_HI: f32 = 1.0 - 1.0 / 65535.0;
/// The tighter clip the noise-scale factors use (§2.5).
pub const NOISE_CLIP_LO: f32 = 2.0 / 65535.0;
pub const NOISE_CLIP_HI: f32 = 1.0 - 2.0 / 65535.0;
/// Stride of the stratified sample: one pixel in `STRIDE²` = 1/16.
pub const SAMPLE_STRIDE: usize = 4;

/// Median with the even-`n` convention "mean of the central two". Reorders
/// `values`; the caller has removed NaNs. NaN for an empty slice.
pub fn median_in_place(values: &mut [f32]) -> f32 {
    let n = values.len();
    if n == 0 {
        return f32::NAN;
    }
    let mid = n / 2;
    let (_, hi, _) = values.select_nth_unstable_by(mid, |a, b| a.total_cmp(b));
    let hi = *hi;
    if n % 2 == 1 {
        hi
    } else {
        let lo = values[..mid]
            .iter()
            .copied()
            .fold(f32::NEG_INFINITY, f32::max);
        0.5 * (lo + hi)
    }
}

pub fn median_of(values: &[f32]) -> f32 {
    let mut v = values.to_vec();
    median_in_place(&mut v)
}

/// Raw median absolute deviation about `m` (no consistency factor).
pub fn mad_about(values: &[f32], m: f32) -> f32 {
    let mut d: Vec<f32> = values.iter().map(|&x| (x - m).abs()).collect();
    median_in_place(&mut d)
}

/// Mean absolute deviation about `m`.
pub fn avg_dev_about(values: &[f32], m: f32) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    let s: f64 = values.iter().map(|&x| (x - m).abs() as f64).sum();
    (s / values.len() as f64) as f32
}

/// σ-consistent biweight midvariance scale about `m` (Wilcox 2012
/// §3.12.1): `u = (x − m)/(9·MAD)`, `BWMV = n·Σ_{|u|<1}(x − m)²(1 − u²)⁴ /
/// [Σ_{|u|<1}(1 − u²)(1 − 5u²)]²`, scale `= √BWMV · 0.991`. `mad` is the
/// raw MAD about `m`. Falls back to the MAD scale when no sample lies in
/// the biweight window or the denominator is not positive.
pub fn bwmv_scale_about(values: &[f32], m: f32, mad: f32) -> f32 {
    let n = values.len();
    if n == 0 {
        return f32::NAN;
    }
    if !(mad > 0.0) {
        return 0.0;
    }
    let w = 9.0 * mad as f64;
    let (mut num, mut den) = (0.0f64, 0.0f64);
    for &x in values {
        let d = (x - m) as f64;
        let u2 = (d / w) * (d / w);
        if u2 < 1.0 {
            let a = 1.0 - u2;
            num += d * d * a.powi(4);
            den += a * (1.0 - 5.0 * u2);
        }
    }
    if den <= 0.0 || num <= 0.0 {
        return mad * MAD_TO_SIGMA;
    }
    let bwmv = n as f64 * num / (den * den);
    (bwmv.sqrt() as f32) * BWMV_TO_SIGMA
}

fn stddev_about_mean(v: &[f32]) -> f64 {
    let n = v.len();
    if n < 2 {
        return 0.0;
    }
    let mean = v.iter().map(|&x| x as f64).sum::<f64>() / n as f64;
    let ss = v
        .iter()
        .map(|&x| {
            let d = x as f64 - mean;
            d * d
        })
        .sum::<f64>();
    (ss / (n - 1) as f64).sqrt()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "camelCase"))]
pub enum ScaleEstimator {
    #[default]
    Bwmv,
    Mad,
    AvgDev,
}

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "camelCase"))]
pub struct LocationScale {
    pub location: f32,
    pub scale: f32,
}

/// Keep the finite samples strictly inside `(lo, hi)`.
pub fn clip_sample(values: &[f32], lo: f32, hi: f32) -> Vec<f32> {
    values
        .iter()
        .copied()
        .filter(|v| v.is_finite() && *v > lo && *v < hi)
        .collect()
}

fn side_scale(side: &[f32], m: f32, est: ScaleEstimator) -> Option<f32> {
    if side.len() < 2 {
        return None;
    }
    Some(match est {
        ScaleEstimator::Mad => mad_about(side, m) * MAD_TO_SIGMA,
        ScaleEstimator::AvgDev => avg_dev_about(side, m) * AVGDEV_TO_SIGMA,
        ScaleEstimator::Bwmv => bwmv_scale_about(side, m, mad_about(side, m)),
    })
}

/// Location (median) and two-sided scale of an already clipped sample: the
/// estimator runs separately on `x ≤ m` and `x > m` and the two are
/// averaged; an empty side defers to the other. `None` below 4 samples.
pub fn location_scale(values: &[f32], est: ScaleEstimator) -> Option<LocationScale> {
    if values.len() < 4 {
        return None;
    }
    let m = median_of(values);
    let (low, high): (Vec<f32>, Vec<f32>) = values.iter().partition(|&&x| x <= m);
    let scale = match (side_scale(&low, m, est), side_scale(&high, m, est)) {
        (Some(a), Some(b)) => 0.5 * (a + b),
        (Some(a), None) | (None, Some(a)) => a,
        (None, None) => return None,
    };
    Some(LocationScale { location: m, scale })
}

/// Visit the 1/16 stratified pixel indices: rows `0, 4, 8, …`; in row `y`
/// the columns start at `(y / 4) % 4` and step by 4, so the four row phases
/// cover all four column phases.
pub fn for_each_stratified(width: usize, height: usize, mut f: impl FnMut(usize)) {
    let mut y = 0;
    while y < height {
        let mut x = (y / SAMPLE_STRIDE) % SAMPLE_STRIDE;
        while x < width {
            f(y * width + x);
            x += SAMPLE_STRIDE;
        }
        y += SAMPLE_STRIDE;
    }
}

/// The finite values at the stratified indices of a plane.
pub fn stratified_sample(data: &[f32], width: usize, height: usize) -> Vec<f32> {
    let mut out = Vec::with_capacity(data.len() / (SAMPLE_STRIDE * SAMPLE_STRIDE) + 1);
    for_each_stratified(width, height, |i| {
        let v = data[i];
        if v.is_finite() {
            out.push(v);
        }
    });
    out
}

/// Noise scaling factors (math reference §2.5): `(σ_low, σ_high)` — the
/// standard deviation of each side of the clipped sample after two
/// restrictions to `[c − 4σ, c]` / `[c, c + 4σ]` about the median `c`.
/// `None` below 8 samples or when a side empties.
pub fn noise_scale_factors(values: &[f32]) -> Option<(f32, f32)> {
    noise_scale_factors_within(values, NOISE_CLIP_LO, NOISE_CLIP_HI)
}

/// [`noise_scale_factors`] with the clip range `(clip_lo, clip_hi)` given
/// instead of the unit domain's `[NOISE_CLIP_LO, NOISE_CLIP_HI]` — the same
/// arithmetic, so the unit-domain call is bit-identical; a raw-ADU plane
/// passes its own range (`star_measure::sky_level`).
pub fn noise_scale_factors_within(
    values: &[f32],
    clip_lo: f32,
    clip_hi: f32,
) -> Option<(f32, f32)> {
    let v = clip_sample(values, clip_lo, clip_hi);
    if v.len() < 8 {
        return None;
    }
    let c = median_of(&v);
    let mut low: Vec<f32> = v.iter().copied().filter(|&x| x <= c).collect();
    let mut high: Vec<f32> = v.iter().copied().filter(|&x| x >= c).collect();
    for _ in 0..2 {
        let lo = (c as f64 - 4.0 * stddev_about_mean(&low)).max(clip_lo as f64) as f32;
        low.retain(|&x| x >= lo);
        let hi = (c as f64 + 4.0 * stddev_about_mean(&high)).min(clip_hi as f64) as f32;
        high.retain(|&x| x <= hi);
    }
    if low.len() < 2 || high.len() < 2 {
        return None;
    }
    Some((
        stddev_about_mean(&low) as f32,
        stddev_about_mean(&high) as f32,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::add_noise;

    fn gaussian(n: usize, sigma: f32, seed: u64) -> Vec<f32> {
        let mut v = vec![0.5f32; n];
        add_noise(&mut v, sigma, seed);
        v
    }
    fn close(a: f32, b: f32, rel: f32) -> bool {
        (a - b).abs() <= rel * b.abs()
    }

    #[test]
    fn median_conventions() {
        assert_eq!(median_of(&[5.0, 1.0, 4.0, 2.0, 3.0]), 3.0);
        assert_eq!(median_of(&[4.0, 1.0, 3.0, 2.0]), 2.5);
        assert_eq!(median_of(&[7.0]), 7.0);
        assert!(median_of(&[]).is_nan());
    }

    #[test]
    fn mad_and_avg_dev_by_hand() {
        let v = [1.0, 2.0, 3.0, 4.0, 100.0];
        let m = median_of(&v);
        assert_eq!(m, 3.0);
        assert_eq!(mad_about(&v, m), 1.0);
        assert!((avg_dev_about(&v, m) - 101.0 / 5.0).abs() < 1e-5);
    }

    #[test]
    fn estimators_are_sigma_consistent_on_gaussian_noise() {
        let v = gaussian(20_000, 0.01, 1);
        let m = median_of(&v);
        let mad = mad_about(&v, m);
        assert!(close(mad * MAD_TO_SIGMA, 0.01, 0.03));
        assert!(close(avg_dev_about(&v, m) * AVGDEV_TO_SIGMA, 0.01, 0.03));
        assert!(close(bwmv_scale_about(&v, m, mad), 0.01, 0.03));
        for est in [
            ScaleEstimator::Bwmv,
            ScaleEstimator::Mad,
            ScaleEstimator::AvgDev,
        ] {
            let ls = location_scale(&v, est).unwrap();
            assert!(
                close(ls.location, 0.5, 0.001),
                "{est:?} location {}",
                ls.location
            );
            assert!(close(ls.scale, 0.01, 0.03), "{est:?} scale {}", ls.scale);
        }
        assert!(location_scale(&[0.1, 0.2, 0.3], ScaleEstimator::Mad).is_none());
    }

    #[test]
    fn robust_scales_ignore_symmetric_outliers() {
        let mut v = gaussian(20_000, 0.01, 2);
        for i in 0..200 {
            v[i * 100] = 0.5 + 0.2; // 1 % at +20σ
            v[i * 100 + 50] = 0.5 - 0.2; // 1 % at −20σ
        }
        let m = median_of(&v);
        let mad = mad_about(&v, m);
        assert!(close(mad * MAD_TO_SIGMA, 0.01, 0.05));
        assert!(close(bwmv_scale_about(&v, m, mad), 0.01, 0.05));
        let sd = stddev_about_mean(&v) as f32;
        assert!(
            sd > 0.025,
            "plain stddev must be blown up by the outliers: {sd}"
        );
        assert!(close(
            location_scale(&v, ScaleEstimator::Bwmv).unwrap().scale,
            0.01,
            0.05
        ));
        assert!(close(
            location_scale(&v, ScaleEstimator::Mad).unwrap().scale,
            0.01,
            0.05
        ));
    }

    #[test]
    fn two_sided_equals_one_sided_on_symmetric_data() {
        let v = gaussian(20_000, 0.02, 3);
        let m = median_of(&v);
        let one = mad_about(&v, m) * MAD_TO_SIGMA;
        let two = location_scale(&v, ScaleEstimator::Mad).unwrap().scale;
        assert!(close(two, one, 0.03));
    }

    #[test]
    fn clip_drops_zeros_saturation_and_nan() {
        let v = [0.0, 0.5, 1.0, f32::NAN, 0.25];
        assert_eq!(clip_sample(&v, CLIP_LO, CLIP_HI), vec![0.5, 0.25]);
    }

    #[test]
    fn stratified_sample_counts() {
        let d = vec![1.0f32; 16 * 16];
        assert_eq!(stratified_sample(&d, 16, 16).len(), 16);
        let d = vec![1.0f32; 17 * 9];
        assert_eq!(stratified_sample(&d, 17, 9).len(), 13);
        let mut d = vec![1.0f32; 8 * 8];
        d[0] = f32::NAN;
        assert_eq!(stratified_sample(&d, 8, 8).len(), 3);
        let mut seen = Vec::new();
        for_each_stratified(8, 8, |i| seen.push(i));
        assert_eq!(seen, vec![0, 4, 33, 37]);
    }

    #[test]
    fn noise_scale_factors_are_symmetric_on_gaussian_noise() {
        let v = gaussian(50_000, 0.01, 4);
        let (lo, hi) = noise_scale_factors(&v).unwrap();
        assert!(lo > 0.0046 && lo < 0.0062, "{lo}");
        assert!(hi > 0.0046 && hi < 0.0062, "{hi}");
        assert!((lo - hi).abs() < 0.0005);
        assert!(noise_scale_factors(&[0.5; 4]).is_none());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn scale_estimator_serde_names() {
        assert_eq!(
            serde_json::to_string(&ScaleEstimator::AvgDev).unwrap(),
            "\"avgDev\""
        );
        assert_eq!(
            serde_json::to_string(&ScaleEstimator::Bwmv).unwrap(),
            "\"bwmv\""
        );
        assert_eq!(ScaleEstimator::default(), ScaleEstimator::Bwmv);
    }
}
