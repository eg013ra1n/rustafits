# Athenaeum Integration: Advisory Trail Detection (v0.4.5)

## What Changed

Trail detection in the analysis pipeline is now **advisory**. Previously, images
flagged as trailed by the Rayleigh test returned a zero-result (empty star list,
zeroed metrics). Now the pipeline always computes full metrics and reports trail
status via two new fields on `AnalysisResult`.

This is a non-breaking additive change — existing code continues to work. The only
behavioral difference is that previously-rejected trailed images now return full
metrics instead of zeros.

## New Fields on `AnalysisResult`

```rust
pub trail_r_squared: f32,
pub possibly_trailed: bool,
```

### `trail_r_squared`

The Rayleigh R² statistic measuring directional coherence of star position angles.
Computed on detected stars before PSF measurement, using stamp-based second-order
moments on doubled angles (2*theta).

| Value | Meaning |
|-------|---------|
| 0.0 | Perfectly uniform angles (no trail) or fewer than 5 stars detected |
| 0.01-0.15 | Typical for undersampled round stars (grid-induced coherence) |
| 0.15-0.40 | Typical for oversampled stars with field-angle effects (coma, curvature) |
| 0.50-0.70 | Suspicious — possible mild trailing |
| 0.70-1.00 | Strong trail (tracking drift, wind shake, satellite) |

### `possibly_trailed`

Boolean flag set to `true` when either detection path fires:

- **Path A**: `trail_r_squared > threshold` AND Rayleigh p < 0.01
- **Path B**: median detection-stage eccentricity > 0.6 AND Rayleigh p < 0.05

The threshold defaults to 0.5 and is configurable.

## New Builder Method

```rust
pub fn with_trail_threshold(mut self, threshold: f32) -> Self
```

Sets the R² threshold used by Path A. Clamped to [0.0, 1.0]. Default: 0.5.

The raw `trail_r_squared` is always computed regardless of this setting, so
Athenaeum can apply its own threshold logic downstream.

## Integration Recommendations

### 1. Display R² in the subframe inspector

Show `trail_r_squared` alongside existing metrics (FWHM, eccentricity, SNR).
It's a continuous quality indicator independent of the boolean flag — useful for
users evaluating their dataset even when no frames are flagged.

Suggested column header: **R²** or **Trail R²**.
Format: `{:.3}` (e.g., 0.142, 0.703).

### 2. Use `possibly_trailed` as a warning, not auto-reject

Display a warning icon or badge on flagged frames. Do not auto-reject — the
whole point of making this advisory is that the 0.5 threshold produces false
positives on some optical systems (fast Newtonians with coma, wide-field
refractors with curvature).

Let users review flagged frames and decide. If Athenaeum has a "reject" workflow,
make trail-flagged frames candidates for rejection but require user confirmation.

### 3. Consider a user-configurable threshold in preferences

Expose the R² threshold as a preference (e.g., "Trail sensitivity"). Map it to
`with_trail_threshold()`. Users with optical systems that produce high baseline
R² (coma, curvature) can raise the threshold to reduce false positives.

Suggested default: 0.5 (matches rustafits default).
Suggested range: 0.3 to 0.8 (slider or numeric input).

### 4. Sorting and filtering

Allow sorting the subframe list by `trail_r_squared` to quickly find the worst
frames. Allow filtering by `possibly_trailed == true` to isolate flagged frames.

### 5. Previously-zero results now have data

Images that previously returned `stars.len() == 0` and zeroed medians due to trail
rejection will now return full metrics. If Athenaeum was treating zero-star results
as "bad frame" signals, verify that logic still works correctly — zero stars now
only means genuinely no detectable stars (very cloudy, completely trailed with no
detectable knots, or blank frame), not "trailed but detectable stars present."

The `possibly_trailed` flag is the new way to identify trailed frames.

## Migration Checklist

- [ ] Update rustafits dependency to 0.4.5
- [ ] Add `trail_r_squared` column to subframe table/inspector
- [ ] Add `possibly_trailed` warning indicator (icon/badge) on flagged frames
- [ ] Review any logic that interprets `stars.len() == 0` as "bad frame" — trailed
      frames now return stars; use `possibly_trailed` instead
- [ ] (Optional) Add trail threshold preference mapped to `with_trail_threshold()`
- [ ] (Optional) Add sort-by-R² and filter-by-trailed to subframe list

## Code Example

```rust
use rustafits::ImageAnalyzer;

// Default threshold (0.5)
let analyzer = ImageAnalyzer::new();
let result = analyzer.analyze("frame.fits")?;

// Access new fields
let r_sq = result.trail_r_squared;       // f32, always computed
let trailed = result.possibly_trailed;   // bool, based on threshold

// Full metrics are always available, even for trailed frames
let fwhm = result.median_fwhm;           // non-zero if stars detected
let stars = result.stars.len();           // non-zero if stars detected

// Custom threshold
let strict = ImageAnalyzer::new()
    .with_trail_threshold(0.4)            // more aggressive
    .analyze("frame.fits")?;

// Or ignore the flag entirely and use raw R² with your own logic
if result.trail_r_squared > 0.6 {
    // custom threshold
}
```
