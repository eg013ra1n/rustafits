# Trail Detection (Rayleigh Test)

Image-level detection of satellite trails, tracking errors, and wind shake using
circular statistics on star position angles. Reports an advisory flag and raw R²
statistic — the caller decides whether to reject.

> **Changed in v0.4.5:** The Rayleigh test is now **advisory**. It no longer
> rejects images (returns zero-result). Instead it always computes full metrics
> and exposes `trail_r_squared` and `possibly_trailed` on `AnalysisResult`.
> The R² threshold is configurable via `with_trail_threshold()` (default 0.5).

## The Problem

Trailed images (from tracking drift, wind, or cable snag) produce elongated star
profiles where all stars are stretched in the same direction. Traditional per-star
eccentricity thresholds can miss mild trailing where individual stars look only
slightly elongated, but the ensemble clearly has a preferred direction.

## Approach: Rayleigh Test on 2-Theta

The [Rayleigh test](https://en.wikipedia.org/wiki/Rayleigh_test) detects non-uniformity
in a circular distribution. We use it on the doubled position angles (2 * theta)
because theta is an axial quantity with period pi (an angle of 0 and pi represent
the same orientation).

### Statistic

Given n detected stars with position angles theta_1 ... theta_n:

```
sum_cos = sum( cos(2 * theta_i) )
sum_sin = sum( sin(2 * theta_i) )

R^2 = (sum_cos^2 + sum_sin^2) / n^2

p = exp(-n * R^2)
```

R^2 is the squared mean resultant length — it measures how concentrated the angles
are. R^2 = 0 means perfectly uniform (random angles), R^2 = 1 means all angles
identical (strong trail).

The p-value gives the probability of observing R^2 this large or larger under the
null hypothesis of uniform random angles.

### Minimum Star Count

The test requires at least 5 detected stars. Below this, there are too few samples
for meaningful circular statistics and the test is skipped (`trail_r_squared = 0.0`,
`possibly_trailed = false`).

---

## Dual-Path Detection

A single threshold does not work for all images because undersampled stars
(FWHM < 3 px) produce grid-induced theta coherence — the pixel grid geometry
biases the moment tensor, making round stars appear to have a preferred angle.

This creates a regime-dependent false positive problem:

| Regime | Stars look like | Theta behavior | Risk |
|--------|----------------|----------------|------|
| Undersampled, round | Small, boxy blobs | Grid-coherent (R^2 ~ 0.15) | False positive |
| Undersampled, trailed | Elongated blobs | Trail-coherent (R^2 > 0.5) | Correctly detected |
| Oversampled, round | Smooth, circular | Random (R^2 ~ 1/n) | No risk |
| Oversampled, trailed | Smooth, elongated knots | Trail-coherent (R^2 > 0.5) | Missed by ecc gate |

The solution is two detection paths:

### Path A — Strong Directional Coherence (R^2 > threshold)

```
Flag if:  R^2 > threshold  AND  p < 0.01
```

The threshold defaults to 0.5 and is configurable via `with_trail_threshold()`.

Real satellite/tracking trails produce R^2 > 0.7 (strong coherence). Grid-induced
coherence on non-trailed stars varies by regime:
- Undersampled round stars (FWHM ~ 2 px): R^2 ~ 0.15
- Oversampled stars with field-angle effects (FWHM ~ 5 px): R^2 up to 0.40

The default 0.5 threshold sits above both regimes while still catching real trails.

This path catches **trails** where individual CCL knots have low eccentricity
(~0.3-0.5) but consistent theta. Without this path, the eccentricity gate
(Path B) would never fire and such trails would slip through.

### Path B — Eccentricity-Gated Rayleigh (median ecc > 0.6)

```
Flag if:  median_eccentricity > 0.6  AND  p < 0.05
```

For undersampled stars, stamp-based moments are noisy and theta has grid bias.
This path only runs the Rayleigh test when blobs are **genuinely elongated**
(median eccentricity > 0.6), which exceeds the undersampled baseline (~0.3-0.4).

This is the original gate — it catches **undersampled trails** where eccentricity
is clearly above the noise floor.

### Combined Decision

```
possibly_trailed = (R^2 > threshold AND p < 0.01)    -- Path A
                OR (median_ecc > 0.6 AND p < 0.05)   -- Path B
```

Both paths compute the same Rayleigh statistic from the same star list; they
differ only in their activation criteria and p-value threshold.

---

## API Surface

The Rayleigh test result is exposed via two fields on `AnalysisResult`:

| Field | Type | Description |
|-------|------|-------------|
| `trail_r_squared` | `f32` | Raw R² statistic. 0.0 = uniform, 1.0 = perfectly aligned. |
| `possibly_trailed` | `bool` | True if either detection path fired. |

The R² threshold for Path A is configurable:

```rust
let analyzer = ImageAnalyzer::new()
    .with_trail_threshold(0.4);  // more aggressive (more false positives)
```

The `trail_r_squared` value is always computed regardless of the threshold setting,
so callers can implement their own logic:

```rust
let result = ImageAnalyzer::new().analyze("frame.fits")?;

// Use the built-in flag
if result.possibly_trailed {
    println!("Warning: possible trail (R²={:.3})", result.trail_r_squared);
}

// Or apply a custom threshold
if result.trail_r_squared > 0.6 {
    // reject with stricter threshold
}
```

---

## How It Works: Scenarios

### Good image, undersampled stars (FWHM ~ 2 px)

```
Detected: 200 stars
Median eccentricity: ~0.48
R^2: ~0.15 (grid-induced)
p: ~0 (statistically significant but weak)

Path A: R^2 = 0.15 < 0.5 --> does not fire
Path B: median_ecc = 0.48 < 0.6 --> does not fire

Result: possibly_trailed = false, trail_r_squared = 0.15
        Full metrics computed normally.
```

### Good image, oversampled stars (FWHM ~ 5 px)

```
Detected: 200 stars
Median eccentricity: ~0.55
R^2: ~0.37 (field-angle effects: coma, curvature)
p: ~0

Path A: R^2 = 0.37 < 0.5 --> does not fire
Path B: median_ecc = 0.55 < 0.6 --> does not fire

Result: possibly_trailed = false, trail_r_squared = 0.37
        Full metrics computed normally.
```

### Trailed image, oversampled (tracking drift)

```
Detected: 100 stars
Median eccentricity: ~0.45
R^2: ~0.70 (strong coherence — all stars elongated same direction)
p: ~0 (exp(-70))

Path A: R^2 = 0.70 > 0.5 AND p < 0.01 --> FIRES
Path B: median_ecc = 0.45 < 0.6 --> does not fire

Result: possibly_trailed = true, trail_r_squared = 0.70
        Full metrics still computed — caller decides whether to reject.
```

### Trailed image, undersampled (wind shake)

```
Detected: 80 stars
Median eccentricity: ~0.72
R^2: ~0.65
p: ~0

Path A: R^2 = 0.65 > 0.5 --> FIRES
Path B: median_ecc = 0.72 > 0.6 AND p < 0.05 --> ALSO FIRES

Result: possibly_trailed = true, trail_r_squared = 0.65
        Full metrics still computed.
```

---

## Why R^2 = 0.5 as the Default Threshold?

The threshold needs to be:
- **Above** non-trail coherence sources (grid effects, field-angle aberrations)
- **Below** real trail coherence (> 0.7 for tracking/wind trails)

Non-trail coherence sources measured on real data:
- Undersampled round stars (FWHM ~ 2 px): R^2 ~ 0.04-0.24 (grid-induced)
- Oversampled stars with field-angle effects (FWHM ~ 5 px): R^2 ~ 0.32-0.40
  (coma, field curvature create systematic elongation patterns that are NOT trailing)

The value 0.5 sits above both regimes. It corresponds to roughly 70% of stars
having coherent theta, which only occurs with real tracking errors or wind shake.

For reference, the expected R^2 under uniformity is 1/n (e.g., 0.005 for n=200).

---

## Phase 2: Per-Star Eccentricity Filter

After PSF measurement (downstream of the Rayleigh test), individual stars with
eccentricity above `max_eccentricity` (default 0.5) are removed from the final
result. This catches:

- Individual trailed stars in an otherwise good image
- Cosmic rays that survived the 3:1 aspect ratio filter
- Optical artifacts at field edges (coma, astigmatism)

This is configured via `ImageAnalyzer::with_max_eccentricity()`. Set to 1.0 to
disable.

---

## Position Angle (Theta) Source

The theta values used in the Rayleigh test come from stamp-based intensity-weighted
second-order moments computed during detection (see [detection.md](detection.md)).

For PSF measurement (metrics), theta is also computed from moments rather than the
Gaussian fit — the moments-based theta is defined for all star shapes (including
nearly round ones) and carries directional signal even for barely-elongated trail
knots, while the Gaussian fitter only fits theta for obviously elliptical stars.

---

## Constants

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Minimum stars | 5 | Below this, circular statistics unreliable |
| Path A R^2 threshold | 0.5 (default, configurable) | Above non-trail coherence (~0.40), below real trails (>0.7) |
| Path A p threshold | 0.01 | Strict — high confidence required |
| Path B ecc threshold | 0.6 | Above undersampled baseline (~0.3-0.4) |
| Path B p threshold | 0.05 | Standard significance level |
| Phase 2 max ecc | 0.5 (default) | Per-star filter, configurable |
