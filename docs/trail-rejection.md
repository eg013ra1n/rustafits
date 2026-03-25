# Trail Detection (Two-Stage Rayleigh + Eccentricity)

Image-level detection of satellite trails, tracking errors, wind shake, and vibration
using circular statistics on star position angles and PSF-fit eccentricity. Reports an
advisory flag and raw R² statistic — the caller decides whether to reject.

> **Changed in v0.7.2:** Trail detection is now **two-stage**. Stage 1 uses the Rayleigh
> test on detection-stage moments (angle coherence). Stage 2 refines the result after
> PSF measurement using the accurate fitted eccentricity — catching non-coherent
> guiding issues (wind shake, vibration) that the angle-based test cannot detect.
> The minimum star count was raised from 5 to 20 for statistical reliability.

## The Problem

Trailed images (from tracking drift, wind, or cable snag) produce elongated star
profiles. **Coherent** trailing (RA drift) stretches all stars in the same direction.
**Non-coherent** issues (wind shake, vibration, periodic error) elongate stars in
random or semi-random directions. A single detection method cannot catch both.

## Stage 1: Rayleigh Test on Detection-Stage Moments

The [Rayleigh test](https://en.wikipedia.org/wiki/Rayleigh_test) detects non-uniformity
in a circular distribution. We use it on the doubled position angles (2θ)
because theta is an axial quantity with period π (an angle of 0 and π represent
the same orientation).

### Statistic

Given n detected stars with position angles θ₁ … θₙ:

```
sum_cos = Σ cos(2θᵢ)
sum_sin = Σ sin(2θᵢ)

R̄² = (sum_cos² + sum_sin²) / n²

p = exp(-n · R̄²)
```

R̄² is the squared mean resultant length — it measures how concentrated the angles
are. R̄² = 0 means perfectly uniform (random angles), R̄² = 1 means all angles
identical (strong trail). Note: R̄² = 0.5 corresponds to R̄ ≈ 0.71 (71% directional
coherence).

The p-value gives the probability of observing R̄² this large or larger under the
null hypothesis of uniform random angles. The asymptotic formula `exp(-n·R̄²)` is
reliable for n ≥ 20.

### Minimum Star Count

The test requires at least **20** detected stars. Below this, the asymptotic p-value
is unreliable and the test has insufficient statistical power. When skipped:
`trail_r_squared = 0.0`, `possibly_trailed = false`.

### Median Eccentricity

The detection-stage median eccentricity is computed from all detected stars using a
proper median (average of the two middle values for even-length arrays).

---

## Dual-Path Rayleigh Detection

A single threshold does not work for all images because undersampled stars
(FWHM < 3 px) produce grid-induced theta coherence — the pixel grid geometry
biases the moment tensor, making round stars appear to have a preferred angle.

This creates a regime-dependent false positive problem:

| Regime | Stars look like | Theta behavior | Risk |
|--------|----------------|----------------|------|
| Undersampled, round | Small, boxy blobs | Grid-coherent (R̄² ~ 0.15) | False positive |
| Undersampled, trailed | Elongated blobs | Trail-coherent (R̄² > 0.5) | Correctly detected |
| Oversampled, round | Smooth, circular | Random (R̄² ~ 1/n) | No risk |
| Oversampled, trailed | Smooth, elongated knots | Trail-coherent (R̄² > 0.5) | Missed by ecc gate |

The solution is two detection paths:

### Path A — Strong Directional Coherence (R̄² > threshold)

```
Flag if:  R̄² > threshold  AND  p < 0.01
```

The threshold defaults to 0.5 and is configurable via `with_trail_threshold()`.

Real satellite/tracking trails produce R̄² > 0.7 (strong coherence). Grid-induced
coherence on non-trailed stars varies by regime:
- Undersampled round stars (FWHM ~ 2 px): R̄² ~ 0.15
- Oversampled stars with field-angle effects (FWHM ~ 5 px): R̄² up to 0.40

The default 0.5 threshold sits above both regimes while still catching real trails.

### Path B — Eccentricity-Gated Rayleigh (median ecc > 0.6)

```
Flag if:  R̄² > 0.05  AND  median_eccentricity > 0.6  AND  p < 0.05
```

For undersampled stars, stamp-based moments are noisy and theta has grid bias.
This path only runs when blobs are **genuinely elongated** (median eccentricity > 0.6)
AND there is at least modest directional coherence (R̄² > 0.05). The R̄² floor
prevents false triggers at very high star counts where even tiny spurious coherence
yields p < 0.05.

### Combined Rayleigh Decision (Stage 1)

```
rayleigh_trailed = (R̄² > threshold AND p < 0.01)     -- Path A
                OR (R̄² > 0.05 AND median_ecc > 0.6 AND p < 0.05)  -- Path B
```

Both paths compute the same Rayleigh statistic from the same star list; they
differ only in their activation criteria and p-value threshold.

---

## Stage 2: Post-Measurement Eccentricity Refinement

After PSF fitting, the pipeline computes the median of PSF-fit eccentricities
(which are much more accurate than detection-stage moments). If the fitted median
eccentricity exceeds **0.55**, the frame is flagged as trailed regardless of the
Rayleigh test result.

```
fit_eccs = [measured_star.eccentricity for all PSF-fit stars]
fit_median_ecc = median(fit_eccs)

possibly_trailed = rayleigh_trailed OR (fit_median_ecc > 0.55)
```

This catches **non-coherent guiding issues** (wind shake, vibration, periodic error)
where star elongation angles are random but all stars are clearly elongated. The
Rayleigh test fundamentally cannot detect these because there is no preferred direction.

### Why 0.55?

Good frames typically have PSF-fit median eccentricity of 0.30–0.50. Frames with
guiding issues show 0.55+. The threshold sits above the normal range while catching
genuinely problematic frames. This was calibrated against both dense M42 frames
(77 files, FWHM R²=0.995 vs PixInsight) and sparse EDPH data (117 files).

---

## Full Decision Logic

```
Stage 1 (detection-stage, before PSF measurement):
    if n_detected >= 20:
        compute R̄², p, detection-stage median_ecc
        rayleigh_trailed = (R̄² > 0.5 && p < 0.01)
                        || (R̄² > 0.05 && median_ecc > 0.6 && p < 0.05)
    else:
        rayleigh_trailed = false

```

**Note (v0.8.2):** The previous Stage 2 PSF-ecc override (`fit_median_ecc > 0.55`)
has been removed. Trail detection now uses the Rayleigh angle coherence test only.
High eccentricity without directional coherence indicates optical aberration (coma,
tilt) or wind shake — not tracking drift — and should not trigger trail mode.

When `possibly_trailed` is true, the statistics computation bypasses the
eccentricity ≤ 0.8 filter for FWHM and HFR medians. This ensures that
trailed frames report accurate (high) values rather than being suppressed.

---

## Effect on Statistics

When `possibly_trailed = true`:

| Statistic | Normal frames | Trailed frames |
|-----------|--------------|----------------|
| FWHM | ecc ≤ 0.8 filter + residual weighting | No ecc filter, residual weighting only |
| Eccentricity | Residual weighting | Residual weighting |
| HFR | ecc ≤ 0.8 filter + residual weighting | No ecc filter, residual weighting only |
| SNR, beta | Always unfiltered | Always unfiltered |

This prevents the ecc ≤ 0.8 cutoff from silently suppressing the FWHM/HFR
signal on genuinely trailed frames.

---

## API Surface

The Rayleigh test result is exposed via two fields on `AnalysisResult`:

| Field | Type | Description |
|-------|------|-------------|
| `trail_r_squared` | `f32` | Raw R̄² statistic. 0.0 = uniform, 1.0 = perfectly aligned. |
| `possibly_trailed` | `bool` | True if any detection stage fired (Rayleigh OR fit ecc > 0.55). |

The R̄² threshold for Path A is configurable:

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
Detection-stage median ecc: ~0.48
R̄²: ~0.15 (grid-induced)
p: ~0 (statistically significant but weak)

Path A: R̄² = 0.15 < 0.5 --> does not fire
Path B: R̄² = 0.15 > 0.05 but median_ecc = 0.48 < 0.6 --> does not fire
Stage 2: PSF-fit median ecc = 0.36 < 0.55 --> does not fire

Result: possibly_trailed = false, trail_r_squared = 0.15
        Statistics use ecc ≤ 0.8 filter as normal.
```

### Trailed image, coherent (tracking drift)

```
Detected: 100 stars
Detection-stage median ecc: ~0.45
R̄²: ~0.70 (strong coherence — all stars elongated same direction)
p: ~0 (exp(-70))

Path A: R̄² = 0.70 > 0.5 AND p < 0.01 --> FIRES
Stage 2: not needed (already flagged)

Result: possibly_trailed = true, trail_r_squared = 0.70
        Statistics bypass ecc filter — reports true eccentricity.
```

### Trailed image, non-coherent (wind shake)

```
Detected: 500 stars
Detection-stage median ecc: ~0.55
R̄²: ~0.02 (random angles — no coherence)
p: ~0 (due to large n, but R̄² is tiny)

Path A: R̄² = 0.02 < 0.5 --> does not fire
Path B: R̄² = 0.02 < 0.05 --> does not fire
Stage 2: PSF-fit median ecc = 0.58 > 0.55 --> FIRES

Result: possibly_trailed = true, trail_r_squared = 0.02
        Statistics bypass ecc filter — reports true eccentricity.
```

### Borderline image (mild guiding wobble)

```
Detected: 300 stars
Detection-stage median ecc: ~0.52
R̄²: ~0.08
p: ~0

Path A: R̄² = 0.08 < 0.5 --> does not fire
Path B: R̄² = 0.08 > 0.05 but median_ecc = 0.52 < 0.6 --> does not fire
Stage 2: PSF-fit median ecc = 0.47 < 0.55 --> does not fire

Result: possibly_trailed = false, trail_r_squared = 0.08
        Slightly elevated R̄² hints at some wobble, but not enough to flag.
```

---

## Why R̄² = 0.5 as the Default Rayleigh Threshold?

The threshold needs to be:
- **Above** non-trail coherence sources (grid effects, field-angle aberrations)
- **Below** real trail coherence (> 0.7 for tracking/wind trails)

Non-trail coherence sources measured on real data:
- Undersampled round stars (FWHM ~ 2 px): R̄² ~ 0.04-0.24 (grid-induced)
- Oversampled stars with field-angle effects (FWHM ~ 5 px): R̄² ~ 0.32-0.40
  (coma, field curvature create systematic elongation patterns that are NOT trailing)

The value 0.5 sits above both regimes. It corresponds to R̄ ≈ 0.71 — roughly 71%
of stars having coherent theta, which only occurs with real tracking errors.

For reference, the expected R̄² under uniformity is 1/n (e.g., 0.005 for n=200).

---

## Position Angle (Theta) Source

The theta values used in the Rayleigh test come from stamp-based intensity-weighted
second-order moments computed during detection (see [detection.md](detection.md)).
The range is (−π/2, π/2] — an axial orientation, not a full-circle direction.
The 2θ doubling maps this to (−π, π] for the circular Rayleigh test.

For PSF measurement (metrics), theta is also computed from moments rather than the
Gaussian fit — the moments-based theta is defined for all star shapes (including
nearly round ones) and carries directional signal even for barely-elongated trail
knots, while the Gaussian fitter only fits theta for obviously elliptical stars.

---

## Constants

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Minimum stars | 20 | Below this, asymptotic p-value unreliable |
| Path A R̄² threshold | 0.5 (default, configurable) | Above non-trail coherence (~0.40), below real trails (>0.7) |
| Path A p threshold | 0.01 | Strict — high confidence required |
| Path B R̄² floor | 0.05 | Prevents false triggers at high n with spurious coherence |
| Path B ecc threshold | 0.6 | Above undersampled baseline (~0.3-0.4) |
| Path B p threshold | 0.05 | Standard significance level |
| Stage 2 fit ecc threshold | 0.55 | Above normal PSF-fit range (0.30-0.50), catches wind shake |
