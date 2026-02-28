# Trail Rejection (Rayleigh Test)

Image-level detection of satellite trails, tracking errors, and wind shake using
circular statistics on star position angles. Rejects entire images where stars
show coherent directional elongation.

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
for meaningful circular statistics and the test is skipped.

---

## Dual-Path Rejection

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

The solution is two rejection paths:

### Path A — Strong Directional Coherence (R^2 > 0.3)

```
Reject if:  R^2 > 0.3  AND  p < 0.01
```

Real satellite/tracking trails produce R^2 > 0.5 (strong coherence). Grid-induced
coherence on undersampled round stars gives R^2 ~ 0.15, well below the 0.3 threshold.

This path catches **oversampled trails** where individual CCL knots have low
eccentricity (~0.3-0.5) but consistent theta. Without this path, the eccentricity
gate (Path B) would never fire and such trails would slip through.

### Path B — Eccentricity-Gated Rayleigh (median ecc > 0.6)

```
Reject if:  median_eccentricity > 0.6  AND  p < 0.05
```

For undersampled stars, stamp-based moments are noisy and theta has grid bias.
This path only runs the Rayleigh test when blobs are **genuinely elongated**
(median eccentricity > 0.6), which exceeds the undersampled baseline (~0.3-0.4).

This is the original gate — it catches **undersampled trails** where eccentricity
is clearly above the noise floor.

### Combined Decision

```
reject = (R^2 > 0.3 AND p < 0.01)           -- Path A
      OR (median_ecc > 0.6 AND p < 0.05)    -- Path B
```

Both paths compute the same Rayleigh statistic from the same star list; they
differ only in their activation criteria and p-value threshold.

---

## How It Works: Scenarios

### Good image, undersampled stars (FWHM ~ 2 px)

```
Detected: 200 stars
Median eccentricity: ~0.39
R^2: ~0.15 (grid-induced)
p: ~0 (statistically significant but weak)

Path A: R^2 = 0.15 < 0.3 --> does not fire
Path B: median_ecc = 0.39 < 0.6 --> does not fire

Result: PASSES (correct — this is a good image)
```

### Good image, oversampled stars (FWHM ~ 5 px)

```
Detected: 150 stars
Median eccentricity: ~0.15
R^2: ~0.01 (random theta)
p: ~0.22

Path A: R^2 = 0.01 < 0.3 --> does not fire
Path B: median_ecc = 0.15 < 0.6 --> does not fire

Result: PASSES (correct)
```

### Trailed image, oversampled (tracking drift)

```
Detected: 100 stars
Median eccentricity: ~0.45
R^2: ~0.70 (strong coherence — all stars elongated same direction)
p: ~0 (exp(-70))

Path A: R^2 = 0.70 > 0.3 AND p < 0.01 --> FIRES
Path B: median_ecc = 0.45 < 0.6 --> does not fire

Result: REJECTED by Path A (correct — oversampled trail caught)
```

### Trailed image, undersampled (wind shake)

```
Detected: 80 stars
Median eccentricity: ~0.72
R^2: ~0.55
p: ~0

Path A: R^2 = 0.55 > 0.3 --> FIRES
Path B: median_ecc = 0.72 > 0.6 AND p < 0.05 --> ALSO FIRES

Result: REJECTED by both paths (correct — obvious trail)
```

---

## Why R^2 = 0.3 as the Path A Threshold?

The threshold needs to be:
- **Above** grid-induced coherence (~0.15 for typical undersampled images)
- **Below** real trail coherence (> 0.5 for even moderate trails)

The value 0.3 sits in the gap between these regimes. It corresponds to roughly 55%
of stars having coherent theta, which is unlikely to occur from noise or grid effects
alone but is easily achieved by a real tracking error.

For reference, the expected R^2 under uniformity is 1/n (e.g., 0.005 for n=200),
and grid-induced bias adds ~0.10-0.15 on top of that.

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
| Path A R^2 threshold | 0.3 | Above grid noise (~0.15), below real trails (>0.5) |
| Path A p threshold | 0.01 | Strict — high confidence required |
| Path B ecc threshold | 0.6 | Above undersampled baseline (~0.3-0.4) |
| Path B p threshold | 0.05 | Standard significance level |
| Phase 2 max ecc | 0.5 (default) | Per-star filter, configurable |
