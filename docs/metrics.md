# PSF Metrics: FWHM, Eccentricity, HFR

Per-star shape measurements computed from extracted stamps around each detected star.

## Stamp Extraction

```
Input: detected star (centroid, area), luminance image, background
       |
       v
+------------------------------+
| Estimate Star Scale          |
|   sigma_est = field_fwhm/2.3548  (if available)
|   fallback: sqrt(area/pi)*0.5    (from detection)
|                              |
|   stamp_radius = ceil(4*sigma_est)
|   clamp: 8 <= radius <= 50  |  Supports FWHM up to ~25px
+------------------------------+
       |
       v
+------------------------------+
| Extract Stamp                |  (2*radius+1) x (2*radius+1) pixel cutout
|                              |  centered on star centroid
|                              |  Background-subtract all pixels
|                              |  (uses bg_map if available, else global bg)
+------------------------------+
       |
       v
  stamp[], stamp_width, stamp_height, relative centroid (cx, cy)
```

## OSC (Bayer) Handling

For OSC (Bayer) images, the analysis pipeline uses green-channel interpolation to
produce a native-resolution mono image. All pixels in the interpolated image are
used for PSF fitting, moments computation, and HFR calculation. No green mask
filtering is applied.

## Adaptive Moments Screening

Before running the expensive Moffat/Gaussian LM fitter (~50µs per star), each
candidate is screened using cheap image moments (~1µs) to reject non-stellar objects.
The thresholds are **adaptive** — computed from the field's own statistics rather than
fixed constants.

### Two-Pass Architecture

```
Pass 0: Compute moments for ALL candidates (cheap, parallelized)
  For each star:
    - Extract stamp, compute 2-iteration intensity-weighted centroid
    - Compute second-order moments (Mxx, Myy, Mxy) → eigenvalues → FWHM, ecc
    - Compute sharpness = peak / mean_flux
  Collect: moment_fwhm[], moment_ecc[], moment_sharpness[]

  Derive adaptive thresholds from field statistics:
    median_ecc    = median(moment_ecc[])
    MAD_ecc       = 1.4826 × median(|ecc_i - median_ecc|), floored at 0.03
    median_sharp  = median(moment_sharpness[])
    MAD_sharp     = 1.4826 × median(|sharp_i - median_sharp|), floored at 0.1

  ScreeningThresholds:
    fwhm_lo   = 0.7 × field_fwhm
    fwhm_hi   = 2.5 × field_fwhm
    ecc_max   = max(0.85, median_ecc + 3 × MAD_ecc)
    sharp_lo  = max(0.15, median_sharp - 3 × MAD_sharp)
    sharp_hi  = min(8.0, median_sharp + 3 × MAD_sharp)

Pass 1: Measure with adaptive gates + Moffat/Gaussian fitting
  For each star:
    - Compute stamp moments (same as Pass 0)
    - Reject if moment_fwhm outside [fwhm_lo, fwhm_hi]
    - Reject if moment_ecc > ecc_max
    - Reject if sharpness outside [sharp_lo, sharp_hi]
    - Survivors → Moffat → Gaussian → Moments fallback chain
```

### Why Adaptive?

Fixed thresholds (e.g., ecc > 0.8, sharpness 0.3–5.0) fail on:

- **Trailed images** — median ecc ~0.77, so a fixed 0.8 threshold rejects most real stars.
  The adaptive ecc_max expands to ~0.77 + 3×0.03 = 0.86 (or higher if MAD is larger).
- **Large-FWHM frames** (>6px) — peak/mean_flux ratio is naturally lower for broad stars,
  falling below the fixed 0.3 sharpness floor. Adaptive sharp_lo follows the field median.
- **Uniform fields** — MAD can be very small (~0.01). The MAD floor prevents thresholds
  from collapsing: MAD_ecc floored at 0.03 (minimum ecc range ±0.09), MAD_sharp floored
  at 0.1 (minimum sharpness range ±0.3).

### Minimum Star Guarantee

After adaptive filtering, if fewer than 10 candidates survive PSF fitting (or fewer
than the total candidate count if < 10 candidates exist), the pipeline re-runs
measurement **without any screening**, then keeps the 10 stars closest to the field
median FWHM. This ensures PSF fitting always produces output — no zero-star results.

### Example: Normal Frame (cocoon.fits)

```
Candidates: 500
Moment statistics: median_ecc=0.35, MAD_ecc=0.12, median_sharp=1.8, MAD_sharp=0.4

Adaptive thresholds:
  ecc_max   = max(0.85, 0.35 + 0.36) = 0.85
  sharp_lo  = max(0.15, 1.8 - 1.2) = 0.60
  sharp_hi  = min(8.0, 1.8 + 1.2) = 3.00

Result: ~490/500 pass screening → measure with Moffat/Gaussian
```

### Example: Trailed Frame (m82.fits)

```
Candidates: 314
Moment statistics: median_ecc=0.78, MAD_ecc=0.06, median_sharp=0.9, MAD_sharp=0.3

Adaptive thresholds:
  ecc_max   = max(0.85, 0.78 + 0.18) = 0.96  ← expanded for trailing
  sharp_lo  = max(0.15, 0.9 - 0.9) = 0.15
  sharp_hi  = min(8.0, 0.9 + 0.9) = 1.80

Result: ~310/314 pass screening → measure with Moffat/Gaussian
```

---

## Half-Flux Radius (HFR)

A non-parametric measure of star size: the intensity-weighted mean distance from centroid.

```
        sum( w_i * d_i )
HFR = -------------------
           sum( w_i )

where:
  w_i = max(0, pixel_i - background)     intensity weight
  d_i = sqrt((x_i - cx)^2 + (y_i - cy)^2)   distance from centroid
```

**OSC handling:** For OSC images, all pixels in the green-interpolated image
contribute to the HFR sum (see [OSC Handling](#osc-bayer-handling)).

**Interpretation:** For a Gaussian PSF with FWHM=F, HFR ~ 0.67 * FWHM/2 ~ 0.34*F.
HFR is more robust than FWHM for non-Gaussian profiles (e.g. coma, defocused stars).

---

## FWHM Estimation

Two methods are available:

### Method 1: 2D Moffat Fit (Default — Two-Pass Calibration)

Uses full elliptical Moffat fitting for accurate PSF measurement. The Moffat profile
models power-law wings better than a Gaussian, giving more accurate FWHM especially
for well-sampled stars (FWHM > 3 px). See [fitting.md](fitting.md) for the full model.

The pipeline uses a two-pass calibration approach:

- **Pass 1** — free-beta Moffat (8 parameters) on top ~100 bright calibration stars
  (eccentricity < 0.5, area >= 5). Derives `field_beta` and `field_fwhm` via
  sigma-clipped median.
- **Pass 2** — fixed-beta Moffat (7 parameters, using `field_beta`) on all detected stars.
  This removes the beta/axis-ratio tradeoff, improving FWHM stability.

When a Moffat fit fails, the pipeline falls back to Gaussian, then to windowed moments.
Each star's `fit_method` field records which method produced the result (`FreeMoffat`,
`FixedMoffat`, `Gaussian`, or `Moments`). Each star also carries a `fit_residual` —
a normalized LM cost used as a quality weight for frame-level statistics
(see [fitting.md](fitting.md#fit-residual-quality-metric)).

### Method 2: 2D Gaussian Fit

Uses full elliptical Gaussian fitting (see [fitting.md](fitting.md)):

```
Star Stamp
    |
    v
+-----------------------------+
| Collect PixelSamples        |  (x, y, intensity) for stamp pixels > 0
| with background subtraction |  All pixels used (including OSC green-interpolated)
|                             |
+-----------------------------+
    |
    v
+-----------------------------+
| fit_gaussian_2d()           |  Levenberg-Marquardt, 2-stage
|   -> sigma_x, sigma_y      |  (see fitting.md for full algorithm)
+-----------------------------+
    |
    v
+-----------------------------+
| Canonicalize & Convert      |  Ensure fwhm_x >= fwhm_y, theta along major:
|                             |
|   if sigma_x >= sigma_y:   |
|     FWHM_x = 2.3548 * sx   |
|     FWHM_y = 2.3548 * sy   |
|     theta  = theta          |
|   else:                     |
|     FWHM_x = 2.3548 * sy   |  swap axes
|     FWHM_y = 2.3548 * sx   |
|     theta  = theta + pi/2  |  rotate to match
|                             |
|   FWHM = sqrt(Fx * Fy)     |  geometric mean (invariant)
+-----------------------------+
    |
    v
+-----------------------------+
| Eccentricity                |
|   e = sqrt(1 - Fy^2/Fx^2)  |  fwhm_x is guaranteed major
+-----------------------------+
```

**Accuracy:** Within ~1-2% of professional tools on typical data.

**Note on axis canonicalization:** The LM optimizer has a degeneracy — the pair
(sigma_x, sigma_y, theta) produces an identical PSF to (sigma_y, sigma_x,
theta + π/2). For nearly-round stars the theta gradient vanishes and the
optimizer freely drifts across the degeneracy ridge. The canonicalization step
ensures a consistent invariant (fwhm_x ≥ fwhm_y, theta along major axis) for
all downstream consumers (annotation ellipses, direction ticks). The same
canonicalization applies to Moffat fits. The moments path already produces
canonical output because eigenvalues are sorted by construction.

### Method 3: Windowed Moments (Fast Path)

Non-parametric method using second-order intensity moments:

```
Star Stamp
    |
    v
+-----------------------------+
| Centroid Refinement (2 iter)|  Two-pass I-weighted centroid:
|                             |    cx = sum(val * x) / sum(val)
|                             |    cy = sum(val * y) / sum(val)
|                             |  All pixels used (including OSC green-interpolated)
+-----------------------------+
    |
    v
+-----------------------------+
| Second Moments              |  Intensity-weighted covariance matrix:
|                             |
|   Mxx = sum(val*(x-cx)^2) / sum(val)    variance along x
|   Myy = sum(val*(y-cy)^2) / sum(val)    variance along y
|   Mxy = sum(val*(x-cx)*(y-cy)) / sum(val)  covariance
+-----------------------------+
    |
    v
+-----------------------------+
| Eigenvalue Decomposition    |  2x2 symmetric matrix -> analytical solution:
|                             |
|   trace = Mxx + Myy        |
|   det   = Mxx*Myy - Mxy^2  |
|   disc  = trace^2 - 4*det  |
|                             |
|   lambda1 = (trace + sqrt(disc)) / 2   major axis variance
|   lambda2 = (trace - sqrt(disc)) / 2   minor axis variance
+-----------------------------+
    |
    v
+-----------------------------+
| FWHM Conversion             |
|   FWHM_x = 2.3548 * sqrt(lambda1)
|   FWHM_y = 2.3548 * sqrt(lambda2)
|   FWHM   = sqrt(Fx * Fy)   |
+-----------------------------+
    |
    v
+-----------------------------+
| Eccentricity                |
|   e = sqrt(1 - lambda2/lambda1)
+-----------------------------+
```

**Trade-offs:**
- ~10x faster than Gaussian fit (no iterative solver)
- Systematically overestimates FWHM (noise in wings inflates moments)
- Less accurate eccentricity (no rotation model)

### Method Comparison

**mono.fits:**

| Metric           | Moffat (default) | Gaussian | Moments |
|------------------|-------------------|----------|---------|
| Median FWHM      | ~2.16 px          | ~2.14 px | ~4.91 px |
| Median Ecc       | ~0.34             | ~0.34    | ~0.30   |
| Reports beta     | Yes (~3.5)        | No       | No      |

**osc.fits (green-pixel-only fitting):**

| Metric           | Moffat (default) | Gaussian |
|------------------|-------------------|----------|
| Median FWHM      | ~2.66 px          | ~2.61 px |

---

## FWHM Conversion Factor

The constant 2.3548 converts Gaussian sigma to Full Width at Half Maximum:

```
For a Gaussian f(x) = exp(-x^2 / (2*sigma^2)):

  f(x) = 0.5  when  x = sigma * sqrt(2 * ln(2))

  FWHM = 2 * sigma * sqrt(2 * ln(2))
       = 2 * sigma * 1.17741
       = 2.3548 * sigma
```

---

## Eccentricity

Measures PSF elongation (tracking errors, optical aberrations):

```
                    +--------------------+
                    |                    |
  e = 0            |   e = 0.5          |   e = 0.9
  (circular)       |   (mildly oval)    |   (highly elongated)
                    |                    |
    * *            |     * *            |       * * *
  *     *          |   *     *          |   *           *
  *     *          |  *       *         |  *             *
  *     *          |   *     *          |   *           *
    * *            |     * *            |       * * *
                    |                    |
                    +--------------------+

Formula:
  e = sqrt(1 - (minor_axis / major_axis)^2)

  For Gaussian fit:  axes = sigma_x, sigma_y
  For moments:       axes = sqrt(lambda1), sqrt(lambda2)

  e = 0:  perfectly circular
  e -> 1: infinitely elongated
```

---

## Constants

| Parameter        | Value     | Rationale                               |
|------------------|-----------|-----------------------------------------|
| fit_residual     | f32       | Normalized LM residual (lower = better fit) |
| Stamp radius min | 8 px      | Enough pixels for reliable moments      |
| Stamp radius max | 50 px     | Limit memory and computation            |
| Stamp scale      | 4 * sigma | Captures >99.9% of Gaussian energy      |
| FWHM factor      | 2.3548    | Gaussian sigma-to-FWHM conversion       |
| Centroid passes   | 2         | Convergence to subpixel accuracy        |
