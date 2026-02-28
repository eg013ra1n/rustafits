# PSF Metrics: FWHM, Eccentricity, HFR

Per-star shape measurements computed from extracted stamps around each detected star.

## Stamp Extraction

```
Input: detected star (centroid, area), luminance image, background
       |
       v
+------------------------------+
| Estimate Star Scale          |
|   sigma_est = sqrt(area/pi)  |  area -> equivalent circular sigma
|              * 0.5           |
|                              |
|   stamp_radius = ceil(4*sigma_est)
|   clamp: 8 <= radius <= 50  |
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

**Interpretation:** For a Gaussian PSF with FWHM=F, HFR ~ 0.67 * FWHM/2 ~ 0.34*F.
HFR is more robust than FWHM for non-Gaussian profiles (e.g. coma, defocused stars).

---

## FWHM Estimation

Two methods are available:

### Method 1: 2D Gaussian Fit (Default)

Uses full elliptical Gaussian fitting (see [fitting.md](fitting.md)):

```
Star Stamp
    |
    v
+-----------------------------+
| Collect PixelSamples        |  (x, y, intensity) for all stamp pixels > 0
| with background subtraction |
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
| Convert to FWHM             |
|   FWHM_x = 2.3548 * sigma_x|
|   FWHM_y = 2.3548 * sigma_y|
|   FWHM = sqrt(Fx * Fy)     |  geometric mean
+-----------------------------+
    |
    v
+-----------------------------+
| Eccentricity                |
|   min_s = min(sx, sy)       |
|   max_s = max(sx, sy)       |
|   e = sqrt(1 - min^2/max^2) |
+-----------------------------+
    |
    v
+-----------------------------+
| Theta (moments-based)       |  Always from I-weighted second moments,
|                             |  NOT from the Gaussian fitter's theta.
|                             |  See note below.
+-----------------------------+
```

**Accuracy:** Matches PixInsight FWHMEccentricity within ~2% on typical data.

**Note on theta:** The position angle is always computed from intensity-weighted
second-order moments over the stamp, even when using Gaussian fitting for FWHM.
The Gaussian fitter only fits theta for stars with ellipticity > 0.1 and may
return 0 for round sources. The moments-based theta is defined for all shapes
and carries directional signal even for barely-elongated trail knots — this is
critical for the Rayleigh trail rejection test to work on oversampled trails.

### Method 2: Windowed Moments (Fast Path)

Non-parametric method using second-order intensity moments:

```
Star Stamp
    |
    v
+-----------------------------+
| Centroid Refinement (2 iter)|  Two-pass I-weighted centroid:
|                             |    cx = sum(val * x) / sum(val)
|                             |    cy = sum(val * y) / sum(val)
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

### Method Comparison (mono.fits)

| Metric           | Gaussian Fit | Moments | PixInsight |
|------------------|-------------|---------|------------|
| Median FWHM      | 2.18 px     | 4.91 px | 2.16 px    |
| Median Ecc       | 0.341       | 0.302   | 0.480      |

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
| Stamp radius min | 8 px      | Enough pixels for reliable moments      |
| Stamp radius max | 50 px     | Limit memory and computation            |
| Stamp scale      | 4 * sigma | Captures >99.9% of Gaussian energy      |
| FWHM factor      | 2.3548    | Gaussian sigma-to-FWHM conversion       |
| Centroid passes   | 2         | Convergence to subpixel accuracy        |
