# SNR Computations

Two levels of signal-to-noise measurement: per-star aperture photometry
and derived quality metrics (SNR weight, PSF signal).

## Per-Star SNR (Aperture Photometry)

Classical CCD noise equation applied per star with a local sky annulus.

```
                        r_outer
                    .  .  .  .  .
                 .  .  .  .  .  .  .
              .  .  . r_inner .  .  .
              .  .  . . . . . .  .  .
              .  . +-----------+ .  .
              .  . |  r_ap     | .  .       Aperture: star flux
              .  . |    [*]    | .  .       Annulus:  local sky
              .  . |           | .  .
              .  . +-----------+ .  .
              .  .  . . . . . .  .  .
              .  .  .  .  .  .  .  .
                 .  .  .  .  .  .  .
                    .  .  .  .  .

  r_ap    = max(1.5 * median_FWHM, 3.0)    star aperture
  r_inner = 3.0 * median_FWHM              sky annulus inner edge
  r_outer = 5.0 * median_FWHM              sky annulus outer edge
```

### Algorithm

```
Input: luminance image, star centroid (cx, cy), aperture radii
       |
       v
+--------------------------------+
| Collect Pixels                 |
|   Aperture: dist <= r_ap       |  -> n_ap pixels
|   Annulus:  r_inner <= dist <= r_outer  -> annulus values
+--------------------------------+
       |
       v
+--------------------------------+
| Sigma-Clip Annulus             |  2 rounds, 3-sigma
|   median -> MAD -> sigma       |  (same MAD-based method as background)
|   clip outliers (stars in annulus)
+--------------------------------+
       |
       v
+--------------------------------+
| Local Sky Estimation           |
|   bg_local  = median(clipped annulus)
|   sigma_local = 1.4826 * MAD(clipped annulus)
+--------------------------------+
       |
       v
+--------------------------------+
| Star Flux                      |
|   F_star = sum( max(0, pixel - bg_local) )   over aperture pixels
+--------------------------------+
       |
       v
+--------------------------------+
| CCD Noise Equation             |
|                                |
|   noise^2 = F_star + n_ap * sigma_local^2
|             ~~~~~~   ~~~~~~~~~~~~~~~~~~~~
|             photon   sky background noise
|             noise    (n pixels * variance per pixel)
|                                |
|   SNR = F_star / sqrt(noise^2) |
+--------------------------------+
       |
       v
  per-star SNR (linear, dimensionless)
```

### CCD Noise Equation Breakdown

```
SNR = F_star / sqrt(F_star + n_ap * sigma_bg^2)

Term 1: F_star           Photon (shot) noise from the star itself.
                          Variance = F for Poisson statistics.

Term 2: n_ap * sigma^2   Sky background noise contribution.
                          Each of the n_ap aperture pixels contributes
                          sigma^2 of background variance.

Not included (assumed negligible for modern CMOS):
  - Read noise (R^2 * n_ap)
  - Dark current (D * t * n_ap)
```

### Validation

Returns SNR = 0 if:
- Star is too close to image edge (aperture truncated)
- Annulus has fewer than 5 pixels after clipping
- Star flux is zero or negative

### Constants

| Parameter       | Value          | Rationale                              |
|-----------------|----------------|----------------------------------------|
| Aperture radius | 1.5 * FWHM    | Captures ~90% of Gaussian flux         |
| Annulus inner   | 3.0 * FWHM    | Far enough to avoid star wings         |
| Annulus outer   | 5.0 * FWHM    | Wide enough for statistical sampling   |
| Min aperture    | 3.0 px         | Floor for very tight PSFs              |
| Clip rounds     | 2              | Remove stars contaminating the annulus |
| Clip threshold  | 3.0 sigma      | Standard astronomical clipping         |
| Min annulus px  | 5              | Need minimum samples for statistics    |

---

## SNR Weight

A star-based frame-ranking metric for subframe evaluation and stacking weight
computation. Immune to background gradients (signal measured from stars only).

```
SNR_weight = median(star_flux)^2 / (noise^2 * background)
```

### Algorithm

```
Input: measured stars (with flux), background, noise
       |
       v
+--------------------------------+
| Collect Star Fluxes            |  background-subtracted integrated flux per star
| Compute Median                 |  robust to outlier stars (saturated, blended)
+--------------------------------+
       |
       v
+--------------------------------+
| Compute Weight                 |
|   SNR_weight = median_flux^2   |
|              / (noise^2 * bg)  |
+--------------------------------+
```

Higher values indicate better frame quality. Returns 0 when no stars are
detected (correct behavior for unusable frames like clouded-out exposures).

The numerator uses stellar flux (not pixel deviation), so background gradients
from light pollution or vignetting cannot inflate the metric. The denominator
penalizes both noisy frames (noise^2) and elevated backgrounds from moonlight
or poor transparency (background).

---

## PSF Signal

Measures typical star brightness relative to the noise floor:

```
PSF_signal = median(star_peaks) / noise
```

### Algorithm

```
Input: measured stars (with peak values), noise
       |
       v
+--------------------------------+
| Collect Peaks                  |  background-subtracted peak per star
| Compute Median                 |
| Divide by Noise                |
+--------------------------------+
```

**Interpretation:** How many noise-sigmas above background is a typical star peak.
Higher = brighter stars relative to noise = better data for photometry and astrometry.

Some subframe evaluation tools normalize this with a user-defined divisor to produce
a [0, 1] weight. That divisor is a user preference, not computed from data.
