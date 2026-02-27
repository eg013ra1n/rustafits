# SNR Computations

Three levels of signal-to-noise measurement: per-star aperture photometry,
image-wide SNR in decibels, and derived quality metrics.

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

## Image-Wide SNR in Decibels

Whole-frame signal-to-noise, comparable to PixInsight's SNRViews script.

```
SNR_dB = 20 * log10(mean_signal / noise)
```

### Algorithm

```
Input: luminance image, noise (from background estimation)
       |
       v
+--------------------------------+
| Subsample Image                |  stride to ~500k pixels
|   mean = sum(pixels) / count   |
+--------------------------------+
       |
       v
+--------------------------------+
| SNR in Decibels                |
|   SNR_dB = 20 * log10(mean / noise)
+--------------------------------+
```

### Interpretation

| SNR dB | Linear Ratio | Image Quality        |
|--------|-------------|----------------------|
| 10 dB  | 3.2         | Very noisy           |
| 20 dB  | 10          | Noisy                |
| 30 dB  | 31.6        | Good                 |
| 40 dB  | 100         | Very clean           |
| 50 dB  | 316         | Exceptional          |

PixInsight SNRViews on `mono.fits`: **28.70 dB**
Our computation on `mono.fits`: **27.65 dB** (~1 dB difference, within methodology variance)

---

## SNR Weight (SubframeSelector-Compatible)

PixInsight SubframeSelector uses this for frame ranking:

```
SNR_weight = (MeanDeviation / noise)^2
```

### Algorithm

```
Input: luminance image, background, noise
       |
       v
+--------------------------------+
| Subsample Image                |  stride to ~500k pixels
|   MeanDev = mean(|pixel - background|)
+--------------------------------+
       |
       v
+--------------------------------+
| Compute Weight                 |
|   ratio = MeanDev / noise      |
|   SNR_weight = ratio^2         |
+--------------------------------+
```

For a pure Gaussian background (no signal), MeanDev ~ 0.8 * sigma, so
SNR_weight ~ 0.64. Higher values indicate more signal content.

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

In PixInsight SubframeSelector, the user sets a "PSF Signal Divisor" to normalize
this into a [0, 1] weight. That divisor is a user preference, not computed from data.
