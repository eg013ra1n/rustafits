# Star Detection

Two-stage pipeline: DAOFIND-style matched-filter convolution followed by connected
component labeling (CCL) with union-find.

## Stage 1: Matched Filter & Peak Detection

```
Input: luminance image, background, noise
       |
       v
+-------------------------------+
| Build Zero-Sum Gaussian Kernel|
|                               |  FWHM_est = 3.0 px
|   sigma_k = FWHM / 2.3548    |  sigma_k ~ 1.274
|   radius = ceil(2 * sigma_k)  |  kernel size = 2*r + 1
|                               |
|   K(x,y) = exp(-(x^2+y^2) / (2*sigma_k^2))
|   K(x,y) -= mean(K)          |  zero-sum: removes DC bias
|   energy = sqrt(sum(K^2))     |  normalization factor
+-------------------------------+
       |
       v
+-------------------------------+
| Compute Detection Threshold   |
|                               |  threshold = detection_sigma * noise * energy
|   default detection_sigma=5.0 |
+-------------------------------+
       |
       v
+-------------------------------+
| Convolve Image with Kernel    |  conv(x,y) = sum( lum(x+dx, y+dy) * K(dx,dy) )
|                               |  Skip border pixels (width = radius)
+-------------------------------+
       |
       v
+-------------------------------+
| Local Maximum Detection       |  For each pixel:
|                               |    conv(x,y) > threshold  AND
|                               |    conv(x,y) > all 8 neighbors
|                               |  Skip border pixels (radius + 1)
+-------------------------------+
       |
       v
+-------------------------------+
| Non-Maximum Suppression       |  Sort peaks by conv value descending
|                               |  For each peak, suppress others within kernel_radius:
|                               |    reject if sqrt(dx^2 + dy^2) < kernel_radius
+-------------------------------+
       |
       v
  List of peak positions (pixel coordinates)
```

### Why a Zero-Sum Kernel?

A standard Gaussian kernel responds to both stars and flat backgrounds. Subtracting
the kernel mean makes it a bandpass filter that responds only to point-source-like
features at the expected FWHM scale:

```
Standard:  [0.05  0.12  0.20  0.12  0.05]   responds to constant signal
Zero-sum:  [-0.06  0.01  0.09  0.01  -0.06]  rejects constant, detects peaks
```

---

## Stage 2: Connected Component Labeling

```
Input: luminance image, peak positions, background, noise
       |
       v
+-------------------------------+
| Binary Threshold              |  threshold_low = background + 1.5 * noise
|                               |  Foreground: pixel > threshold_low
|                               |  (lower than detection threshold to capture star wings)
+-------------------------------+
       |
       v
+-------------------------------+
| Two-Pass CCL (Union-Find)     |
|                               |  Pass 1: scan L->R, T->B
|   +---+---+---+               |    - If left or top neighbor labeled: inherit or union
|   | . | T | . |               |    - Else: new label
|   +---+---+---+               |
|   | L | * |   |               |  Pass 2: resolve all labels to roots
|   +---+---+---+               |
+-------------------------------+
       |
       v
+-------------------------------+
| Link Peaks to Components      |  Each peak matched to its component via 3x3 neighborhood
|                               |  Only components containing a detected peak are kept
+-------------------------------+
       |
       v
+-------------------------------+
| Per-Component Metrics         |
|                               |  area = pixel count
|   Intensity-squared centroid: |
|     w_i = max(0, pixel - bg)^2
|     cx = sum(w_i * x_i) / sum(w_i)
|     cy = sum(w_i * y_i) / sum(w_i)
|                               |
|   peak = max(pixel - bg)      |
|   flux = sum(pixel - bg)      |
|   bbox = bounding rectangle   |
+-------------------------------+
       |
       v
+-------------------------------+
| Filtering                     |
|                               |  Reject if:
|   [x] area < min_star_area    |    (hot pixels, default 5 px)
|   [x] area > max_star_area    |    (galaxies/nebulae, default 2000 px)
|   [x] touches image border    |    (truncated profiles)
|   [x] peak > saturation_limit |    (default 0.95 * 65535 = 62163 ADU)
|   [x] aspect_ratio > 3.0      |    (cosmic rays, satellite trails)
|                               |
|   aspect_ratio = max(bbox_w, bbox_h) / min(bbox_w, bbox_h)
+-------------------------------+
       |
       v
+-------------------------------+
| Sort by Flux Descending       |  Keep top max_stars (default 200)
+-------------------------------+
       |
       v
  List of DetectedStar { x, y, peak, flux, area, bbox }
```

### Why I^2-Weighted Centroid?

Standard intensity-weighted centroids (`w = I`) are biased by faint wings and background
noise. Squaring the weights (`w = I^2`) concentrates the centroid estimate on the bright
core, giving more accurate subpixel positions for undersampled stars.

### Constants

| Parameter         | Value        | Rationale                                 |
|-------------------|--------------|-------------------------------------------|
| Initial FWHM      | 3.0 px       | Typical well-sampled star size            |
| Detection sigma   | 5.0 (default)| Configurable; 5-sigma = very few false positives |
| Low CCL threshold | 1.5 * noise  | Captures star wings for area measurement  |
| Min star area     | 5 px         | Rejects hot pixels and cosmic ray hits    |
| Max star area     | 2000 px      | Rejects extended objects                  |
| Max aspect ratio  | 3.0          | Rejects elongated artifacts               |
| Saturation        | 0.95 * 65535 | Saturated stars have unreliable profiles  |
| Max stars         | 200          | Performance cap (configurable)            |
