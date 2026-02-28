# Star Detection

Two-stage pipeline: DAOFIND-style matched-filter convolution followed by connected
component labeling (CCL) with union-find. Each detected star carries stamp-based
position angle (theta) and eccentricity for downstream trail rejection.

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

### Non-Maximum Suppression via Spatial Hashing

Peaks are sorted by convolution value (brightest first) and inserted into a spatial
grid with cell size equal to the kernel radius. For each candidate peak, only neighboring
grid cells need to be checked for suppression, giving O(n) instead of O(n^2) pairwise
comparisons.

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
| Stamp-Based Theta & Ecc      |  (see section below)
+-------------------------------+
       |
       v
+-------------------------------+
| Sort by Flux Descending       |  Keep top max_stars (default 200)
+-------------------------------+
       |
       v
  List of DetectedStar { x, y, peak, flux, area, theta, eccentricity }
```

---

## Stamp-Based Position Angle and Eccentricity

After computing the centroid, a stamp region around each star is used to compute
intensity-weighted second-order moments. These give the position angle (theta) and
eccentricity of the blob shape.

### Why Stamps Instead of CCL Pixels?

CCL blobs have too few pixels for reliable moments — especially for undersampled stars
(FWHM < 3 px), where the blob might be only 5-15 pixels. This causes strong
grid-induced theta coherence: the pixel geometry dominates the moment tensor.

Using a continuous stamp region (e.g. 11x11 or larger) over the background-subtracted
image includes many more pixels, and the intensity weighting naturally focuses on the
star's true shape rather than the CCL boundary geometry.

### Adaptive Stamp Radius

The stamp radius scales with blob area so that larger stars get appropriately sized
moment windows:

```
stamp_r = max(5, floor(2 * sqrt(area / pi)))
```

| CCL area | Effective FWHM | stamp_r | Stamp size |
|----------|----------------|---------|------------|
| 10 px    | ~2 px          | 5       | 11 x 11   |
| 50 px    | ~5 px          | 8       | 17 x 17   |
| 100 px   | ~7 px          | 11      | 23 x 23   |
| 200 px   | ~10 px         | 16      | 33 x 33   |

The minimum of 5 ensures undersampled stars still get the standard 11x11 window.
Larger stamps for oversampled stars capture the full PSF wings, giving more accurate
moments.

### Moment Computation

For each pixel (px, py) in the stamp around centroid (cx, cy):

```
v = max(0, data[py, px] - bg[py, px])    intensity weight

sf  += v               total flux
six += v * px          first moment X
siy += v * py          first moment Y
sixx += v * px^2       second moment XX
siyy += v * py^2       second moment YY
sixy += v * px * py    second moment XY
```

Central moments:

```
icx = six / sf        intensity-weighted centroid X
icy = siy / sf        intensity-weighted centroid Y

Mxx = sixx/sf - icx^2     central moment XX (variance along X)
Myy = siyy/sf - icy^2     central moment YY (variance along Y)
Mxy = sixy/sf - icx*icy   central moment XY (covariance)
```

### Position Angle (Theta)

The angle of the major axis of the intensity distribution:

```
theta = 0.5 * atan2(2 * Mxy, Mxx - Myy)
```

Range: [-pi/2, pi/2] radians. Measures the direction of elongation counter-clockwise
from the +X axis.

For round stars theta is dominated by noise — this is expected and handled by the
Rayleigh test's statistical framework (random angles cancel out).

### Eccentricity

From the eigenvalues of the 2x2 moment matrix:

```
trace = Mxx + Myy
det   = Mxx * Myy - Mxy^2
disc  = trace^2 - 4 * det

lambda_1 = (trace + sqrt(disc)) / 2     major axis variance
lambda_2 = (trace - sqrt(disc)) / 2     minor axis variance

eccentricity = sqrt(max(0, 1 - lambda_2 / lambda_1))
```

Values: 0 = perfectly circular, approaching 1 = highly elongated.

### Why I^2-Weighted Centroid?

Standard intensity-weighted centroids (`w = I`) are biased by faint wings and background
noise. Squaring the weights (`w = I^2`) concentrates the centroid estimate on the bright
core, giving more accurate subpixel positions for undersampled stars.

---

## Output: DetectedStar

Each detection carries everything needed for downstream trail rejection and PSF
measurement:

| Field         | Type    | Description                                    |
|---------------|---------|------------------------------------------------|
| `x`           | f32     | Intensity-weighted centroid X (subpixel)       |
| `y`           | f32     | Intensity-weighted centroid Y (subpixel)       |
| `peak`        | f32     | Background-subtracted peak value               |
| `flux`        | f32     | Total background-subtracted flux               |
| `area`        | usize   | Number of connected pixels above threshold     |
| `theta`       | f32     | Position angle from stamp moments (radians)    |
| `eccentricity`| f32     | Eccentricity from stamp moments (0=round)      |

---

## Constants

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
| Min stamp radius  | 5            | 11x11 stamp for smallest blobs           |
| Stamp scale       | 2 * sqrt(area/pi) | Scales with equivalent circular radius |
