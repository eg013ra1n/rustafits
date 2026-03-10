# Background Estimation

Background estimation uses a mesh-grid approach with auto-tuned cell size. Noise is
estimated via iterative MRS wavelet (B3-spline à trous, default 4 layers with significance
masking). After pass 1 calibration, a source-masked re-estimation is always performed for
cleaner statistics.

## Global Background Estimation

> **Note:** This algorithm is used internally within each mesh cell. It is no longer
> available as a standalone pipeline option.

**Algorithm:** Iterative sigma-clipped statistics with SExtractor-style mode estimation.

```
Input: luminance image (f32)
       |
       v
+---------------------+
| Subsample           |  stride = sqrt(n_pixels / 500k), skip 2px border
| ~500k pixels        |  only finite values
+---------------------+
       |
       v
+---------------------+
| Sigma-Clip          |  3 rounds, each:
| (3 rounds, 3-sigma) |    median = median(samples)
|                     |    MAD = median(|xi - median|)
|                     |    sigma = 1.4826 * MAD
|                     |    keep only: median - 3*sigma <= x <= median + 3*sigma
+---------------------+
       |
       v
+---------------------+
| Mode Estimate       |  if |mean - median| < 0.3 * sigma:
| (SExtractor)        |    mode = 2.5 * median - 1.5 * mean
|                     |  else (skewed):
|                     |    mode = median   (asymmetry fallback)
|                     |  noise = max(1.4826 * MAD, 0.001)
+---------------------+
       |
       v
  background (ADU), noise (ADU)
```

### Why SExtractor Mode?

For astronomical images the pixel distribution is asymmetric: a Gaussian background peak
with a long tail from stars and nebulosity. The mean is pulled high by bright pixels.
The SExtractor formula `mode = 2.5*median - 1.5*mean` corrects for this skew, estimating
the true background peak more accurately than either mean or median alone.

### Asymmetry Fallback

When `|mean - median| >= 0.3 * sigma`, the distribution is too skewed for the mode formula
to be reliable (e.g., a mesh cell dominated by bright nebulosity or a galaxy core). In this
case, the formula can extrapolate to unphysical values. Following SExtractor, we fall back
to the median, which is more robust under heavy contamination.

### MAD to Sigma Conversion

The Median Absolute Deviation (MAD) is a robust scale estimator with 50% breakdown point
(up to 50% of data can be corrupted before the estimate is affected):

```
MAD = median(|xi - median(x)|)
```

For a Gaussian distribution, `sigma = 1.4826 * MAD`. This factor (1/Phi^-1(3/4)) makes
MAD consistent with standard deviation while being insensitive to outliers (stars).

### Constants

| Parameter        | Value   | Rationale                            |
|------------------|---------|--------------------------------------|
| Target samples   | 500,000 | Balance speed vs. statistical power  |
| Border exclusion | 2 px    | Avoid edge artifacts                 |
| Clip rounds      | 3       | Converges for typical star densities |
| Clip threshold   | 3.0 sigma | Rejects >99.7% outliers per round  |
| MAD scale factor | 1.4826  | Gaussian consistency factor          |
| Asymmetry limit  | 0.3 sigma | SExtractor fallback threshold      |
| Minimum noise    | 0.001   | Prevent division by zero             |

---

## Mesh-Grid Background Estimation

The primary background estimation method. Handles large-scale gradients (vignetting,
light pollution, differential dust) by estimating the background locally per cell and
interpolating to full resolution.

Cell size is determined automatically: `auto_cell_size = max(16, max(w, h) / 32)`.

```
Input: luminance image, auto_cell_size
       |
       v
+---------------------------+
| Divide into nx * ny cells |  cell_size x cell_size (last row/col may be smaller)
+---------------------------+
       |
       v
+---------------------------+
| Per-Cell Sigma-Clip       |  Same 3-round, 3-sigma algorithm as global
|                           |  -> cell_bg, cell_sigma per cell
|                           |  If >30% of pixels clipped -> mark cell INVALID
+---------------------------+
       |
       v
+---------------------------+
| Fill Invalid Cells        |  Nearest-neighbor search (expanding rings)
|                           |  Copy bg & sigma from closest valid cell
+---------------------------+
       |
       v
+---------------------------+
| 3x3 Median Filter on Grid|  Suppress residual star contamination
|                           |  Each cell = median of 3x3 neighborhood
+---------------------------+
       |
       v
+---------------------------+
| Bicubic Spline Interp    |  Catmull-Rom spline: grid -> full resolution
|                           |  C1-continuous, uses 4x4 neighborhood per pixel
|                           |  Passes through grid values, no grid artifacts
|                           |  Parallelized per row with rayon
+---------------------------+
       |
       v
  background_map[y * width + x],
  noise_map[y * width + x],
  global_background = median(grid),
  global_noise = median(cell_sigmas)
```

### Invalid Cell Detection

A cell is marked invalid when sigma-clipping removes more than 30% of its pixels.
This indicates heavy star contamination — the remaining "background" pixels are
unreliable. The nearest valid neighbor provides a better estimate.

### Bicubic Spline Interpolation

Grid cell values represent the background at cell centers. Between centers, Catmull-Rom
bicubic spline interpolation produces a smooth background map with continuous first
derivatives (C1 continuity). This matches SExtractor/SEP's approach and eliminates the
grid-boundary artifacts that bilinear interpolation would produce.

For each pixel (x, y), the interpolation uses a 4x4 neighborhood of grid cells with
Catmull-Rom weights:

```
w(-1) = -0.5*t^3 +     t^2 - 0.5*t
w( 0) =  1.5*t^3 - 2.5*t^2 + 1.0
w(+1) = -1.5*t^3 + 2.0*t^2 + 0.5*t
w(+2) =  0.5*t^3 - 0.5*t^2
```

where t is the fractional position between cell centers. The 2D interpolation is
separable: compute 4 horizontal interpolations, then 1 vertical.

### Source-Masked Re-estimation

The pipeline always performs one source-masked re-estimation after pass 1 calibration.
After the initial background estimate and star detection, star pixels are masked and
the mesh-grid background is re-estimated for cleaner statistics:

1. Pass 1: estimate background normally (sources contaminate statistics).
2. Detect stars using the initial background.
3. Build a source mask: circle of r = 2.5 * FWHM around each detected star.
4. Re-estimate background excluding masked pixels (cleaner statistics).

This is particularly valuable for crowded fields where bright stars bias the per-cell
background estimate.

### Constants

| Parameter            | Value | Rationale                              |
|----------------------|-------|----------------------------------------|
| Min cell size        | 16 px | Smaller cells lack statistical samples |
| Invalid threshold    | 30%   | High clip fraction = star-contaminated |
| Grid filter          | 3x3 median | Smooth grid, reject cell-scale artifacts |
| Interpolation        | Catmull-Rom bicubic | C1-continuous, matches SExtractor/SEP |

---

## MRS Wavelet Noise Estimation

Always-on (default `noise_layers=4`). The MRS (Multiresolution Support) wavelet
method isolates the finest-scale fluctuations using the à trous (with holes) wavelet
transform with iterative significance masking. This avoids the problem where
sigma-clipped MAD overestimates noise by 30-40% on nebula-rich fields (M42, NGC 7000,
etc.) by conflating nebulosity with background fluctuations.

### Algorithm

```
Input: luminance image, noise_layers (default 4)
       |
       v
+-------------------------------+
| Layer 1: B3-Spline Smooth     |  Separable 5-tap kernel: [1/16, 1/4, 3/8, 1/4, 1/16]
|                               |  w1[i] = data[i] - smoothed[i]
|                               |  Finest-scale coefficients (noise + point sources)
+-------------------------------+
       |
       v
+-------------------------------+
| Initial sigma from w1         |  sigma1 = 1.4826 * MAD(w1)
| Mark significant pixels       |  |w1[i]| > 3*sigma1 → sig_mask[i] = true
+-------------------------------+
       |
       v
+-------------------------------+
| Layers 2..N: Dilated Smooth   |  Spacing = 2^(layer-1) between taps
|                               |  w_j = c_{j-1} - c_j (coarser wavelet coefficients)
|                               |  sigma_j from unmasked pixels only
|                               |  |w_j[i]| > 3*sigma_j → add to sig_mask
+-------------------------------+
       |
       v
+-------------------------------+
| Re-estimate from w1           |  Use only unmasked pixels (structure-free)
| Iterate until convergence     |  sigma-clip on w1, 4 rounds, <0.1% change
+-------------------------------+
       |
       v
  noise estimate (ADU), floor clamped to 0.001
```

### Why B3-Spline?

The B3-spline kernel `[1/16, 1/4, 3/8, 1/4, 1/16]` is the standard smoothing function
in the à trous wavelet transform. It is:

- **DC-preserving**: coefficients sum to exactly 1.0, so constant signals pass through unchanged
- **Separable**: the 2D smooth decomposes into horizontal + vertical 1D passes
- **Compact**: only 5 taps, so the smooth is fast even without SIMD
- **Smooth**: sufficient support to separate noise (1-2 px scale) from signal (stars, nebulae)

The first wavelet layer (`w1 = data - smoothed`) captures structure at the 1-2 pixel
scale — exactly where Poisson/read noise lives. Nebulosity, gradients, and extended
objects have much larger spatial scales and are almost entirely removed by the subtraction.

MRS wavelet noise is always-on (default 4 layers). The layer count is configurable via
`with_mrs_layers(n)`. Higher layers (2+) build a significance mask using dilated
convolution to identify and exclude structure at progressively larger scales, producing
a purer noise estimate from the layer-1 coefficients.

### Constants

| Parameter        | Value     | Rationale                            |
|------------------|-----------|--------------------------------------|
| B3-spline kernel | [1/16, 1/4, 3/8, 1/4, 1/16] | Standard à trous smoothing function |
| Kernel radius    | 2 px      | 5-tap separable convolution          |
| Boundary         | Reflected | Avoids edge discontinuities          |
| Target samples   | 500,000   | Balance speed vs. statistical power  |
| Border exclusion | 2 px      | Avoid edge artifacts                 |
| MAD scale factor | 1.4826    | Gaussian consistency factor          |
| Minimum noise    | 0.001     | Prevent division by zero             |

---

## Comparison with Professional Tools

| Feature | SExtractor/SEP | photutils | rustafits |
|---------|---------------|-----------|-----------|
| Mode estimator | 2.5m - 1.5mu | Pluggable (SExtractor, biweight, MMM) | 2.5m - 1.5mu |
| Asymmetry fallback | Yes (0.3sigma) | Depends on estimator | Yes (0.3sigma) |
| Background method | Mesh-grid | Mesh-grid | Mesh-grid (always, auto cell size) |
| Noise estimator | Clipped sigma | MAD or biweight scale | MRS wavelet (always-on, 4 layers) |
| Wavelet noise | — | — | B3-spline à trous (default 4 layers, significance masking) |
| Grid filtering | 3x3 median | Configurable median | 3x3 median |
| Interpolation | Bicubic spline | Bicubic spline (scipy zoom) | Catmull-Rom bicubic |
| Source masking | No (single pass) | Manual mask input | Always-on (1 re-estimation, auto-mask) |
| Parallelization | No | No (Python) | Yes (rayon per-row) |
