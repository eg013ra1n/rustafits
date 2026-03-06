# Background Estimation

Two methods are available: global (default) and mesh-grid (for images with gradients).

## Global Background Estimation

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

For images with large-scale gradients (vignetting, light pollution, differential dust).

```
Input: luminance image, cell_size (e.g. 64 px)
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

### Iterative Source-Masked Re-estimation

When `with_iterative_background(n)` is used with `with_background_mesh(cell_size)`,
the pipeline performs n rounds of background estimation:

1. First pass: estimate background normally (sources contaminate statistics).
2. Detect stars using the initial background.
3. Build a source mask: circle of r = 2.5 * FWHM around each detected star.
4. Re-estimate background excluding masked pixels (cleaner statistics).
5. Re-detect stars with improved background.

One iteration (n=2) is usually sufficient. This is particularly valuable for crowded
fields where bright stars bias the per-cell background estimate.

### Constants

| Parameter            | Value | Rationale                              |
|----------------------|-------|----------------------------------------|
| Min cell size        | 16 px | Smaller cells lack statistical samples |
| Invalid threshold    | 30%   | High clip fraction = star-contaminated |
| Grid filter          | 3x3 median | Smooth grid, reject cell-scale artifacts |
| Interpolation        | Catmull-Rom bicubic | C1-continuous, matches SExtractor/SEP |

---

## Comparison with Professional Tools

| Feature | SExtractor/SEP | photutils | rustafits |
|---------|---------------|-----------|-----------|
| Mode estimator | 2.5m - 1.5mu | Pluggable (SExtractor, biweight, MMM) | 2.5m - 1.5mu |
| Asymmetry fallback | Yes (0.3sigma) | Depends on estimator | Yes (0.3sigma) |
| Noise estimator | Clipped sigma | MAD or biweight scale | MAD (1.4826 * MAD) |
| Grid filtering | 3x3 median | Configurable median | 3x3 median |
| Interpolation | Bicubic spline | Bicubic spline (scipy zoom) | Catmull-Rom bicubic |
| Source masking | No (single pass) | Manual mask input | Iterative (auto-mask) |
| Parallelization | No | No (Python) | Yes (rayon per-row) |
