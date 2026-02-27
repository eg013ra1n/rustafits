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
| Mode Estimate       |  mode = 2.5 * median - 1.5 * mean   (SExtractor formula)
| (SExtractor)        |  noise = max(1.4826 * MAD, 0.001)
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

### MAD to Sigma Conversion

The Median Absolute Deviation (MAD) is a robust scale estimator:

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
| Clip threshold   | 3.0 sigma | Rejects >99.7% outliers per round    |
| MAD scale factor | 1.4826  | Gaussian consistency factor          |
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
| Bilinear Interpolation    |  Map grid -> full resolution background map
|                           |
|                           |  For pixel (x, y):
|                           |    fx = (x - cell_size/2) / cell_size
|                           |    fy = (y - cell_size/2) / cell_size
|                           |    bg(x,y) = bilinear(grid, fx, fy)
+---------------------------+
       |
       v
  background_map[y * width + x],
  global_background = median(grid),
  global_noise = median(cell_sigmas)
```

### Invalid Cell Detection

A cell is marked invalid when sigma-clipping removes more than 30% of its pixels.
This indicates heavy star contamination — the remaining "background" pixels are
unreliable. The nearest valid neighbor provides a better estimate.

### Bilinear Interpolation

Grid cell values represent the background at cell centers. Between centers, bilinear
interpolation produces a smooth, continuous background map:

```
bg(x,y) = v00*(1-tx)*(1-ty) + v10*tx*(1-ty) + v01*(1-tx)*ty + v11*tx*ty
```

where `tx, ty` are fractional positions between cell centers, and `v00..v11` are the
four surrounding cell values.

### Constants

| Parameter            | Value | Rationale                              |
|----------------------|-------|----------------------------------------|
| Min cell size        | 16 px | Smaller cells lack statistical samples |
| Invalid threshold    | 30%   | High clip fraction = star-contaminated |
| Grid filter          | 3x3 median | Smooth grid, reject cell-scale artifacts |
