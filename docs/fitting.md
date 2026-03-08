# PSF Fitting (Levenberg-Marquardt)

Non-linear least-squares fitting of 2D Gaussian and Moffat PSF models using the
Levenberg-Marquardt algorithm. Used for accurate FWHM and eccentricity measurement.

Two PSF models are available:

| Model | Parameters | Use Case |
|-------|------------|----------|
| **Gaussian** | 7 (B, A, x0, y0, sigma_x, sigma_y, theta) | Fast, good for well-sampled stars |
| **Moffat** (default) | 8 free / 7 fixed-beta | Better wing modeling, reports beta shape parameter |

Moffat fitting is enabled by default. See [Moffat section](#2d-elliptical-moffat-fit) below.

## 2D Elliptical Gaussian Fit (Primary — Used for Stars)

### Model

```
f(x, y) = B + A * exp(-Q(x,y) / 2)
```

where Q is the rotated quadratic form:

```
u = (x - x0) * cos(theta) + (y - y0) * sin(theta)
v = -(x - x0) * sin(theta) + (y - y0) * cos(theta)

Q(x, y) = u^2 / sigma_x^2  +  v^2 / sigma_y^2
```

**7 Parameters:** `[B, A, x0, y0, sigma_x, sigma_y, theta]`

- `B` — background pedestal
- `A` — peak amplitude above background
- `x0, y0` — centroid position
- `sigma_x, sigma_y` — Gaussian widths along major/minor axes
- `theta` — rotation angle of the ellipse

### Two-Stage Fitting Strategy

```
Input: star stamp (background-subtracted pixels + positions)
       |
       v
+------------------------------+
| Initial Parameter Estimates  |
|                              |
|   B0 = 0 (pre-subtracted)   |
|   A0 = max(stamp pixels)    |
|   x0, y0 = input centroid   |
|   sigma0 = sqrt(area/pi)*0.5|
|   theta0 = 0                |
+------------------------------+
       |
       v
+------------------------------+
| Stage 1: Axis-Aligned Fit   |  6 parameters (theta frozen at 0)
| Levenberg-Marquardt          |
|   max 50 iterations          |
|   convergence: 1e-6          |
+------------------------------+
       |
       v
+------------------------------+
| Check Ellipticity            |  e = |sigma_x - sigma_y| / max(sigma_x, sigma_y)
|                              |
|   if e > 0.1:               |
|     -> Stage 2               |
|   else:                     |
|     -> Done (theta = 0)     |
+------------------------------+
       |  (e > 0.1)
       v
+------------------------------+
| Stage 2: Full Elliptical Fit |  7 parameters (theta unfrozen)
| Levenberg-Marquardt          |  Initialized from Stage 1 result
|   max 50 iterations          |
+------------------------------+
       |
       v
  Gaussian2DResult { B, A, x0, y0, sigma_x, sigma_y, theta }
```

### Levenberg-Marquardt Algorithm

Hybrid between gradient descent (far from minimum) and Gauss-Newton (near minimum).

```
At each iteration k:

  1. Compute residuals:  r_i = observed_i - model_i

  2. Build Jacobian J (n_pixels x n_params):
       J[i][0] = df/dB     = 1
       J[i][1] = df/dA     = exp(-Q/2)
       J[i][2] = df/dx0    = A * exp(-Q/2) * (cos(th)*u/sx^2 - sin(th)*v/sy^2)
       J[i][3] = df/dy0    = A * exp(-Q/2) * (sin(th)*u/sx^2 + cos(th)*v/sy^2)
       J[i][4] = df/dsigx  = A * exp(-Q/2) * u^2 / sx^3
       J[i][5] = df/dsigy  = A * exp(-Q/2) * v^2 / sy^3
       J[i][6] = df/dtheta = A * exp(-Q/2) * u*v * (1/sy^2 - 1/sx^2)

  3. Normal equations:
       H = J^T * J                     (approximate Hessian)
       g = J^T * r                     (gradient)

  4. Damped solve:
       (H + lambda * diag(H)) * delta = g

  5. Trial step:
       params_trial = params + delta

  6. Evaluate gain ratio:
       rho = (cost_old - cost_new) / predicted_decrease

  7. Accept/reject:
       if rho > 0:  accept step, lambda <- lambda / 3
       if rho <= 0: reject step, lambda <- lambda * 2

  8. Check convergence:
       ||delta|| / ||params|| < 1e-6  ->  converged
```

### Linear System Solver

The damped normal equations are solved via **Cholesky decomposition** (the matrix
`H + lambda*diag(H)` is symmetric positive-definite). All arithmetic is in f64 for
numerical stability.

### Input Filtering (OSC)

For OSC (Bayer) images, only green CFA pixels are provided as input samples.
The caller filters the stamp using a `green_mask` before building the
`PixelSample` vector. This eliminates ~6% PSF broadening from interpolated
R/B positions while retaining ~50% of pixels — more than sufficient for the
7-parameter fit. The fitting algorithm itself is unchanged.

### Validation

The fit result is rejected (returns `None`) if:
- `sigma_x < 0.3` or `sigma_y < 0.3` (sub-pixel — unphysical)
- `A <= 0` (no star detected)
- Fewer than 10 pixel samples (underdetermined)

### FWHM and Eccentricity from Fit

The LM optimizer can converge with sigma_x < sigma_y (the PSF is identical with
swapped axes and theta rotated by π/2). The metrics layer canonicalizes the output
so that **fwhm_x ≥ fwhm_y** and **theta points along the major axis**:

```
if sigma_x >= sigma_y:
    FWHM_x = 2.3548 * sigma_x        (major)
    FWHM_y = 2.3548 * sigma_y        (minor)
    theta  = theta                    (unchanged)
else:
    FWHM_x = 2.3548 * sigma_y        (major — was the y-axis)
    FWHM_y = 2.3548 * sigma_x        (minor — was the x-axis)
    theta  = theta + π/2             (rotated to match new major)

FWHM   = sqrt(FWHM_x * FWHM_y)      (geometric mean, invariant to swap)

eccentricity = sqrt(1 - FWHM_y^2 / FWHM_x^2)
```

The same canonicalization applies to Moffat fits (alpha_x/alpha_y instead of
sigma_x/sigma_y).

---

## 1D Gaussian Fit (Utility)

**Model:** `f(x) = B + A * exp(-(x - mu)^2 / (2 * sigma^2))`

**Parameters:** `[B, A, mu, sigma]`

Same Levenberg-Marquardt algorithm but with 4 parameters and a 4x4 system solved via
Gaussian elimination with partial pivoting. Used for profile slicing and diagnostics.

### Initial Estimates

```
B0    = (y[0] + y[n-1]) / 2          (average of endpoints)
A0    = max(y) - B0                   (peak above background)
mu0   = argmax(y)                     (peak position)
sigma0 = half-max-width / 2.3548      (from left/right half-maximum crossings)
```

---

## 2D Elliptical Moffat Fit

The Moffat profile models PSF wings more accurately than a Gaussian. Real stellar PSFs
have power-law wings from atmospheric scattering and optical diffraction that a Gaussian
cannot capture. The Moffat function interpolates between Gaussian (beta → ∞) and
Lorentzian (beta = 1) profiles.

Moffat fitting is **always enabled** in the pipeline.

### Model

```
M(x, y) = B + A × (1 + Q(x,y))^(-β)
```

where Q is the rotated quadratic form:

```
u = (x - x0) × cos(θ) + (y - y0) × sin(θ)
v = -(x - x0) × sin(θ) + (y - y0) × cos(θ)

Q(x, y) = (u / α_x)² + (v / α_y)²
```

**8 Parameters (free beta):** `[B, A, x0, y0, alpha_x, alpha_y, theta, beta]`

- `B` — background pedestal
- `A` — peak amplitude above background
- `x0, y0` — centroid position
- `alpha_x, alpha_y` — Moffat scale parameters along major/minor axes
- `theta` — rotation angle of the ellipse
- `beta` — wing slope exponent (higher = steeper = more Gaussian-like)

**7 Parameters (fixed beta):** Same as above but `beta` is held constant during
optimization. The fixed beta value is automatically derived from the calibration
pass (sigma-clipped median of free-beta fits on bright calibration stars). This
removes one degree of freedom, improving FWHM stability when beta and axis ratio
trade off against each other.

### FWHM from Moffat Parameters

```
FWHM = 2α × √(2^(1/β) - 1)
```

This is the analytic half-maximum width of the Moffat profile. For beta → ∞,
the factor `√(2^(1/β) - 1)` approaches `√(ln 2)` and recovers the Gaussian FWHM.

### Fitting Strategy

```
Input: star stamp (background-subtracted pixels + positions)
       |
       v
+------------------------------+
| Initial Parameter Estimates  |
|                              |
|   B0 = 0 (pre-subtracted)   |
|   A0 = max(stamp pixels)    |
|   x0, y0 = input centroid   |
|   α0 = σ_gauss × FWHM / (2 × √(2^(1/β₀) - 1))
|   β0 = fixed_beta or 3.0    |
|   theta0 = 0                |
+------------------------------+
       |
       v
+------------------------------+
| Levenberg-Marquardt          |  np = 7 (fixed beta) or 8 (free beta)
|   max 50 iterations          |  Same LM engine as Gaussian
|   convergence: 1e-6          |  Cholesky solver (np × np system)
+------------------------------+
       |
       v
  Moffat2DResult { b, a, x0, y0, alpha_x, alpha_y, theta, beta, converged }
```

### Jacobian (8-parameter, free beta)

```
At each pixel (x, y):

  base = 1 + Q(x,y)
  power = base^(-β)
  dpower/dQ = -β × base^(-β - 1)

  J[0] = 1                        dM/dB
  J[1] = power                    dM/dA
  J[2] = A × dpower/dQ × dQ/dx0  dM/dx0
  J[3] = A × dpower/dQ × dQ/dy0  dM/dy0
  J[4] = A × dpower/dQ × (-2u²/α_x³)    dM/dα_x
  J[5] = A × dpower/dQ × (-2v²/α_y³)    dM/dα_y
  J[6] = A × dpower/dQ × dQ/dθ   dM/dθ
  J[7] = -A × power × ln(base)   dM/dβ
```

For fixed beta (7-parameter): `J[7]` (dM/dβ) is omitted. The Cholesky solver
operates on a 7×7 system instead of 8×8.

### Beta Interpretation

| Beta | Wing Shape | Equivalent |
|------|-----------|------------|
| 1.0 | Very broad wings | Lorentzian |
| 2.0 | Moderate wings | Typical poor seeing |
| 3.0 | Default initial guess | Typical good optics |
| 4.0 | Gentle wings | Well-corrected optics, common fixed value |
| 4.765 | Gaussian-equivalent | Exactly matches Gaussian shape |
| > 6 | Very steep falloff | Near-Gaussian, unusual |

### Validation

The fit result is rejected (returns `None`) if:
- `alpha_x < 0.3` or `alpha_y < 0.3` (sub-pixel — unphysical)
- `A <= 0` (no star detected)
- `beta < 1.0` or `beta > 20.0` (unphysical wing parameter)
- Fewer than 12 pixel samples (underdetermined for 8-parameter system)

### FWHM and Eccentricity from Moffat Fit

The same axis canonicalization as the Gaussian fit applies:

```
if alpha_x >= alpha_y:
    FWHM_x = 2 × alpha_x × √(2^(1/β) - 1)    (major)
    FWHM_y = 2 × alpha_y × √(2^(1/β) - 1)    (minor)
    theta  = theta                              (unchanged)
else:
    FWHM_x = 2 × alpha_y × √(2^(1/β) - 1)    (major — was y-axis)
    FWHM_y = 2 × alpha_x × √(2^(1/β) - 1)    (minor — was x-axis)
    theta  = theta + π/2                        (rotated to match)

FWHM   = √(FWHM_x × FWHM_y)     (geometric mean)

eccentricity = √(1 - FWHM_y² / FWHM_x²)
```

### Fitting Fallback Chain

The pipeline attempts PSF models in order of decreasing accuracy. Each star is
measured by the first method that succeeds:

```
Moffat fit (primary)
  |
  |-- success --> StarMetrics (FreeMoffat or FixedMoffat)
  |
  |-- fail (non-convergence / unphysical params)
  v
Gaussian fit (fallback)
  |
  |-- success --> StarMetrics (Gaussian)
  |
  |-- fail
  v
Windowed moments (last resort)
  --> StarMetrics (Moments) — flagged as lowest accuracy
```

The `FitMethod` enum on each `StarMetrics` records which method produced the result.

| Aspect | Moffat (primary) | Gaussian (fallback) | Moments (last resort) |
|--------|-------------------|---------------------|-----------------------|
| Wing accuracy | Accurate power-law wings | Underestimates wings | N/A |
| Parameters | 8 (free) or 7 (fixed beta) | 7 | N/A |
| Min samples | 12 | 10 | 1 |
| Reports beta | Yes | No | No |
| Accuracy | Highest | Good | Lowest |

---

## Two-Pass Calibration

The analysis pipeline uses a two-pass strategy to derive per-field PSF parameters
before measuring all stars.

### Pass 1: Calibration (free-beta Moffat)

1. Detect stars with an initial FWHM estimate.
2. Select up to **100 bright calibration stars** filtered by:
   - Eccentricity < 0.5 (reject elongated/blended sources)
   - Area >= 5 pixels (reject hot pixels and noise spikes)
3. Fit each calibration star with a **free-beta Moffat** (8 parameters).
4. Derive field-wide PSF parameters via **sigma-clipped median** of the
   calibration results:
   - `field_beta` -- the characteristic Moffat beta for the image
   - `field_fwhm` -- used to size source masks and refine the detection kernel

If fewer than 3 calibration stars yield valid beta values, `field_beta` falls
back to a simple median (or `None` if no betas are available), and the pipeline
proceeds without a fixed beta constraint.

### Pass 2: Measurement (fixed-beta Moffat)

1. Re-estimate background with source masks sized by `field_fwhm`.
2. Re-detect stars with the refined FWHM kernel.
3. Measure every detected star using the fallback chain:
   - **Fixed-beta Moffat** (7 params, using `field_beta`) -- primary
   - Gaussian -- fallback
   - Moments -- last resort

When `field_beta` is `None` (calibration did not produce one), pass 2 falls
back to free-beta Moffat as the primary fitter.

### FitMethod Enum

Each `StarMetrics` carries a `FitMethod` value indicating which model produced
the measurement:

| Variant | Description |
|---------|-------------|
| `FreeMoffat` | Free-beta Moffat (8 params) -- highest accuracy |
| `FixedMoffat` | Fixed-beta Moffat (7 params) -- field median beta |
| `Gaussian` | Gaussian fallback (7 params) |
| `Moments` | Windowed moments -- lowest accuracy, flagged unreliable |

---

## Constants

| Parameter           | Value  | Rationale                                     |
|---------------------|--------|-----------------------------------------------|
| Max iterations      | 50     | Sufficient for well-conditioned star profiles |
| Convergence tol     | 1e-6   | Relative parameter change threshold           |
| Min sigma (Gaussian)| 0.3 px | Below this, fit is unphysical (noise artifact)|
| Min alpha (Moffat)  | 0.3 px | Below this, fit is unphysical (noise artifact)|
| Max beta            | 20.0   | Above this, effectively Gaussian              |
| Min beta            | 1.0    | Below this, unphysical                        |
| Ellipticity trigger | 0.1    | Below this, rotation angle adds no value      |
| Min pixel samples   | 10/12  | Gaussian / Moffat minimum for overdetermined system |
| Arithmetic          | f64    | Avoid float32 precision loss in normal eqns   |
