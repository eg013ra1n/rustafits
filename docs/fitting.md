# Gaussian Fitting (Levenberg-Marquardt)

Non-linear least-squares fitting of 1D and 2D Gaussian models using the
Levenberg-Marquardt algorithm. Used for accurate FWHM and eccentricity measurement.

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
|   max 30 iterations          |
|   convergence: 1e-7          |
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
|   max 30 iterations          |
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
       ||delta|| / ||params|| < 1e-7  ->  converged
```

### Linear System Solver

The damped normal equations are solved via **Cholesky decomposition** (the matrix
`H + lambda*diag(H)` is symmetric positive-definite). All arithmetic is in f64 for
numerical stability.

### Validation

The fit result is rejected (returns `None`) if:
- `sigma_x < 0.3` or `sigma_y < 0.3` (sub-pixel — unphysical)
- `A <= 0` (no star detected)
- Fewer than 10 pixel samples (underdetermined)

### FWHM and Eccentricity from Fit

```
FWHM_x = 2.3548 * sigma_x
FWHM_y = 2.3548 * sigma_y
FWHM   = sqrt(FWHM_x * FWHM_y)     (geometric mean)

min_s  = min(sigma_x, sigma_y)
max_s  = max(sigma_x, sigma_y)
eccentricity = sqrt(1 - min_s^2 / max_s^2)
```

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

## Constants

| Parameter           | Value  | Rationale                                     |
|---------------------|--------|-----------------------------------------------|
| Max iterations      | 30     | Sufficient for well-conditioned star profiles |
| Convergence tol     | 1e-7   | Relative parameter change threshold           |
| Min sigma           | 0.3 px | Below this, fit is unphysical (noise artifact)|
| Ellipticity trigger | 0.1    | Below this, rotation angle adds no value      |
| Min pixel samples   | 10     | Need overdetermined system for robust fit     |
| Arithmetic          | f64    | Avoid float32 precision loss in normal eqns   |
