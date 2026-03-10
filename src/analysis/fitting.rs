/// Levenberg-Marquardt 2D elliptical Gaussian fitting (7-param).
/// All internal computations in f64 for numerical stability.

/// Default LM solver parameters for calibration pass (full precision).
/// Used by the debug binary and tests; the library API passes values explicitly.
#[allow(dead_code)]
pub const CALIBRATION_MAX_ITER: usize = 50;
#[allow(dead_code)]
pub const CALIBRATION_CONV_TOL: f64 = 1e-6;
#[allow(dead_code)]
pub const CALIBRATION_MAX_REJECTS: usize = 5;

/// 2D Gaussian model: f(x,y) = B + A * exp(-0.5 * Q(x,y))
/// Q(x,y) = u^2/sx^2 + v^2/sy^2
/// u = (x-x0)*cos(theta) + (y-y0)*sin(theta)
/// v = -(x-x0)*sin(theta) + (y-y0)*cos(theta)
/// Params: [B, A, x0, y0, sigma_x, sigma_y, theta]
#[allow(dead_code)]
pub struct Gaussian2DResult {
    pub b: f64,
    pub a: f64,
    pub x0: f64,
    pub y0: f64,
    pub sigma_x: f64,
    pub sigma_y: f64,
    pub theta: f64,
    pub converged: bool,
}

/// Pixel coordinate + value for 2D fitting.
pub struct PixelSample {
    pub x: f64,
    pub y: f64,
    pub value: f64,
}

/// Fit 2D elliptical Gaussian with 7 parameters (including theta).
/// Caller provides per-axis sigma and theta initialisation from moments.
///
/// `max_iter`, `conv_tol`, `max_rejects` control LM solver convergence.
/// Use `CALIBRATION_*` constants for full-precision calibration fits.
pub fn fit_gaussian_2d(
    pixels: &[PixelSample],
    init_b: f64,
    init_a: f64,
    init_x0: f64,
    init_y0: f64,
    init_sigma_x: f64,
    init_sigma_y: f64,
    init_theta: f64,
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
) -> Option<Gaussian2DResult> {
    if pixels.len() < 10 {
        return None;
    }

    let mut params = [init_b, init_a, init_x0, init_y0, init_sigma_x, init_sigma_y, init_theta];
    let converged = lm_solve_2d(pixels, &mut params, true, max_iter, conv_tol, max_rejects);

    let sigma_x = params[4].abs();
    let sigma_y = params[5].abs();
    if sigma_x < 0.3 || sigma_y < 0.3 || params[1] <= 0.0 {
        return None;
    }

    // Normalize theta to [-π/2, π/2]
    let mut theta = params[6] % std::f64::consts::PI;
    if theta > std::f64::consts::FRAC_PI_2 {
        theta -= std::f64::consts::PI;
    } else if theta < -std::f64::consts::FRAC_PI_2 {
        theta += std::f64::consts::PI;
    }

    Some(Gaussian2DResult {
        b: params[0],
        a: params[1],
        x0: params[2],
        y0: params[3],
        sigma_x,
        sigma_y,
        theta,
        converged,
    })
}

/// LM solver for 2D Gaussian. If `fit_theta` is true, uses 7 params; else 6 (theta=0).
fn lm_solve_2d(pixels: &[PixelSample], params: &mut [f64], fit_theta: bool, max_iter: usize, conv_tol: f64, max_rejects: usize) -> bool {
    let np = if fit_theta { 7 } else { 6 };
    let mut lambda = 1e-3_f64;
    let mut nu = 2.0_f64;
    let mut best_cost = residual_cost_2d(pixels, params, fit_theta);
    let mut prev_cost: f64;
    let mut converged = false;
    let mut consecutive_rejects: usize = 0;

    // Scratch space for normal equations
    let mut jtj = vec![0.0_f64; np * np];
    let mut jtr = vec![0.0_f64; np];
    let mut j = vec![0.0_f64; np];
    let mut mat = vec![0.0_f64; np * np];
    let mut new_params = [0.0_f64; 7];

    for _ in 0..max_iter {
        // Zero out
        jtj.fill(0.0);
        jtr.fill(0.0);

        let theta = if fit_theta { params[6] } else { 0.0 };
        let (cos_t, sin_t) = (theta.cos(), theta.sin());
        let sx = params[4];
        let sy = params[5];
        let inv_sx2 = 1.0 / (sx * sx);
        let inv_sy2 = 1.0 / (sy * sy);

        for px in pixels {
            let dx = px.x - params[2];
            let dy = px.y - params[3];
            let u = dx * cos_t + dy * sin_t;
            let v = -dx * sin_t + dy * cos_t;
            let q = u * u * inv_sx2 + v * v * inv_sy2;
            let e = (-0.5 * q).exp();
            let model = params[0] + params[1] * e;
            let r = px.value - model;

            // Jacobian
            j[0] = 1.0; // dF/dB
            j[1] = e; // dF/dA
            // dF/dx0
            j[2] = params[1] * e * (cos_t * u * inv_sx2 - sin_t * v * inv_sy2);
            // dF/dy0
            j[3] = params[1] * e * (sin_t * u * inv_sx2 + cos_t * v * inv_sy2);
            // dF/dsigma_x
            j[4] = params[1] * e * u * u / (sx * sx * sx);
            // dF/dsigma_y
            j[5] = params[1] * e * v * v / (sy * sy * sy);
            if fit_theta {
                // dF/dtheta
                j[6] = params[1] * e * u * v * (inv_sy2 - inv_sx2);
            }

            for p in 0..np {
                jtr[p] += j[p] * r;
                for q in p..np {
                    jtj[p * np + q] += j[p] * j[q];
                }
            }
        }

        // Fill symmetric lower triangle
        for p in 0..np {
            for q in 0..p {
                jtj[p * np + q] = jtj[q * np + p];
            }
        }

        // Damped normal equations
        mat.copy_from_slice(&jtj);
        for p in 0..np {
            mat[p * np + p] += lambda * jtj[p * np + p].max(1e-12);
        }

        let delta = {
            let mut result = cholesky_solve(&mat, &jtr, np);
            if result.is_none() {
                for _ in 0..3 {
                    lambda *= 10.0;
                    for p in 0..np {
                        mat[p * np + p] = jtj[p * np + p] + lambda * jtj[p * np + p].max(1e-12);
                    }
                    result = cholesky_solve(&mat, &jtr, np);
                    if result.is_some() { break; }
                }
            }
            match result {
                Some(d) => d,
                None => break,
            }
        };

        new_params[..np].copy_from_slice(&params[..np]);
        for p in 0..np {
            new_params[p] += delta[p];
        }
        // Keep sigmas positive
        if new_params[4] <= 0.0 {
            new_params[4] = params[4] * 0.5;
        }
        if new_params[5] <= 0.0 {
            new_params[5] = params[5] * 0.5;
        }

        // Reject step if sigma too small or centroid drifts too far
        if new_params[4] < 0.5 || new_params[5] < 0.5
            || {
                let cdx = new_params[2] - params[2];
                let cdy = new_params[3] - params[3];
                cdx * cdx + cdy * cdy > 4.0
            }
        {
            lambda *= nu;
            nu *= 2.0;
            consecutive_rejects += 1;
            if consecutive_rejects >= max_rejects {
                break;
            }
            continue;
        }

        let new_cost = residual_cost_2d(pixels, &new_params[..np], fit_theta);

        // Nielsen gain ratio
        let predicted: f64 = delta
            .iter()
            .enumerate()
            .map(|(i, d)| d * (lambda * jtj[i * np + i].max(1e-12) * d + jtr[i]))
            .sum();

        prev_cost = best_cost;
        let mut step_accepted = false;
        if predicted > 0.0 {
            let rho = (best_cost - new_cost) / predicted;
            if rho > 0.0 {
                params[..np].copy_from_slice(&new_params[..np]);
                best_cost = new_cost;
                lambda *= (1.0_f64 / 3.0).max(1.0 - (2.0 * rho - 1.0).powi(3));
                nu = 2.0;
                step_accepted = true;
                consecutive_rejects = 0;
            } else {
                lambda *= nu;
                nu *= 2.0;
                consecutive_rejects += 1;
                if consecutive_rejects >= max_rejects {
                    break;
                }
            }
        } else {
            lambda *= nu;
            nu *= 2.0;
            consecutive_rejects += 1;
            if consecutive_rejects >= max_rejects {
                break;
            }
        }

        // Dual convergence criterion
        let param_norm = params[..np].iter().map(|p| p * p).sum::<f64>().sqrt();
        let delta_norm = delta.iter().map(|d| d * d).sum::<f64>().sqrt();
        let param_converged = delta_norm / param_norm.max(1e-12) < conv_tol;
        let residual_converged = if step_accepted && best_cost > 0.0 {
            (prev_cost - best_cost).abs() / best_cost < 1e-4
        } else {
            false
        };
        if param_converged || residual_converged {
            converged = true;
            break;
        }
    }

    converged
}

fn residual_cost_2d(pixels: &[PixelSample], params: &[f64], fit_theta: bool) -> f64 {
    let theta = if fit_theta { params[6] } else { 0.0 };
    let (cos_t, sin_t) = (theta.cos(), theta.sin());
    let inv_sx2 = 1.0 / (params[4] * params[4]);
    let inv_sy2 = 1.0 / (params[5] * params[5]);

    pixels
        .iter()
        .map(|px| {
            let dx = px.x - params[2];
            let dy = px.y - params[3];
            let u = dx * cos_t + dy * sin_t;
            let v = -dx * sin_t + dy * cos_t;
            let q = u * u * inv_sx2 + v * v * inv_sy2;
            let model = params[0] + params[1] * (-0.5 * q).exp();
            let r = px.value - model;
            r * r
        })
        .sum()
}

/// Cholesky decomposition solver for symmetric positive-definite system.
/// Matrix stored as flat array, row-major, size np×np.
fn cholesky_solve(mat: &[f64], rhs: &[f64], np: usize) -> Option<Vec<f64>> {
    // Cholesky: A = L * L^T
    let mut l = vec![0.0_f64; np * np];

    for i in 0..np {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l[i * np + k] * l[j * np + k];
            }
            if i == j {
                let diag = mat[i * np + i] - sum;
                if diag <= 0.0 {
                    return None; // Not positive definite
                }
                l[i * np + j] = diag.sqrt();
            } else {
                l[i * np + j] = (mat[i * np + j] - sum) / l[j * np + j];
            }
        }
    }

    // Solve L * y = rhs (forward substitution)
    let mut y = vec![0.0_f64; np];
    for i in 0..np {
        let mut sum = 0.0;
        for j in 0..i {
            sum += l[i * np + j] * y[j];
        }
        y[i] = (rhs[i] - sum) / l[i * np + i];
    }

    // Solve L^T * x = y (back substitution)
    let mut x = vec![0.0_f64; np];
    for i in (0..np).rev() {
        let mut sum = 0.0;
        for j in (i + 1)..np {
            sum += l[j * np + i] * x[j]; // L^T[i][j] = L[j][i]
        }
        x[i] = (y[i] - sum) / l[i * np + i];
    }

    Some(x)
}

// ── Moffat PSF fitting ──────────────────────────────────────────────────────
//
// Elliptical Moffat: M(x,y) = B + A × (1 + Q(x,y))^(-β)
// Q(x,y) = (u/α_x)² + (v/α_y)²   (rotated coordinates)
// FWHM = 2α√(2^(1/β) - 1)
//
// 8 parameters: [B, A, x0, y0, alpha_x, alpha_y, theta, beta]

/// Result of 2D elliptical Moffat fit.
#[allow(dead_code)]
pub struct Moffat2DResult {
    pub b: f64,
    pub a: f64,
    pub x0: f64,
    pub y0: f64,
    pub alpha_x: f64,
    pub alpha_y: f64,
    pub theta: f64,
    pub beta: f64,
    pub converged: bool,
}

impl Moffat2DResult {
    /// Compute FWHM along each axis from Moffat parameters.
    /// FWHM = 2α√(2^(1/β) - 1)
    pub fn fwhm_x(&self) -> f64 {
        2.0 * self.alpha_x * (2.0_f64.powf(1.0 / self.beta) - 1.0).sqrt()
    }
    pub fn fwhm_y(&self) -> f64 {
        2.0 * self.alpha_y * (2.0_f64.powf(1.0 / self.beta) - 1.0).sqrt()
    }
}

/// Fit 2D elliptical Moffat with 8 parameters (free beta).
///
/// Initial alpha values are derived from Gaussian sigma: α = σ × √(2^(1/β₀) - 1) / √(2 ln 2)
/// where β₀ = 3.0 (typical for well-corrected optics).
pub fn fit_moffat_2d(
    pixels: &[PixelSample],
    init_b: f64,
    init_a: f64,
    init_x0: f64,
    init_y0: f64,
    init_sigma_x: f64,
    init_sigma_y: f64,
    init_theta: f64,
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
) -> Option<Moffat2DResult> {
    fit_moffat_2d_impl(pixels, init_b, init_a, init_x0, init_y0, init_sigma_x, init_sigma_y, init_theta, None, max_iter, conv_tol, max_rejects)
}

/// Fit 2D elliptical Moffat with fixed beta (7 free parameters).
///
/// Beta is held constant during optimization, removing the beta/axis-ratio
/// tradeoff that can destabilize FWHM. Common choice: beta=4 for
/// well-corrected optics.
pub fn fit_moffat_2d_fixed_beta(
    pixels: &[PixelSample],
    init_b: f64,
    init_a: f64,
    init_x0: f64,
    init_y0: f64,
    init_sigma_x: f64,
    init_sigma_y: f64,
    init_theta: f64,
    beta: f64,
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
) -> Option<Moffat2DResult> {
    fit_moffat_2d_impl(pixels, init_b, init_a, init_x0, init_y0, init_sigma_x, init_sigma_y, init_theta, Some(beta), max_iter, conv_tol, max_rejects)
}

/// Evaluate the 2D Moffat model at a single point.
///
/// `params`: [B, A, x0, y0, alpha_x, alpha_y, theta, beta]
/// Returns the model value at (x, y).
#[cfg(test)]
pub fn evaluate_moffat_2d(params: &[f64], x: f64, y: f64) -> f64 {
    let (cos_t, sin_t) = (params[6].cos(), params[6].sin());
    let dx = x - params[2];
    let dy = y - params[3];
    let u = dx * cos_t + dy * sin_t;
    let v = -dx * sin_t + dy * cos_t;
    let ax = params[4];
    let ay = params[5];
    let q = (u / ax).powi(2) + (v / ay).powi(2);
    params[0] + params[1] * (1.0 + q).powf(-params[7])
}

fn fit_moffat_2d_impl(
    pixels: &[PixelSample],
    init_b: f64,
    init_a: f64,
    init_x0: f64,
    init_y0: f64,
    init_sigma_x: f64,
    init_sigma_y: f64,
    init_theta: f64,
    fixed_beta: Option<f64>,
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
) -> Option<Moffat2DResult> {
    if pixels.len() < 12 {
        return None;
    }

    // Convert Gaussian sigma to Moffat alpha
    let beta_init = fixed_beta.unwrap_or(3.0);
    let moffat_scale = (2.0_f64.powf(1.0 / beta_init) - 1.0).sqrt();
    let alpha_x = init_sigma_x * 2.3548 / (2.0 * moffat_scale);
    let alpha_y = init_sigma_y * 2.3548 / (2.0 * moffat_scale);

    let mut params = [init_b, init_a, init_x0, init_y0, alpha_x, alpha_y, init_theta, beta_init];
    let converged = lm_solve_moffat_impl(pixels, &mut params, fixed_beta, max_iter, conv_tol, max_rejects);

    let alpha_x = params[4].abs();
    let alpha_y = params[5].abs();
    let beta = if fixed_beta.is_some() { beta_init } else { params[7] };

    // Reject non-physical results
    if alpha_x < 0.3 || alpha_y < 0.3 || params[1] <= 0.0
        || beta < 1.0 || beta > 20.0
    {
        return None;
    }

    // Normalize theta to [-π/2, π/2]
    let mut theta = params[6] % std::f64::consts::PI;
    if theta > std::f64::consts::FRAC_PI_2 {
        theta -= std::f64::consts::PI;
    } else if theta < -std::f64::consts::FRAC_PI_2 {
        theta += std::f64::consts::PI;
    }

    Some(Moffat2DResult {
        b: params[0],
        a: params[1],
        x0: params[2],
        y0: params[3],
        alpha_x,
        alpha_y,
        theta,
        beta,
        converged,
    })
}

/// LM solver for 2D Moffat. When `fixed_beta` is `Some(b)`, beta is held constant
/// and only 7 parameters are optimized. When `None`, all 8 parameters are free.
fn lm_solve_moffat_impl(pixels: &[PixelSample], params: &mut [f64], fixed_beta: Option<f64>, max_iter: usize, conv_tol: f64, max_rejects: usize) -> bool {
    let np = if fixed_beta.is_some() { 7 } else { 8 };
    let mut lambda = 1e-3_f64;
    let mut nu = 2.0_f64;
    let mut best_cost = residual_cost_moffat(pixels, params);
    let mut prev_cost: f64;
    let mut converged = false;
    let mut consecutive_rejects: usize = 0;

    let mut jtj = vec![0.0_f64; np * np];
    let mut jtr = vec![0.0_f64; np];
    let mut j = vec![0.0_f64; np];
    let mut mat = vec![0.0_f64; np * np];
    let mut new_params = [0.0_f64; 8];

    for _ in 0..max_iter {
        jtj.fill(0.0);
        jtr.fill(0.0);

        let theta = params[6];
        let (cos_t, sin_t) = (theta.cos(), theta.sin());
        let ax = params[4];
        let ay = params[5];
        let beta = fixed_beta.unwrap_or(params[7]);
        let inv_ax2 = 1.0 / (ax * ax);
        let inv_ay2 = 1.0 / (ay * ay);

        for px in pixels {
            let dx = px.x - params[2];
            let dy = px.y - params[3];
            let u = dx * cos_t + dy * sin_t;
            let v = -dx * sin_t + dy * cos_t;
            let q = u * u * inv_ax2 + v * v * inv_ay2;
            let base = 1.0 + q;
            let power = base.powf(-beta);
            let model = params[0] + params[1] * power;
            let r = px.value - model;

            // Jacobian of M(x,y) = B + A × (1+Q)^(-β)
            j[0] = 1.0; // dM/dB
            j[1] = power; // dM/dA

            let dpower_dq = -beta * base.powf(-beta - 1.0);
            let a_dpq = params[1] * dpower_dq;

            let dq_dx0 = -2.0 * (cos_t * u * inv_ax2 - sin_t * v * inv_ay2);
            j[2] = a_dpq * dq_dx0;

            let dq_dy0 = -2.0 * (sin_t * u * inv_ax2 + cos_t * v * inv_ay2);
            j[3] = a_dpq * dq_dy0;

            j[4] = a_dpq * (-2.0 * u * u / (ax * ax * ax));
            j[5] = a_dpq * (-2.0 * v * v / (ay * ay * ay));

            let dq_dtheta = 2.0 * u * v * (inv_ay2 - inv_ax2);
            j[6] = a_dpq * dq_dtheta;

            // Only include dM/dbeta when beta is free
            if fixed_beta.is_none() {
                j[7] = params[1] * power * (-base.ln());
            }

            for p in 0..np {
                jtr[p] += j[p] * r;
                for q in p..np {
                    jtj[p * np + q] += j[p] * j[q];
                }
            }
        }

        // Fill symmetric lower triangle
        for p in 0..np {
            for q in 0..p {
                jtj[p * np + q] = jtj[q * np + p];
            }
        }

        // Damped normal equations
        mat[..np * np].copy_from_slice(&jtj);
        for p in 0..np {
            mat[p * np + p] += lambda * jtj[p * np + p].max(1e-12);
        }

        let delta = {
            let mut result = cholesky_solve(&mat[..np * np], &jtr, np);
            if result.is_none() {
                for _ in 0..3 {
                    lambda *= 10.0;
                    for p in 0..np {
                        mat[p * np + p] = jtj[p * np + p] + lambda * jtj[p * np + p].max(1e-12);
                    }
                    result = cholesky_solve(&mat[..np * np], &jtr, np);
                    if result.is_some() { break; }
                }
            }
            match result {
                Some(d) => d,
                None => break,
            }
        };

        new_params[..8].copy_from_slice(&params[..8]);
        for p in 0..np {
            new_params[p] += delta[p];
        }
        // Keep alphas positive
        if new_params[4] <= 0.0 { new_params[4] = params[4] * 0.5; }
        if new_params[5] <= 0.0 { new_params[5] = params[5] * 0.5; }
        // Keep beta in reasonable range (only when free)
        if fixed_beta.is_none() {
            if new_params[7] < 1.5 { new_params[7] = 1.5; }
            if new_params[7] > 20.0 { new_params[7] = 20.0; }
        }

        // Reject step if alpha too small or centroid drifts too far
        if new_params[4] < 0.5 || new_params[5] < 0.5
            || {
                let cdx = new_params[2] - params[2];
                let cdy = new_params[3] - params[3];
                cdx * cdx + cdy * cdy > 4.0
            }
        {
            lambda *= nu;
            nu *= 2.0;
            consecutive_rejects += 1;
            if consecutive_rejects >= max_rejects {
                break;
            }
            continue;
        }

        let new_cost = residual_cost_moffat(pixels, &new_params);

        let predicted: f64 = delta
            .iter()
            .enumerate()
            .map(|(i, d)| d * (lambda * jtj[i * np + i].max(1e-12) * d + jtr[i]))
            .sum();

        prev_cost = best_cost;
        let mut step_accepted = false;
        if predicted > 0.0 {
            let rho = (best_cost - new_cost) / predicted;
            if rho > 0.0 {
                params[..8].copy_from_slice(&new_params[..8]);
                best_cost = new_cost;
                lambda *= (1.0_f64 / 3.0).max(1.0 - (2.0 * rho - 1.0).powi(3));
                nu = 2.0;
                step_accepted = true;
                consecutive_rejects = 0;
            } else {
                lambda *= nu;
                nu *= 2.0;
                consecutive_rejects += 1;
                if consecutive_rejects >= max_rejects {
                    break;
                }
            }
        } else {
            lambda *= nu;
            nu *= 2.0;
            consecutive_rejects += 1;
            if consecutive_rejects >= max_rejects {
                break;
            }
        }

        // Dual convergence criterion
        let param_norm = params[..np].iter().map(|p| p * p).sum::<f64>().sqrt();
        let delta_norm = delta.iter().map(|d| d * d).sum::<f64>().sqrt();
        let param_converged = delta_norm / param_norm.max(1e-12) < conv_tol;
        let residual_converged = if step_accepted && best_cost > 0.0 {
            (prev_cost - best_cost).abs() / best_cost < 1e-4
        } else {
            false
        };
        if param_converged || residual_converged {
            converged = true;
            break;
        }
    }

    converged
}

fn residual_cost_moffat(pixels: &[PixelSample], params: &[f64]) -> f64 {
    let theta = params[6];
    let (cos_t, sin_t) = (theta.cos(), theta.sin());
    let inv_ax2 = 1.0 / (params[4] * params[4]);
    let inv_ay2 = 1.0 / (params[5] * params[5]);
    let beta = params[7];

    pixels
        .iter()
        .map(|px| {
            let dx = px.x - params[2];
            let dy = px.y - params[3];
            let u = dx * cos_t + dy * sin_t;
            let v = -dx * sin_t + dy * cos_t;
            let q = u * u * inv_ax2 + v * v * inv_ay2;
            let model = params[0] + params[1] * (1.0 + q).powf(-beta);
            let r = px.value - model;
            r * r
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_2d_gaussian_fit_isotropic() {
        // Generate isotropic 2D Gaussian: B=100, A=5000, x0=10, y0=10, sigma=3
        let size = 21;
        let mut pixels = Vec::new();
        for y in 0..size {
            for x in 0..size {
                let dx = x as f64 - 10.0;
                let dy = y as f64 - 10.0;
                let v = 100.0 + 5000.0 * (-0.5 * (dx * dx + dy * dy) / 9.0).exp();
                pixels.push(PixelSample {
                    x: x as f64,
                    y: y as f64,
                    value: v,
                });
            }
        }

        let result = fit_gaussian_2d(&pixels, 100.0, 5000.0, 10.0, 10.0, 3.0, 3.0, 0.0, CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS).unwrap();
        assert!(result.converged);
        assert!((result.sigma_x - 3.0).abs() < 0.1, "sx: {}", result.sigma_x);
        assert!((result.sigma_y - 3.0).abs() < 0.1, "sy: {}", result.sigma_y);
        assert!((result.x0 - 10.0).abs() < 0.05, "x0: {}", result.x0);
        assert!((result.y0 - 10.0).abs() < 0.05, "y0: {}", result.y0);
    }

    #[test]
    fn test_2d_gaussian_fit_elliptical() {
        // Elongated: sigma_x=2, sigma_y=5, theta=0.3 rad
        let size = 25;
        let theta = 0.3_f64;
        let (ct, st) = (theta.cos(), theta.sin());
        let mut pixels = Vec::new();
        for y in 0..size {
            for x in 0..size {
                let dx = x as f64 - 12.0;
                let dy = y as f64 - 12.0;
                let u = dx * ct + dy * st;
                let v = -dx * st + dy * ct;
                let q = u * u / 4.0 + v * v / 25.0;
                let val = 50.0 + 8000.0 * (-0.5 * q).exp();
                pixels.push(PixelSample {
                    x: x as f64,
                    y: y as f64,
                    value: val,
                });
            }
        }

        let result = fit_gaussian_2d(&pixels, 50.0, 8000.0, 12.0, 12.0, 3.0, 3.0, 0.0, CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS).unwrap();
        assert!(result.converged);
        let min_s = result.sigma_x.min(result.sigma_y);
        let max_s = result.sigma_x.max(result.sigma_y);
        assert!((min_s - 2.0).abs() < 0.2, "min_sigma: {}", min_s);
        assert!((max_s - 5.0).abs() < 0.3, "max_sigma: {}", max_s);
    }

    #[test]
    fn test_2d_gaussian_fit_elliptical_45deg() {
        // Elongated star at 45 degrees — the exact failure case for two-stage fitting.
        // sigma_x=2, sigma_y=5, theta=π/4
        let size = 25;
        let theta = std::f64::consts::FRAC_PI_4;
        let (ct, st) = (theta.cos(), theta.sin());
        let mut pixels = Vec::new();
        for y in 0..size {
            for x in 0..size {
                let dx = x as f64 - 12.0;
                let dy = y as f64 - 12.0;
                let u = dx * ct + dy * st;
                let v = -dx * st + dy * ct;
                let q = u * u / 4.0 + v * v / 25.0;
                let val = 50.0 + 8000.0 * (-0.5 * q).exp();
                pixels.push(PixelSample {
                    x: x as f64,
                    y: y as f64,
                    value: val,
                });
            }
        }

        // Provide moments-based init: theta=π/4, per-axis sigmas close to truth
        let result = fit_gaussian_2d(&pixels, 50.0, 8000.0, 12.0, 12.0, 2.0, 5.0, theta, CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS).unwrap();
        assert!(result.converged, "should converge");
        let min_s = result.sigma_x.min(result.sigma_y);
        let max_s = result.sigma_x.max(result.sigma_y);
        assert!(
            (min_s - 2.0).abs() < 0.2,
            "min_sigma: {} expected ~2.0",
            min_s
        );
        assert!(
            (max_s - 5.0).abs() < 0.3,
            "max_sigma: {} expected ~5.0",
            max_s
        );
        // Eccentricity: sqrt(1 - 4/25) ≈ 0.917
        let ecc = (1.0 - (min_s * min_s) / (max_s * max_s)).max(0.0).sqrt();
        assert!(
            (ecc - 0.917).abs() < 0.05,
            "eccentricity: {} expected ~0.917",
            ecc
        );
    }

    #[test]
    fn test_moffat_isotropic() {
        // Generate isotropic Moffat: B=100, A=5000, x0=15, y0=15, alpha=3, beta=3
        let size = 31;
        let alpha = 3.0_f64;
        let beta = 3.0_f64;
        let mut pixels = Vec::new();
        for y in 0..size {
            for x in 0..size {
                let dx = x as f64 - 15.0;
                let dy = y as f64 - 15.0;
                let r2 = dx * dx + dy * dy;
                let v = 100.0 + 5000.0 * (1.0 + r2 / (alpha * alpha)).powf(-beta);
                pixels.push(PixelSample {
                    x: x as f64,
                    y: y as f64,
                    value: v,
                });
            }
        }

        // Init from Gaussian sigma ≈ FWHM / 2.3548, where Moffat FWHM = 2α√(2^(1/β)-1)
        let moffat_fwhm = 2.0 * alpha * (2.0_f64.powf(1.0 / beta) - 1.0).sqrt();
        let init_sigma = moffat_fwhm / 2.3548;
        let result = fit_moffat_2d(&pixels, 100.0, 5000.0, 15.0, 15.0, init_sigma, init_sigma, 0.0, CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS).unwrap();

        assert!(result.converged, "should converge");
        assert!((result.alpha_x - alpha).abs() < 0.3, "alpha_x: {} expected ~{}", result.alpha_x, alpha);
        assert!((result.alpha_y - alpha).abs() < 0.3, "alpha_y: {} expected ~{}", result.alpha_y, alpha);
        assert!((result.beta - beta).abs() < 0.5, "beta: {} expected ~{}", result.beta, beta);
        assert!((result.x0 - 15.0).abs() < 0.1, "x0: {}", result.x0);
        assert!((result.y0 - 15.0).abs() < 0.1, "y0: {}", result.y0);
    }

    #[test]
    fn test_moffat_fwhm_vs_gaussian() {
        // A Moffat profile should yield ~5-15% smaller FWHM than a Gaussian fit
        // to the same data, since Moffat models the wings better.
        let size = 41;
        let alpha = 4.0_f64;
        let beta = 2.5_f64;
        let mut pixels = Vec::new();
        for y in 0..size {
            for x in 0..size {
                let dx = x as f64 - 20.0;
                let dy = y as f64 - 20.0;
                let r2 = dx * dx + dy * dy;
                let v = 50.0 + 8000.0 * (1.0 + r2 / (alpha * alpha)).powf(-beta);
                pixels.push(PixelSample { x: x as f64, y: y as f64, value: v });
            }
        }

        let true_fwhm = 2.0 * alpha * (2.0_f64.powf(1.0 / beta) - 1.0).sqrt();
        let init_sigma = true_fwhm / 2.3548;

        let gauss = fit_gaussian_2d(&pixels, 50.0, 8000.0, 20.0, 20.0, init_sigma, init_sigma, 0.0, CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS).unwrap();
        let moffat = fit_moffat_2d(&pixels, 50.0, 8000.0, 20.0, 20.0, init_sigma, init_sigma, 0.0, CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS).unwrap();

        let gauss_fwhm = 2.3548 * (gauss.sigma_x * gauss.sigma_y).sqrt();
        let moffat_fwhm = (moffat.fwhm_x() * moffat.fwhm_y()).sqrt();

        // Moffat should recover true FWHM more accurately
        assert!(
            (moffat_fwhm - true_fwhm).abs() < (gauss_fwhm - true_fwhm).abs(),
            "Moffat FWHM {:.3} should be closer to true {:.3} than Gaussian {:.3}",
            moffat_fwhm, true_fwhm, gauss_fwhm
        );
    }

    #[test]
    fn test_moffat_fixed_beta() {
        // Generate isotropic Moffat with beta=4, then fit with fixed beta=4.
        // Should converge and recover alpha accurately.
        let size = 31;
        let alpha = 3.0_f64;
        let beta = 4.0_f64;
        let mut pixels = Vec::new();
        for y in 0..size {
            for x in 0..size {
                let dx = x as f64 - 15.0;
                let dy = y as f64 - 15.0;
                let r2 = dx * dx + dy * dy;
                let v = 100.0 + 5000.0 * (1.0 + r2 / (alpha * alpha)).powf(-beta);
                pixels.push(PixelSample {
                    x: x as f64,
                    y: y as f64,
                    value: v,
                });
            }
        }

        let moffat_fwhm = 2.0 * alpha * (2.0_f64.powf(1.0 / beta) - 1.0).sqrt();
        let init_sigma = moffat_fwhm / 2.3548;

        let result = fit_moffat_2d_fixed_beta(
            &pixels, 100.0, 5000.0, 15.0, 15.0, init_sigma, init_sigma, 0.0, beta,
            CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS,
        ).unwrap();

        assert!(result.converged, "should converge");
        assert!((result.beta - beta).abs() < 0.01, "beta should be fixed at {}, got {}", beta, result.beta);
        assert!((result.alpha_x - alpha).abs() < 0.3, "alpha_x: {} expected ~{}", result.alpha_x, alpha);
        assert!((result.alpha_y - alpha).abs() < 0.3, "alpha_y: {} expected ~{}", result.alpha_y, alpha);
    }

    #[test]
    fn test_cholesky_identity() {
        // 3×3 identity system: I * x = [1, 2, 3]
        let mat = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let rhs = vec![1.0, 2.0, 3.0];
        let x = cholesky_solve(&mat, &rhs, 3).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 2.0).abs() < 1e-10);
        assert!((x[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_cholesky_retry_on_faint_star() {
        let mut samples = Vec::new();
        let sigma = 2.0_f64;
        let amplitude = 5.0_f64;
        let noise_bg = 100.0_f64;
        for y in -8..=8 {
            for x in -8..=8 {
                let r2 = (x * x + y * y) as f64;
                let val = noise_bg + amplitude * (-0.5 * r2 / (sigma * sigma)).exp();
                samples.push(PixelSample { x: x as f64, y: y as f64, value: val });
            }
        }
        let result = fit_gaussian_2d(&samples, noise_bg, amplitude, 0.0, 0.0, sigma, sigma, 0.0, CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS);
        assert!(result.is_some(), "Faint star should produce a result");
        let result = result.unwrap();
        assert!(result.converged, "Faint star should converge");
        assert!((result.sigma_x - sigma).abs() < 0.5);
    }

    #[test]
    fn test_evaluate_moffat_2d_center() {
        let params = [100.0, 5000.0, 0.0, 0.0, 3.0, 3.0, 0.0, 3.0];
        let val = evaluate_moffat_2d(&params, 0.0, 0.0);
        assert!((val - 5100.0).abs() < 0.01, "At center: B + A = 5100, got {}", val);
    }

    #[test]
    fn test_evaluate_moffat_2d_off_center() {
        let params = [100.0, 5000.0, 0.0, 0.0, 3.0, 3.0, 0.0, 3.0];
        let val_center = evaluate_moffat_2d(&params, 0.0, 0.0);
        let val_off = evaluate_moffat_2d(&params, 3.0, 0.0);
        assert!(val_off < val_center, "Off-center should be less");
        assert!(val_off > 100.0, "Should be above background");
    }

    #[test]
    fn test_evaluate_moffat_2d_far_away() {
        let params = [100.0, 5000.0, 0.0, 0.0, 3.0, 3.0, 0.0, 3.0];
        let val = evaluate_moffat_2d(&params, 100.0, 100.0);
        assert!((val - 100.0).abs() < 1.0, "Far from center should approach background, got {}", val);
    }
}
