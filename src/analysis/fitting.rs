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

/// Result from Levenberg-Marquardt solver: convergence flag + final cost.
///
/// `#[non_exhaustive]`: this crate is published, and the solvers grow a
/// reporting field now and then (`iters` arrived with the C3 arithmetic
/// pins). Downstream code reads the fields and never builds or exhaustively
/// destructures the struct, so sealing construction costs nothing and keeps
/// the next field from being a breaking change.
#[allow(dead_code)]
#[non_exhaustive]
pub struct LmResult {
    pub converged: bool,
    pub final_cost: f64,
    /// Iterations of the LM loop actually entered. Not used by the fitters
    /// themselves; the C3 arithmetic pins assert it is unchanged, which is
    /// what "the control flow is untouched" means in practice.
    pub iters: usize,
}

/// Largest normal-equation system either solver builds: the 8 free Moffat
/// parameters. Both solvers keep their Cholesky scratch on the stack at this
/// size — the allocating shape this replaced took three `Vec`s per LM
/// iteration, on the per-star hot path of Measure, LN and Register.
const MAX_LM_PARAMS: usize = 8;

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
    /// Normalized fit residual: sqrt(cost / n_pixels) / amplitude.
    /// Lower = better fit. 0 = perfect, ~1 = moments fallback.
    pub fit_residual: f64,
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
    let LmResult {
        converged,
        final_cost,
        ..
    } = lm_solve_2d(pixels, &mut params, true, max_iter, conv_tol, max_rejects);

    let sigma_x = params[4].abs();
    let sigma_y = params[5].abs();
    if sigma_x < 0.3 || sigma_y < 0.3 || params[1] <= 0.0 {
        return None;
    }

    // Normalized residual: RMS per pixel / amplitude (dimensionless)
    let fit_residual = (final_cost / pixels.len() as f64).sqrt() / params[1].max(1e-10);

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
        fit_residual,
    })
}

/// LM solver for 2D Gaussian. If `fit_theta` is true, uses 7 params; else 6 (theta=0).
fn lm_solve_2d(pixels: &[PixelSample], params: &mut [f64], fit_theta: bool, max_iter: usize, conv_tol: f64, max_rejects: usize) -> LmResult {
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
    let mut chol_l = [0.0_f64; MAX_LM_PARAMS * MAX_LM_PARAMS];
    let mut chol_y = [0.0_f64; MAX_LM_PARAMS];
    let mut delta_buf = [0.0_f64; MAX_LM_PARAMS];
    let mut iters = 0_usize;

    for _ in 0..max_iter {
        iters += 1;
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

        let solved = {
            let mut ok =
                cholesky_solve_into(&mat, &jtr, np, &mut chol_l, &mut chol_y, &mut delta_buf);
            if !ok {
                for _ in 0..3 {
                    lambda *= 10.0;
                    for p in 0..np {
                        mat[p * np + p] = jtj[p * np + p] + lambda * jtj[p * np + p].max(1e-12);
                    }
                    ok = cholesky_solve_into(
                        &mat,
                        &jtr,
                        np,
                        &mut chol_l,
                        &mut chol_y,
                        &mut delta_buf,
                    );
                    if ok {
                        break;
                    }
                }
            }
            ok
        };
        if !solved {
            break;
        }
        let delta = &delta_buf[..np];

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

    LmResult {
        converged,
        final_cost: best_cost,
        iters,
    }
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
///
/// Writes the solution into `out[..np]` and returns `true`; returns `false`
/// with `out` untouched when the matrix is not positive definite. `l` and `y`
/// are caller-owned scratch of at least `np*np` and `np` elements — the
/// arithmetic is term-for-term the allocating version this replaced, only the
/// three `Vec`s moved to the caller's stack.
fn cholesky_solve_into(
    mat: &[f64],
    rhs: &[f64],
    np: usize,
    l: &mut [f64],
    y: &mut [f64],
    out: &mut [f64],
) -> bool {
    // Cholesky: A = L * L^T
    let l = &mut l[..np * np];
    l.fill(0.0);

    for i in 0..np {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l[i * np + k] * l[j * np + k];
            }
            if i == j {
                let diag = mat[i * np + i] - sum;
                if diag <= 0.0 {
                    return false; // Not positive definite
                }
                l[i * np + j] = diag.sqrt();
            } else {
                l[i * np + j] = (mat[i * np + j] - sum) / l[j * np + j];
            }
        }
    }

    // Solve L * y = rhs (forward substitution)
    for i in 0..np {
        let mut sum = 0.0;
        for j in 0..i {
            sum += l[i * np + j] * y[j];
        }
        y[i] = (rhs[i] - sum) / l[i * np + i];
    }

    // Solve L^T * x = y (back substitution)
    for i in (0..np).rev() {
        let mut sum = 0.0;
        for j in (i + 1)..np {
            sum += l[j * np + i] * out[j]; // L^T[i][j] = L[j][i]
        }
        out[i] = (y[i] - sum) / l[i * np + i];
    }

    true
}

/// Allocating wrapper, kept so the decomposition keeps its own unit pin and so
/// the frozen pre-C3 reference solver can be written exactly as it was.
#[cfg(test)]
fn cholesky_solve(mat: &[f64], rhs: &[f64], np: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0_f64; np * np];
    let mut y = vec![0.0_f64; np];
    let mut out = vec![0.0_f64; np];
    cholesky_solve_into(mat, rhs, np, &mut l, &mut y, &mut out).then_some(out)
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
    /// Normalized fit residual: sqrt(cost / n_pixels) / amplitude.
    pub fit_residual: f64,
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
    fit_moffat_2d_impl(pixels, init_b, init_a, init_x0, init_y0, init_sigma_x, init_sigma_y, init_theta, BetaParam::FreeLinear, max_iter, conv_tol, max_rejects, ThetaColumn::Legacy)
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
    fit_moffat_2d_impl(pixels, init_b, init_a, init_x0, init_y0, init_sigma_x, init_sigma_y, init_theta, BetaParam::Fixed(beta), max_iter, conv_tol, max_rejects, ThetaColumn::Legacy)
}

/// Fit 2D elliptical Moffat with fixed beta (7 free parameters), the θ column of the Jacobian carrying its correct
/// sign, `∂Q/∂θ = 2uv·(1/α_x² − 1/α_y²)`.
///
/// [`fit_moffat_2d`] and [`fit_moffat_2d_fixed_beta`] build that column with the opposite sign (every release
/// through 1.3.0). With it the solver's θ step points uphill: once the other parameters have absorbed what they can,
/// every step is rejected and the fit stops away from the minimum — on a noiseless elliptical star started 5° off,
/// at 44.8° for a true 30°. Those two entry points keep the column they shipped, so what they return does not
/// change; this one is the corrected fit, otherwise identical (same arguments, same solver, same acceptance rules).
///
/// Named for the conditioned stamps its first caller hands it (centred on the stamp, amplitude-normalised, so every
/// parameter is O(1)); it is correct on any input.
pub fn fit_moffat_2d_fixed_beta_conditioned(
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
    fit_moffat_2d_impl(pixels, init_b, init_a, init_x0, init_y0, init_sigma_x, init_sigma_y, init_theta, BetaParam::Fixed(beta), max_iter, conv_tol, max_rejects, ThetaColumn::Corrected)
}

/// Result of [`fit_moffat_2d_bounded`]: the parameters of [`Moffat2DResult`] plus the solver's iteration count.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Moffat2DBoundedResult {
    pub b: f64,
    pub a: f64,
    pub x0: f64,
    pub y0: f64,
    pub alpha_x: f64,
    pub alpha_y: f64,
    pub theta: f64,
    /// Inside the caller's bounds; exactly a bound when the fit ended there.
    pub beta: f64,
    pub converged: bool,
    /// Normalized fit residual: sqrt(cost / n_pixels) / amplitude.
    pub fit_residual: f64,
    /// LM iterations the solve ran.
    pub iters: usize,
}

impl Moffat2DBoundedResult {
    /// FWHM = 2α√(2^(1/β) - 1) along each axis.
    pub fn fwhm_x(&self) -> f64 {
        2.0 * self.alpha_x * (2.0_f64.powf(1.0 / self.beta) - 1.0).sqrt()
    }
    pub fn fwhm_y(&self) -> f64 {
        2.0 * self.alpha_y * (2.0_f64.powf(1.0 / self.beta) - 1.0).sqrt()
    }
}

/// Fit 2D elliptical Moffat with 8 parameters, β free inside the caller's `beta_bounds` `(lo, hi)`, `0 < lo ≤ hi`.
///
/// The solver carries ln β rather than β: as β grows the profile approaches a Gaussian and β trades against α along
/// a flat ridge, which a step in ln β crosses in a few iterations where a step in β stalls. Each step's ln β is
/// clamped to `[ln lo, ln hi]`. A fit that ends at a bound is returned with `beta` exactly that bound — the caller
/// decides what a bound means; only a non-finite parameter or residual, `α < 0.3` or `A ≤ 0` refuses it (`None`), as
/// does fewer than 12 pixels or bounds that are not `0 < lo ≤ hi < ∞`. `init_beta` is clamped into the bounds and
/// sets the start α from the Gaussian σ as [`fit_moffat_2d`] does from its fixed 3.0. The θ column carries its
/// correct sign, as in [`fit_moffat_2d_fixed_beta_conditioned`].
///
/// [`fit_moffat_2d`] (β linear, bounds [1, 20], the legacy θ column) is unchanged: this is a separate entry point
/// on the same solver.
pub fn fit_moffat_2d_bounded(
    pixels: &[PixelSample],
    init_b: f64,
    init_a: f64,
    init_x0: f64,
    init_y0: f64,
    init_sigma_x: f64,
    init_sigma_y: f64,
    init_theta: f64,
    init_beta: f64,
    beta_bounds: (f64, f64),
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
) -> Option<Moffat2DBoundedResult> {
    let (lo, hi) = beta_bounds;
    if !(lo > 0.0 && lo <= hi && hi.is_finite() && init_beta.is_finite()) {
        return None;
    }
    let (ln_min, ln_max) = (lo.ln(), hi.ln());
    let MoffatSolve { params, lm } = solve_moffat_2d(
        pixels,
        init_b,
        init_a,
        init_x0,
        init_y0,
        init_sigma_x,
        init_sigma_y,
        init_theta,
        BetaParam::FreeLn { ln_min, ln_max },
        init_beta.clamp(lo, hi),
        max_iter,
        conv_tol,
        max_rejects,
        ThetaColumn::Corrected,
    )?;
    if !(params.iter().all(|p| p.is_finite()) && lm.final_cost.is_finite()) {
        return None;
    }
    let alpha_x = params[4].abs();
    let alpha_y = params[5].abs();
    if alpha_x < 0.3 || alpha_y < 0.3 || params[1] <= 0.0 {
        return None;
    }
    // At a bound the step clamp left ln β exactly there; report the bound itself, not `exp(ln hi)`'s last bits.
    let beta = if params[7] >= ln_max {
        hi
    } else if params[7] <= ln_min {
        lo
    } else {
        params[7].exp().clamp(lo, hi)
    };
    Some(Moffat2DBoundedResult {
        b: params[0],
        a: params[1],
        x0: params[2],
        y0: params[3],
        alpha_x,
        alpha_y,
        theta: normalize_theta(params[6]),
        beta,
        converged: lm.converged,
        fit_residual: (lm.final_cost / pixels.len() as f64).sqrt() / params[1].max(1e-10),
        iters: lm.iters,
    })
}

/// How the Moffat solver carries β: what `params[7]` holds, how many parameters are free, the β column of the
/// Jacobian and the step's clamp. The one private difference between [`fit_moffat_2d_bounded`] and the entry points
/// that shipped before it, beside [`ThetaColumn`].
#[derive(Clone, Copy, Debug, PartialEq)]
enum BetaParam {
    /// β held at this value, carried in `params[7]` (7 free parameters): [`fit_moffat_2d_fixed_beta`] and
    /// [`fit_moffat_2d_fixed_beta_conditioned`].
    Fixed(f64),
    /// β free, `params[7] = β`: the column `∂M/∂β = A·(1+Q)^(−β)·(−ln(1+Q))`, each step's β clamped to [1.5, 20] and
    /// the result refused outside [1, 20]. [`fit_moffat_2d`]'s arm, operation for operation what it always ran, so its
    /// output is bit-identical.
    FreeLinear,
    /// β free, `params[7] = ln β`: the column `∂M/∂(ln β) = A·(1+Q)^(−β)·(−ln(1+Q))·β` (the chain rule), each step's
    /// ln β clamped to `[ln_min, ln_max]`. [`fit_moffat_2d_bounded`]'s arm.
    FreeLn { ln_min: f64, ln_max: f64 },
}

impl BetaParam {
    /// The legacy solver's `Option`: `Some(b)` fixed, `None` the linear free arm.
    #[cfg(test)]
    fn fixed_or_linear(fixed_beta: Option<f64>) -> Self {
        fixed_beta.map_or(BetaParam::FreeLinear, BetaParam::Fixed)
    }

    fn is_free(self) -> bool {
        !matches!(self, BetaParam::Fixed(_))
    }

    /// β at the parameter vector's `params[7]`: the value itself, except in ln β.
    #[inline(always)]
    fn beta(self, p7: f64) -> f64 {
        match self {
            BetaParam::Fixed(_) | BetaParam::FreeLinear => p7,
            BetaParam::FreeLn { .. } => p7.exp(),
        }
    }

    /// The β column of one Jacobian row from `a_power = A·(1+Q)^(−β)`, `base = 1 + Q` and β. `FreeLinear` is the
    /// expression the solver always evaluated, operand for operand (`A·power` first, then `·(−ln base)`).
    #[inline(always)]
    fn column(self, a_power: f64, base: f64, beta: f64) -> f64 {
        match self {
            BetaParam::Fixed(_) => 0.0,
            BetaParam::FreeLinear => a_power * (-base.ln()),
            BetaParam::FreeLn { .. } => a_power * (-base.ln()) * beta,
        }
    }

    /// Keep a step's β (`new_p7`) in range: [1.5, 20] in β for `FreeLinear` (as always), `[ln_min, ln_max]` in ln β for
    /// `FreeLn`; a fixed β is never stepped.
    #[inline(always)]
    fn clamp_step(self, new_p7: &mut f64) {
        match self {
            BetaParam::Fixed(_) => {}
            BetaParam::FreeLinear => {
                if *new_p7 < 1.5 { *new_p7 = 1.5; }
                if *new_p7 > 20.0 { *new_p7 = 20.0; }
            }
            BetaParam::FreeLn { ln_min, ln_max } => {
                if *new_p7 < ln_min { *new_p7 = ln_min; }
                if *new_p7 > ln_max { *new_p7 = ln_max; }
            }
        }
    }
}

/// θ normalised to [-π/2, π/2].
fn normalize_theta(t: f64) -> f64 {
    let mut theta = t % std::f64::consts::PI;
    if theta > std::f64::consts::FRAC_PI_2 {
        theta -= std::f64::consts::PI;
    } else if theta < -std::f64::consts::FRAC_PI_2 {
        theta += std::f64::consts::PI;
    }
    theta
}

/// Which sign the Moffat solver's θ column carries — the one private difference between
/// [`fit_moffat_2d_fixed_beta_conditioned`] and the two entry points that shipped before it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ThetaColumn {
    /// `2uv·(1/α_y² − 1/α_x²)`, the NEGATIVE of `∂Q/∂θ`: what [`fit_moffat_2d`] and [`fit_moffat_2d_fixed_beta`]
    /// have always built, kept so their output does not move.
    Legacy,
    /// `∂Q/∂θ = 2uv·(1/α_x² − 1/α_y²)`, from `∂u/∂θ = v`, `∂v/∂θ = −u`.
    Corrected,
}

impl ThetaColumn {
    /// The θ factor of the Jacobian row, `∂Q/∂θ` up to this column's sign. `Legacy` is the expression the solver
    /// always evaluated, operand for operand, so its result is bit-identical.
    #[inline(always)]
    fn dq_dtheta(self, u: f64, v: f64, inv_ax2: f64, inv_ay2: f64) -> f64 {
        match self {
            ThetaColumn::Legacy => 2.0 * u * v * (inv_ay2 - inv_ax2),
            ThetaColumn::Corrected => 2.0 * u * v * (inv_ax2 - inv_ay2),
        }
    }
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

/// One Moffat solve before acceptance: the parameters the LM run ended at (`params[7]` in `beta_param`'s own
/// parametrisation) and the run's outcome.
struct MoffatSolve {
    params: [f64; 8],
    lm: LmResult,
}

/// The start vector every Moffat entry point builds — α from the Gaussian σ at `beta_init`,
/// α = σ × 2.3548 / (2√(2^(1/β₀) − 1)) — then the LM solve. `None` under 12 pixels.
#[allow(clippy::too_many_arguments)]
fn solve_moffat_2d(
    pixels: &[PixelSample],
    init_b: f64,
    init_a: f64,
    init_x0: f64,
    init_y0: f64,
    init_sigma_x: f64,
    init_sigma_y: f64,
    init_theta: f64,
    beta_param: BetaParam,
    beta_init: f64,
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
    theta_column: ThetaColumn,
) -> Option<MoffatSolve> {
    if pixels.len() < 12 {
        return None;
    }

    // Convert Gaussian sigma to Moffat alpha
    let moffat_scale = (2.0_f64.powf(1.0 / beta_init) - 1.0).sqrt();
    let alpha_x = init_sigma_x * 2.3548 / (2.0 * moffat_scale);
    let alpha_y = init_sigma_y * 2.3548 / (2.0 * moffat_scale);
    let p7 = match beta_param {
        BetaParam::FreeLn { .. } => beta_init.ln(),
        BetaParam::Fixed(_) | BetaParam::FreeLinear => beta_init,
    };

    let mut params = [init_b, init_a, init_x0, init_y0, alpha_x, alpha_y, init_theta, p7];
    let lm = lm_solve_moffat_impl(
        pixels,
        &mut params,
        beta_param,
        max_iter,
        conv_tol,
        max_rejects,
        theta_column,
    );
    Some(MoffatSolve { params, lm })
}

/// [`fit_moffat_2d`], [`fit_moffat_2d_fixed_beta`] and [`fit_moffat_2d_fixed_beta_conditioned`]: `beta_param` is
/// `Fixed` or `FreeLinear` (start β 3.0), and the acceptance is the one they always applied — `α < 0.3`, `A ≤ 0`
/// or β outside [1, 20] refuses the fit.
#[allow(clippy::too_many_arguments)]
fn fit_moffat_2d_impl(
    pixels: &[PixelSample],
    init_b: f64,
    init_a: f64,
    init_x0: f64,
    init_y0: f64,
    init_sigma_x: f64,
    init_sigma_y: f64,
    init_theta: f64,
    beta_param: BetaParam,
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
    theta_column: ThetaColumn,
) -> Option<Moffat2DResult> {
    debug_assert!(
        !matches!(beta_param, BetaParam::FreeLn { .. }),
        "the ln-beta arm is fit_moffat_2d_bounded's"
    );
    let beta_init = match beta_param {
        BetaParam::Fixed(b) => b,
        _ => 3.0,
    };
    let MoffatSolve {
        params,
        lm: LmResult {
            converged,
            final_cost,
            ..
        },
    } = solve_moffat_2d(
        pixels,
        init_b,
        init_a,
        init_x0,
        init_y0,
        init_sigma_x,
        init_sigma_y,
        init_theta,
        beta_param,
        beta_init,
        max_iter,
        conv_tol,
        max_rejects,
        theta_column,
    )?;

    let alpha_x = params[4].abs();
    let alpha_y = params[5].abs();
    let beta = if beta_param.is_free() { params[7] } else { beta_init };

    // Reject non-physical results
    if alpha_x < 0.3 || alpha_y < 0.3 || params[1] <= 0.0
        || beta < 1.0 || beta > 20.0
    {
        return None;
    }

    // Normalized residual: RMS per pixel / amplitude (dimensionless)
    let fit_residual = (final_cost / pixels.len() as f64).sqrt() / params[1].max(1e-10);

    Some(Moffat2DResult {
        b: params[0],
        a: params[1],
        x0: params[2],
        y0: params[3],
        alpha_x,
        alpha_y,
        theta: normalize_theta(params[6]),
        beta,
        converged,
        fit_residual,
    })
}

/// `(1 + Q)^(-β)`, with the integral-β case on `powi`.
///
/// The stacking pipeline fits at a fixed β that is always one of
/// `AUTO_BETAS`-style candidates, three of which (4, 6, 10) are integers, and
/// `powi` is a short chain of multiplies and one reciprocal where `powf` is a
/// full `log`/`exp` pair. Resolved once per Jacobian/cost pass, never per
/// sample. The two arms differ in the last bits; the C3 pins record by how
/// much on real star stamps.
#[derive(Clone, Copy)]
enum MoffatExponent {
    Integral(i32),
    Real(f64),
}

impl MoffatExponent {
    fn new(beta: f64) -> Self {
        let n = beta as i32;
        if (1..=20).contains(&n) && n as f64 == beta {
            Self::Integral(-n)
        } else {
            Self::Real(-beta)
        }
    }

    #[inline(always)]
    fn eval(self, base: f64) -> f64 {
        match self {
            Self::Integral(n) => base.powi(n),
            Self::Real(e) => base.powf(e),
        }
    }
}

/// Per-sample state of the parameter vector the Moffat solver currently holds
/// — `(1+Q)^(-β)` and the residual — carried from the cost evaluation that
/// produced it into the next Jacobian pass instead of being recomputed there
/// (C3). That removes the Jacobian's own `powf` per sample per iteration; the
/// second one it used to spend on `(1+Q)^(-β-1)` is gone too, that factor
/// being `power / base` exactly.
///
/// Thread-local so the per-star hot path allocates nothing: one instance per
/// rayon worker, never shared, never `Send`.
#[derive(Default)]
struct MoffatScratch {
    cur_power: Vec<f64>,
    cur_resid: Vec<f64>,
    cand_power: Vec<f64>,
    cand_resid: Vec<f64>,
}

impl MoffatScratch {
    /// Grow every buffer to `n` without clearing: every element is written
    /// before it is read, so retained capacity is invisible to the numbers.
    fn reserve_exact_len(&mut self, n: usize) {
        for v in [
            &mut self.cur_power,
            &mut self.cur_resid,
            &mut self.cand_power,
            &mut self.cand_resid,
        ] {
            if v.len() < n {
                v.resize(n, 0.0);
            }
        }
    }
}

thread_local! {
    static MOFFAT_SCRATCH: std::cell::RefCell<MoffatScratch> =
        std::cell::RefCell::new(MoffatScratch::default());
}

/// LM solver for 2D Moffat. `BetaParam::Fixed(b)` holds β at `b` and optimizes the other 7 parameters; the free
/// arms optimize all 8, with β (`FreeLinear`) or ln β (`FreeLn`) in `params[7]`.
fn lm_solve_moffat_impl(
    pixels: &[PixelSample],
    params: &mut [f64],
    beta_param: BetaParam,
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
    theta_column: ThetaColumn,
) -> LmResult {
    // The solver never re-enters itself, so the borrow always succeeds; the
    // fallback keeps a future caller from panicking rather than degrading.
    MOFFAT_SCRATCH.with(|cell| match cell.try_borrow_mut() {
        Ok(mut scratch) => lm_solve_moffat_with(
            pixels,
            params,
            beta_param,
            max_iter,
            conv_tol,
            max_rejects,
            theta_column,
            &mut scratch,
        ),
        Err(_) => lm_solve_moffat_with(
            pixels,
            params,
            beta_param,
            max_iter,
            conv_tol,
            max_rejects,
            theta_column,
            &mut MoffatScratch::default(),
        ),
    })
}

#[allow(clippy::too_many_arguments)]
fn lm_solve_moffat_with(
    pixels: &[PixelSample],
    params: &mut [f64],
    beta_param: BetaParam,
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
    theta_column: ThetaColumn,
    scratch: &mut MoffatScratch,
) -> LmResult {
    // The reuse seam's contract: `residual_cost_moffat_into` resolves its
    // exponent from `beta_param.beta(params[7])`, while the Jacobian pass
    // below takes the fixed β itself for `Fixed(b)` — a cached `power` is the
    // value the Jacobian wants only while those two agree. They do by
    // construction (`solve_moffat_2d` seeds `params[7] = beta_init = b`, and
    // `np = 7` keeps `delta` from ever reaching index 7), but that is a
    // cross-function rule, so state it rather than leave it unwritten.
    debug_assert!(
        match beta_param {
            BetaParam::Fixed(b) => params[7] == b,
            BetaParam::FreeLinear | BetaParam::FreeLn { .. } => true,
        },
        "a fixed-beta solve must carry its beta in params[7]: {} != {beta_param:?}",
        params[7]
    );
    let np = if beta_param.is_free() { 8 } else { 7 };
    let n = pixels.len();
    scratch.reserve_exact_len(n);
    let MoffatScratch {
        cur_power,
        cur_resid,
        cand_power,
        cand_resid,
    } = scratch;
    let (mut cur_power, mut cur_resid) = (&mut cur_power[..n], &mut cur_resid[..n]);
    let (mut cand_power, mut cand_resid) = (&mut cand_power[..n], &mut cand_resid[..n]);

    let mut lambda = 1e-3_f64;
    let mut nu = 2.0_f64;
    // Seeds the cache: `cur_*` describe `params` from here on.
    let mut best_cost = residual_cost_moffat_into(pixels, params, beta_param, cur_power, cur_resid);
    let mut prev_cost: f64;
    let mut converged = false;
    let mut consecutive_rejects: usize = 0;

    let mut jtj = vec![0.0_f64; np * np];
    let mut jtr = vec![0.0_f64; np];
    let mut j = vec![0.0_f64; np];
    let mut mat = vec![0.0_f64; np * np];
    let mut new_params = [0.0_f64; 8];
    let mut chol_l = [0.0_f64; MAX_LM_PARAMS * MAX_LM_PARAMS];
    let mut chol_y = [0.0_f64; MAX_LM_PARAMS];
    let mut delta_buf = [0.0_f64; MAX_LM_PARAMS];
    let mut iters = 0_usize;

    for _ in 0..max_iter {
        iters += 1;
        jtj.fill(0.0);
        jtr.fill(0.0);

        let theta = params[6];
        let (cos_t, sin_t) = (theta.cos(), theta.sin());
        let ax = params[4];
        let ay = params[5];
        let beta = match beta_param {
            BetaParam::Fixed(b) => b,
            BetaParam::FreeLinear | BetaParam::FreeLn { .. } => beta_param.beta(params[7]),
        };
        let inv_ax2 = 1.0 / (ax * ax);
        let inv_ay2 = 1.0 / (ay * ay);

        for (i, px) in pixels.iter().enumerate() {
            let dx = px.x - params[2];
            let dy = px.y - params[3];
            let u = dx * cos_t + dy * sin_t;
            let v = -dx * sin_t + dy * cos_t;
            let q = u * u * inv_ax2 + v * v * inv_ay2;
            let base = 1.0 + q;
            // Both come from the cost evaluation that accepted this parameter
            // vector (or, on the first pass, from the seeding call above) —
            // the same expressions, at the same `params`, so the value is the
            // one this loop used to recompute.
            let power = cur_power[i];
            let r = cur_resid[i];

            // Jacobian of M(x,y) = B + A × (1+Q)^(-β)
            j[0] = 1.0; // dM/dB
            j[1] = power; // dM/dA

            // (1+Q)^(-β-1) is (1+Q)^(-β) / (1+Q): one division in place of a
            // second `powf`.
            let dpower_dq = -beta * (power / base);
            let a_dpq = params[1] * dpower_dq;

            let dq_dx0 = -2.0 * (cos_t * u * inv_ax2 - sin_t * v * inv_ay2);
            j[2] = a_dpq * dq_dx0;

            let dq_dy0 = -2.0 * (sin_t * u * inv_ax2 + cos_t * v * inv_ay2);
            j[3] = a_dpq * dq_dy0;

            j[4] = a_dpq * (-2.0 * u * u / (ax * ax * ax));
            j[5] = a_dpq * (-2.0 * v * v / (ay * ay * ay));

            let dq_dtheta = theta_column.dq_dtheta(u, v, inv_ax2, inv_ay2);
            j[6] = a_dpq * dq_dtheta;

            // Only include the β column when β is free — the one `ln` in the
            // solver, and the reason the normal system is 8×8 instead of 7×7.
            if beta_param.is_free() {
                j[7] = beta_param.column(params[1] * power, base, beta);
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

        let solved = {
            let mut ok = cholesky_solve_into(
                &mat[..np * np],
                &jtr,
                np,
                &mut chol_l,
                &mut chol_y,
                &mut delta_buf,
            );
            if !ok {
                for _ in 0..3 {
                    lambda *= 10.0;
                    for p in 0..np {
                        mat[p * np + p] = jtj[p * np + p] + lambda * jtj[p * np + p].max(1e-12);
                    }
                    ok = cholesky_solve_into(
                        &mat[..np * np],
                        &jtr,
                        np,
                        &mut chol_l,
                        &mut chol_y,
                        &mut delta_buf,
                    );
                    if ok {
                        break;
                    }
                }
            }
            ok
        };
        if !solved {
            break;
        }
        let delta = &delta_buf[..np];

        new_params[..8].copy_from_slice(&params[..8]);
        for p in 0..np {
            new_params[p] += delta[p];
        }
        // Keep alphas positive
        if new_params[4] <= 0.0 { new_params[4] = params[4] * 0.5; }
        if new_params[5] <= 0.0 { new_params[5] = params[5] * 0.5; }
        // Keep beta in range (only when free)
        beta_param.clamp_step(&mut new_params[7]);

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
            // `params` did not move, so `cur_*` still describe it.
            continue;
        }

        let new_cost =
            residual_cost_moffat_into(pixels, &new_params, beta_param, cand_power, cand_resid);

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
                // The candidate became `params`, so its per-sample state
                // becomes the cache the next Jacobian pass reads.
                std::mem::swap(&mut cur_power, &mut cand_power);
                std::mem::swap(&mut cur_resid, &mut cand_resid);
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

    LmResult {
        converged,
        final_cost: best_cost,
        iters,
    }
}

/// Sum of squared residuals at `params`, recording each sample's
/// `(1+Q)^(-β)` and residual so the next Jacobian pass can read them instead
/// of evaluating the power again. β is `beta_param.beta(params[7])`.
fn residual_cost_moffat_into(
    pixels: &[PixelSample],
    params: &[f64],
    beta_param: BetaParam,
    power_out: &mut [f64],
    resid_out: &mut [f64],
) -> f64 {
    let theta = params[6];
    let (cos_t, sin_t) = (theta.cos(), theta.sin());
    let inv_ax2 = 1.0 / (params[4] * params[4]);
    let inv_ay2 = 1.0 / (params[5] * params[5]);
    let exponent = MoffatExponent::new(beta_param.beta(params[7]));

    let mut cost = 0.0_f64;
    for (i, px) in pixels.iter().enumerate() {
        let dx = px.x - params[2];
        let dy = px.y - params[3];
        let u = dx * cos_t + dy * sin_t;
        let v = -dx * sin_t + dy * cos_t;
        let q = u * u * inv_ax2 + v * v * inv_ay2;
        let power = exponent.eval(1.0 + q);
        let model = params[0] + params[1] * power;
        let r = px.value - model;
        power_out[i] = power;
        resid_out[i] = r;
        cost += r * r;
    }
    cost
}

/// The Moffat LM solver exactly as it stood before the C3 arithmetic change
/// (Tier C, spec `2026-09-20-stacking-compute-tierC-design.md` §4): two
/// `powf` per sample per iteration, the residual recomputed by every
/// Jacobian pass, an allocating Cholesky. Frozen on purpose — it is the
/// oracle the DELTA pins below measure the shipped solver against, so a later
/// regression is caught at the fixture level. Do not "improve" it.
#[cfg(test)]
mod reference {
    use super::{LmResult, PixelSample, ThetaColumn};

    pub fn cholesky_solve_alloc(mat: &[f64], rhs: &[f64], np: usize) -> Option<Vec<f64>> {
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
                        return None;
                    }
                    l[i * np + j] = diag.sqrt();
                } else {
                    l[i * np + j] = (mat[i * np + j] - sum) / l[j * np + j];
                }
            }
        }
        let mut y = vec![0.0_f64; np];
        for i in 0..np {
            let mut sum = 0.0;
            for j in 0..i {
                sum += l[i * np + j] * y[j];
            }
            y[i] = (rhs[i] - sum) / l[i * np + i];
        }
        let mut x = vec![0.0_f64; np];
        for i in (0..np).rev() {
            let mut sum = 0.0;
            for j in (i + 1)..np {
                sum += l[j * np + i] * x[j];
            }
            x[i] = (y[i] - sum) / l[i * np + i];
        }
        Some(x)
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

    pub fn lm_solve_moffat(
        pixels: &[PixelSample],
        params: &mut [f64],
        fixed_beta: Option<f64>,
        max_iter: usize,
        conv_tol: f64,
        max_rejects: usize,
        theta_column: ThetaColumn,
    ) -> LmResult {
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
        let mut iters = 0_usize;

        for _ in 0..max_iter {
            iters += 1;
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

                j[0] = 1.0;
                j[1] = power;

                let dpower_dq = -beta * base.powf(-beta - 1.0);
                let a_dpq = params[1] * dpower_dq;

                let dq_dx0 = -2.0 * (cos_t * u * inv_ax2 - sin_t * v * inv_ay2);
                j[2] = a_dpq * dq_dx0;

                let dq_dy0 = -2.0 * (sin_t * u * inv_ax2 + cos_t * v * inv_ay2);
                j[3] = a_dpq * dq_dy0;

                j[4] = a_dpq * (-2.0 * u * u / (ax * ax * ax));
                j[5] = a_dpq * (-2.0 * v * v / (ay * ay * ay));

                let dq_dtheta = theta_column.dq_dtheta(u, v, inv_ax2, inv_ay2);
                j[6] = a_dpq * dq_dtheta;

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

            for p in 0..np {
                for q in 0..p {
                    jtj[p * np + q] = jtj[q * np + p];
                }
            }

            mat[..np * np].copy_from_slice(&jtj);
            for p in 0..np {
                mat[p * np + p] += lambda * jtj[p * np + p].max(1e-12);
            }

            let delta = {
                let mut result = cholesky_solve_alloc(&mat[..np * np], &jtr, np);
                if result.is_none() {
                    for _ in 0..3 {
                        lambda *= 10.0;
                        for p in 0..np {
                            mat[p * np + p] = jtj[p * np + p] + lambda * jtj[p * np + p].max(1e-12);
                        }
                        result = cholesky_solve_alloc(&mat[..np * np], &jtr, np);
                        if result.is_some() {
                            break;
                        }
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
            if new_params[4] <= 0.0 {
                new_params[4] = params[4] * 0.5;
            }
            if new_params[5] <= 0.0 {
                new_params[5] = params[5] * 0.5;
            }
            if fixed_beta.is_none() {
                if new_params[7] < 1.5 {
                    new_params[7] = 1.5;
                }
                if new_params[7] > 20.0 {
                    new_params[7] = 20.0;
                }
            }

            if new_params[4] < 0.5 || new_params[5] < 0.5 || {
                let cdx = new_params[2] - params[2];
                let cdy = new_params[3] - params[3];
                cdx * cdx + cdy * cdy > 4.0
            } {
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

        LmResult {
            converged,
            final_cost: best_cost,
            iters,
        }
    }
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

    #[test]
    fn test_gaussian_fit_residual_near_zero() {
        // Perfect synthetic Gaussian → fit_residual should be ≈ 0
        let size = 21;
        let mut pixels = Vec::new();
        for y in 0..size {
            for x in 0..size {
                let dx = x as f64 - 10.0;
                let dy = y as f64 - 10.0;
                let v = 100.0 + 5000.0 * (-0.5 * (dx * dx + dy * dy) / 9.0).exp();
                pixels.push(PixelSample { x: x as f64, y: y as f64, value: v });
            }
        }
        let result = fit_gaussian_2d(&pixels, 100.0, 5000.0, 10.0, 10.0, 3.0, 3.0, 0.0,
            CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS).unwrap();
        assert!(result.fit_residual < 0.001,
            "Perfect Gaussian fit_residual should be ~0, got {}", result.fit_residual);
    }

    #[test]
    fn test_moffat_fit_residual_near_zero() {
        // Perfect synthetic Moffat → fit_residual should be ≈ 0
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
                pixels.push(PixelSample { x: x as f64, y: y as f64, value: v });
            }
        }
        let moffat_fwhm = 2.0 * alpha * (2.0_f64.powf(1.0 / beta) - 1.0).sqrt();
        let init_sigma = moffat_fwhm / 2.3548;
        let result = fit_moffat_2d(&pixels, 100.0, 5000.0, 15.0, 15.0,
            init_sigma, init_sigma, 0.0,
            CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS).unwrap();
        assert!(result.fit_residual < 0.001,
            "Perfect Moffat fit_residual should be ~0, got {}", result.fit_residual);
    }

    // ── C3 (Tier C) arithmetic pins ─────────────────────────────────────────
    //
    // The shipped Moffat solver evaluates ONE `powf` per sample per iteration
    // where the pre-C3 one evaluated three — two in the Jacobian pass, one in
    // the candidate's cost pass — takes `(1+Q)^(-β-1)` as `(1+Q)^(-β)/(1+Q)`,
    // and uses `powi` when β is integral. Each of those is the same
    // real-valued function rounded differently, so the LM fixed point moves in
    // its last bits. `reference::lm_solve_moffat` is the frozen pre-C3 code;
    // these pins measure the DELTA against it on realistic star stamps and
    // assert ruling C-4 — the control flow (iteration count, accept/reject,
    // convergence verdict) is untouched.

    /// Deterministic LCG: the fixtures must be the same bytes on every machine.
    struct Lcg(u64);

    impl Lcg {
        fn unit(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 11) as f64) / ((1_u64 << 53) as f64)
        }

        fn gauss(&mut self) -> f64 {
            let u1 = self.unit().max(1e-12);
            let u2 = self.unit();
            (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        }
    }

    /// One noisy elliptical star stamp plus the imperfect seed a detector
    /// would hand the fitter.
    struct Fixture {
        label: String,
        pixels: Vec<PixelSample>,
        seed_x: f64,
        seed_y: f64,
        seed_sx: f64,
        seed_sy: f64,
        seed_theta: f64,
        b0: f64,
        a0: f64,
    }

    /// The initial parameter vector `fit_moffat_2d_impl` builds, so the pins
    /// enter the solver in exactly the state production reaches it in.
    fn moffat_init(f: &Fixture, fixed_beta: Option<f64>) -> [f64; 8] {
        let beta_init = fixed_beta.unwrap_or(3.0);
        let moffat_scale = (2.0_f64.powf(1.0 / beta_init) - 1.0).sqrt();
        [
            f.b0,
            f.a0,
            f.seed_x,
            f.seed_y,
            f.seed_sx * 2.3548 / (2.0 * moffat_scale),
            f.seed_sy * 2.3548 / (2.0 * moffat_scale),
            f.seed_theta,
            beta_init,
        ]
    }

    /// 36 stamps over stamp radius, size, axis ratio, orientation, amplitude
    /// and true β, at two noise levels.
    fn fixture_stars() -> Vec<Fixture> {
        let mut rng = Lcg(0x5EED_1272_A11E_C3D5);
        let radii = [6_i64, 9, 14, 20];
        let alphas = [1.8_f64, 3.0, 5.0];
        let ratios = [1.0_f64, 1.4, 2.2];
        let thetas = [0.0_f64, 0.4, 1.1, -0.7];
        let amps = [60.0_f64, 900.0, 12_000.0];
        let true_betas = [2.2_f64, 3.5, 7.0];
        let bg = 250.0_f64;

        let mut out = Vec::new();
        for (i, &r) in radii.iter().enumerate() {
            for (k, &alpha) in alphas.iter().enumerate() {
                for (m, &amp) in amps.iter().enumerate() {
                    let ratio = ratios[(i + k) % ratios.len()];
                    let theta = thetas[(k + m) % thetas.len()];
                    let beta_true = true_betas[(i + m) % true_betas.len()];
                    let noise = if m == 0 { 1.5 } else { 12.0 };
                    let ax = alpha * ratio;
                    let ay = alpha / ratio;
                    let tx = 0.6 * (rng.unit() - 0.5);
                    let ty = 0.6 * (rng.unit() - 0.5);
                    let (ct, st) = (theta.cos(), theta.sin());

                    let mut pixels = Vec::with_capacity(((2 * r + 1) * (2 * r + 1)) as usize);
                    for y in -r..=r {
                        for x in -r..=r {
                            let dx = x as f64 - tx;
                            let dy = y as f64 - ty;
                            let u = dx * ct + dy * st;
                            let v = -dx * st + dy * ct;
                            let q = (u / ax).powi(2) + (v / ay).powi(2);
                            let value = bg + amp * (1.0 + q).powf(-beta_true) + noise * rng.gauss();
                            pixels.push(PixelSample {
                                x: x as f64,
                                y: y as f64,
                                value,
                            });
                        }
                    }

                    // A seed is never exact: centroid off by up to 0.2 px,
                    // size by ±20 %, orientation by ±0.15 rad.
                    let fwhm =
                        2.0 * (ax * ay).sqrt() * (2.0_f64.powf(1.0 / beta_true) - 1.0).sqrt();
                    let sigma = fwhm / 2.3548;
                    let scale = 1.0 + 0.4 * (rng.unit() - 0.5);
                    let qr = ratio.sqrt();
                    out.push(Fixture {
                        label: format!("r{r}/a{alpha}/amp{amp}/b{beta_true}/th{theta}"),
                        pixels,
                        seed_x: tx + 0.4 * (rng.unit() - 0.5),
                        seed_y: ty + 0.4 * (rng.unit() - 0.5),
                        seed_sx: sigma * scale * qr,
                        seed_sy: sigma * scale / qr,
                        seed_theta: theta + 0.3 * (rng.unit() - 0.5),
                        b0: bg,
                        a0: amp,
                    });
                }
            }
        }
        out
    }

    /// Total Moffat flux `A·π·α_x·α_y/(β-1)` — what a moved fixed point costs
    /// the metric the pipeline actually weights frames by.
    fn moffat_flux(p: &[f64; 8]) -> f64 {
        p[1] * std::f64::consts::PI * p[4].abs() * p[5].abs() / (p[7] - 1.0)
    }

    /// The C3 grid's stamp the golden free-β pin reads: [`fixture_stars`]' `r14/a3/amp900` stamp (noisy, elliptical,
    /// on a 250 background), its pixel coordinates shifted by +15 so the stamp centre sits at (15, 15).
    fn c3_fixture_stamp() -> Vec<PixelSample> {
        let f = fixture_stars()
            .into_iter()
            .find(|f| f.label.starts_with("r14/a3/amp900/"))
            .expect("the C3 grid holds an r14/a3/amp900 stamp");
        f.pixels
            .iter()
            .map(|p| PixelSample { x: p.x + 15.0, y: p.y + 15.0, value: p.value })
            .collect()
    }

    #[test]
    fn c3_moffat_arithmetic_delta_against_the_pre_c3_solver() {
        let fixtures = fixture_stars();
        assert!(fixtures.len() >= 30, "fixture count {}", fixtures.len());
        let modes: [Option<f64>; 5] = [Some(2.5), Some(4.0), Some(6.0), Some(10.0), None];

        let mut cases = 0_usize;
        let mut identical = 0_usize;
        let mut max_pos = 0.0_f64;
        let mut max_flux = 0.0_f64;
        let mut max_beta = 0.0_f64;
        let mut max_cost = 0.0_f64;
        let (mut worst_pos, mut worst_flux) = (String::new(), String::new());

        for f in &fixtures {
            for &fixed_beta in &modes {
                let init = moffat_init(f, fixed_beta);
                let tag = format!("{} fixed_beta={fixed_beta:?}", f.label);

                let mut now = init;
                let r_now = lm_solve_moffat_impl(
                    &f.pixels,
                    &mut now,
                    BetaParam::fixed_or_linear(fixed_beta),
                    CALIBRATION_MAX_ITER,
                    CALIBRATION_CONV_TOL,
                    CALIBRATION_MAX_REJECTS,
                    ThetaColumn::Legacy,
                );
                let mut was = init;
                let r_was = reference::lm_solve_moffat(
                    &f.pixels,
                    &mut was,
                    fixed_beta,
                    CALIBRATION_MAX_ITER,
                    CALIBRATION_CONV_TOL,
                    CALIBRATION_MAX_REJECTS,
                    ThetaColumn::Legacy,
                );

                // Ruling C-4 — arithmetic only, the control flow is untouched.
                assert_eq!(r_now.iters, r_was.iters, "iteration count moved: {tag}");
                assert_eq!(
                    r_now.converged, r_was.converged,
                    "convergence verdict moved: {tag}"
                );

                if let Some(b) = fixed_beta {
                    assert_eq!(now[7], b, "fixed beta left its pin: {tag}");
                    assert_eq!(now[7], was[7], "fixed beta moved: {tag}");
                }

                let pos = (now[2] - was[2]).abs().max((now[3] - was[3]).abs());
                assert!(pos <= 1e-5, "position moved {pos:.3e} px: {tag}");

                let (flux_now, flux_was) = (moffat_flux(&now), moffat_flux(&was));
                let flux = if flux_was != 0.0 {
                    ((flux_now - flux_was) / flux_was).abs()
                } else {
                    0.0
                };
                assert!(flux <= 1e-5, "flux moved {flux:.3e} relative: {tag}");

                let beta = if was[7] != 0.0 {
                    ((now[7] - was[7]) / was[7]).abs()
                } else {
                    0.0
                };
                assert!(beta <= 1e-5, "beta moved {beta:.3e} relative: {tag}");

                let cost = if r_was.final_cost > 0.0 {
                    ((r_now.final_cost - r_was.final_cost) / r_was.final_cost).abs()
                } else {
                    0.0
                };

                if now == was && r_now.final_cost == r_was.final_cost {
                    identical += 1;
                }
                if pos > max_pos {
                    max_pos = pos;
                    worst_pos = tag.clone();
                }
                if flux > max_flux {
                    max_flux = flux;
                    worst_flux = tag.clone();
                }
                max_beta = max_beta.max(beta);
                max_cost = max_cost.max(cost);
                cases += 1;
            }
        }

        // The loop ran: 36 stamps x 5 beta modes. Deliberately NOT an
        // assertion on how many of them moved — which fits land bit-identical
        // depends on the host's `powf`/`powi`.
        assert_eq!(cases, fixtures.len() * modes.len());
        assert_eq!(cases, 180, "fixture x mode grid");

        println!(
            "C3 DELTA over {cases} solves ({identical} bit-identical, \
             {} moved): max |dpos| {max_pos:.3e} px [{worst_pos}], \
             max |dflux|/flux {max_flux:.3e} [{worst_flux}], \
             max |dbeta|/beta {max_beta:.3e}, max |dcost|/cost {max_cost:.3e}",
            cases - identical
        );
    }

    #[test]
    fn c3_cholesky_on_caller_scratch_is_bit_identical() {
        let mut rng = Lcg(0xC401_0C8E_5017_0000);
        // Hoisted, and poisoned up front: the solver is only ever handed
        // scratch that a previous solve (of some other `np`) left dirty, so
        // the pin must exercise that — every element it reads has to be one
        // it wrote in THIS call. A NaN survivor would show up as a mismatch.
        let mut l = [f64::NAN; MAX_LM_PARAMS * MAX_LM_PARAMS];
        let mut y = [f64::NAN; MAX_LM_PARAMS];
        let mut got = [f64::NAN; MAX_LM_PARAMS];
        for np in 2..=MAX_LM_PARAMS {
            for _ in 0..64 {
                // A random SPD matrix: M = R·Rᵀ + n·I.
                let r: Vec<f64> = (0..np * np).map(|_| rng.gauss()).collect();
                let mut mat = vec![0.0_f64; np * np];
                for i in 0..np {
                    for j in 0..np {
                        let mut acc = 0.0;
                        for k in 0..np {
                            acc += r[i * np + k] * r[j * np + k];
                        }
                        mat[i * np + j] = acc + if i == j { np as f64 } else { 0.0 };
                    }
                }
                let rhs: Vec<f64> = (0..np).map(|_| rng.gauss()).collect();

                let want =
                    reference::cholesky_solve_alloc(&mat, &rhs, np).expect("SPD system must solve");
                assert!(cholesky_solve_into(
                    &mat, &rhs, np, &mut l, &mut y, &mut got
                ));
                assert_eq!(&got[..np], &want[..], "np={np}");
            }
        }
    }

    #[test]
    fn c3_integral_exponent_matches_powf_to_a_few_ulp() {
        let mut worst = 0.0_f64;
        // 20 is the free-beta clamp's ceiling and so is reachable; the base
        // runs past the far corner of the widest stamp the fitter builds
        // (radius 48, alpha at its 0.3 floor: 1 + Q reaches ~2.6e4).
        for beta in [2.0_f64, 3.0, 4.0, 6.0, 10.0, 20.0] {
            let e = MoffatExponent::new(beta);
            assert!(matches!(e, MoffatExponent::Integral(_)), "beta {beta}");
            let mut base = 1.0_f64;
            while base < 40_000.0 {
                let got = e.eval(base);
                let want = base.powf(-beta);
                if want > 0.0 {
                    worst = worst.max(((got - want) / want).abs());
                }
                base *= 1.031;
            }
        }
        // A few ulp of f64 — `powi` squares and reciprocates where `powf`
        // takes a log/exp pair.
        assert!(worst < 1e-14, "integral exponent deviation {worst:.3e}");
        println!("C3 powi vs powf: max relative deviation {worst:.3e}");

        // 2.5 is an `AUTO_BETAS` member and must stay on the real path; so
        // must anything outside the 1..=20 window, clamp ceiling included.
        assert!(matches!(MoffatExponent::new(2.5), MoffatExponent::Real(_)));
        assert!(matches!(MoffatExponent::new(0.0), MoffatExponent::Real(_)));
        assert!(matches!(MoffatExponent::new(21.0), MoffatExponent::Real(_)));
        assert!(matches!(MoffatExponent::new(25.0), MoffatExponent::Real(_)));
        assert!(matches!(
            MoffatExponent::new(f64::NAN),
            MoffatExponent::Real(_)
        ));
    }

    // ── The conditioned entry point: the θ column with its correct sign ─────

    /// A noiseless elliptical Moffat (`B = 0`, `A = 1`) on an `n × n` stamp centred on `(cx, cy)`, through the
    /// module's own model.
    fn moffat_stamp(n: usize, cx: f64, cy: f64, ax: f64, ay: f64, theta: f64, beta: f64) -> Vec<PixelSample> {
        let params = [0.0, 1.0, cx, cy, ax, ay, theta, beta];
        let mut px = Vec::with_capacity(n * n);
        for y in 0..n {
            for x in 0..n {
                let (x, y) = (x as f64, y as f64);
                px.push(PixelSample { x, y, value: evaluate_moffat_2d(&params, x, y) });
            }
        }
        px
    }

    fn result_bits(r: &Moffat2DResult) -> [u64; 9] {
        [
            r.b.to_bits(),
            r.a.to_bits(),
            r.x0.to_bits(),
            r.y0.to_bits(),
            r.alpha_x.to_bits(),
            r.alpha_y.to_bits(),
            r.theta.to_bits(),
            r.beta.to_bits(),
            r.fit_residual.to_bits(),
        ]
    }

    /// The θ column the solver builds, `A·(∂(1+Q)^(−β)/∂Q)·∂Q/∂θ` (`a_dpq * column.dq_dtheta(..)` exactly as
    /// `lm_solve_moffat_with` writes it), against a central difference of the model in θ on a small elliptical stamp,
    /// pixel by pixel. The corrected column matches it; the legacy one is its exact negative everywhere.
    #[test]
    fn the_corrected_theta_column_is_the_central_difference_of_the_model() {
        let h = 1e-6;
        let mut checked = 0_usize;
        for &(ax, ay, theta, beta) in &[(2.0_f64, 3.5_f64, 0.6_f64, 4.0_f64), (3.2, 1.7, -1.1, 2.5)] {
            let params = [0.02, 0.9, 4.3, 3.8, ax, ay, theta, beta];
            let (cos_t, sin_t) = (theta.cos(), theta.sin());
            let (inv_ax2, inv_ay2) = (1.0 / (ax * ax), 1.0 / (ay * ay));
            for y in 0..9 {
                for x in 0..9 {
                    let (x, y) = (x as f64, y as f64);
                    let (dx, dy) = (x - params[2], y - params[3]);
                    let u = dx * cos_t + dy * sin_t;
                    let v = -dx * sin_t + dy * cos_t;
                    let base = 1.0 + u * u * inv_ax2 + v * v * inv_ay2;
                    let power = base.powf(-beta);
                    let a_dpq = params[1] * (-beta * (power / base));
                    let corrected = a_dpq * ThetaColumn::Corrected.dq_dtheta(u, v, inv_ax2, inv_ay2);
                    let legacy = a_dpq * ThetaColumn::Legacy.dq_dtheta(u, v, inv_ax2, inv_ay2);

                    let (mut up, mut down) = (params, params);
                    up[6] += h;
                    down[6] -= h;
                    let fd = (evaluate_moffat_2d(&up, x, y) - evaluate_moffat_2d(&down, x, y)) / (2.0 * h);

                    assert!(
                        (corrected - fd).abs() <= 1e-8,
                        "corrected column {corrected:e} vs central difference {fd:e} at ({x}, {y}), params {params:?}"
                    );
                    assert_eq!(legacy.to_bits(), (-corrected).to_bits(), "legacy column at ({x}, {y})");
                    if fd.abs() > 1e-3 {
                        assert!(legacy * fd < 0.0, "the legacy column has the wrong sign at ({x}, {y})");
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 60, "the stamps must carry a real θ gradient: {checked} pixels");
    }

    #[test]
    fn conditioned_convergence_does_not_stop_on_a_theta_step() {
        // a noiseless elliptical Moffat (αx 2.5, αy 4, θ 30°, β 4) sampled on a 31×31 stamp, started with θ off by 5°
        // and everything else exact: the conditioned solver ends within 1° of 30°; the vector-norm solver is recorded
        // (`the_legacy_entry_points_are_output_identical_to_1_3_0`: θ 44.84°, not converged)
        let pixels = moffat_stamp(31, 15.0, 15.0, 2.5, 4.0, 30f64.to_radians(), 4.0);
        let r = fit_moffat_2d_fixed_beta_conditioned(
            &pixels,
            0.0,
            1.0,
            15.0,
            15.0,
            2.5 * 0.6,
            4.0 * 0.6,
            35f64.to_radians(),
            4.0,
            100,
            1e-5,
            5,
        )
        .unwrap();
        assert!((r.theta.to_degrees() - 30.0).abs() < 1.0, "theta {}", r.theta.to_degrees());
    }

    /// A pin, not RED: `fit_moffat_2d_fixed_beta` and `fit_moffat_2d` keep the θ column every release through 1.3.0
    /// shipped, so their output on the input above is bit for bit what 1.3.0 (8321412) gives — recorded there before
    /// the conditioned entry point was added. The numbers are the stall the legacy column causes: from a 35° start the
    /// fit walks AWAY from the true 30°.
    #[test]
    fn the_legacy_entry_points_are_output_identical_to_1_3_0() {
        let pixels = moffat_stamp(31, 15.0, 15.0, 2.5, 4.0, 30f64.to_radians(), 4.0);
        let fixed = fit_moffat_2d_fixed_beta(
            &pixels,
            0.0,
            1.0,
            15.0,
            15.0,
            2.5 * 0.6,
            4.0 * 0.6,
            35f64.to_radians(),
            4.0,
            100,
            1e-5,
            5,
        )
        .unwrap();
        // θ 44.838°, not converged.
        assert_eq!(
            result_bits(&fixed),
            [
                0x3f3597d2f746f108,
                0x3fef2644e43b3e52,
                0x402e000000000000,
                0x402e000000000000,
                0x4003d8f7cfb9d57a,
                0x400f1ab3020d95ce,
                0x3fe90ad307042a7c,
                0x4010000000000000,
                0x3f8177ab7774f2fb,
            ]
        );
        assert!(!fixed.converged);
        let free = fit_moffat_2d(
            &pixels,
            0.0,
            1.0,
            15.0,
            15.0,
            2.5 * 0.6,
            4.0 * 0.6,
            35f64.to_radians(),
            100,
            1e-5,
            5,
        )
        .unwrap();
        // θ 48.396°, β 1.87, not converged.
        assert_eq!(
            result_bits(&free),
            [
                0xbf61565c765f8519,
                0x3fefaedc7e5b15b2,
                0x402e000000000000,
                0x402e000000000000,
                0x3ff8d134924b8a6c,
                0x40039d13198561fa,
                0x3feb078330706328,
                0x3ffdf1d4ebb0f55c,
                0x3f87d78dba90a921,
            ]
        );
        assert!(!free.converged);
    }

    /// Ruling C-4 for the corrected column: the shipped solver and the frozen pre-C3 reference, both on the corrected
    /// column, run the same control flow (iterations, verdict) and land within the C3 tolerances on every stamp of
    /// the C3 grid — the column is threaded through both the same way.
    #[test]
    fn c3_contract_holds_for_the_corrected_theta_column() {
        let fixtures = fixture_stars();
        let modes: [Option<f64>; 5] = [Some(2.5), Some(4.0), Some(6.0), Some(10.0), None];
        let mut cases = 0_usize;
        for f in &fixtures {
            for &fixed_beta in &modes {
                let init = moffat_init(f, fixed_beta);
                let tag = format!("{} fixed_beta={fixed_beta:?}", f.label);
                let mut now = init;
                let r_now = lm_solve_moffat_impl(
                    &f.pixels,
                    &mut now,
                    BetaParam::fixed_or_linear(fixed_beta),
                    CALIBRATION_MAX_ITER,
                    CALIBRATION_CONV_TOL,
                    CALIBRATION_MAX_REJECTS,
                    ThetaColumn::Corrected,
                );
                let mut was = init;
                let r_was = reference::lm_solve_moffat(
                    &f.pixels,
                    &mut was,
                    fixed_beta,
                    CALIBRATION_MAX_ITER,
                    CALIBRATION_CONV_TOL,
                    CALIBRATION_MAX_REJECTS,
                    ThetaColumn::Corrected,
                );
                assert_eq!(r_now.iters, r_was.iters, "iteration count: {tag}");
                assert_eq!(r_now.converged, r_was.converged, "convergence verdict: {tag}");
                let pos = (now[2] - was[2]).abs().max((now[3] - was[3]).abs());
                assert!(pos <= 1e-5, "position {pos:.3e} px: {tag}");
                let (flux_now, flux_was) = (moffat_flux(&now), moffat_flux(&was));
                let flux = if flux_was != 0.0 { ((flux_now - flux_was) / flux_was).abs() } else { 0.0 };
                assert!(flux <= 1e-5, "flux {flux:.3e} relative: {tag}");
                let beta = if was[7] != 0.0 { ((now[7] - was[7]) / was[7]).abs() } else { 0.0 };
                assert!(beta <= 1e-5, "beta {beta:.3e} relative: {tag}");
                cases += 1;
            }
        }
        assert_eq!(cases, 180, "fixture x mode grid");
    }
    // ── The bounded free-β fit (ln β, caller bounds; Athenaeum's per-frame β) ──

    /// A noiseless round Gaussian (`B = 0`, `A = 1`) of width `sigma` on an `n × n` stamp centred on `(cx, cy)`,
    /// point-sampled like [`moffat_stamp`].
    fn gauss_stamp(n: usize, cx: f64, cy: f64, sigma: f64) -> Vec<PixelSample> {
        let s2 = 2.0 * sigma * sigma;
        let mut px = Vec::with_capacity(n * n);
        for y in 0..n {
            for x in 0..n {
                let (x, y) = (x as f64, y as f64);
                let d2 = (x - cx).powi(2) + (y - cy).powi(2);
                px.push(PixelSample { x, y, value: (-d2 / s2).exp() });
            }
        }
        px
    }

    #[test]
    fn bounded_free_beta_reaches_a_gaussian_star_at_the_cap() {
        let pixels = gauss_stamp(41, 20.0, 20.0, 2.0);
        let r = fit_moffat_2d_bounded(&pixels, 0.0, 1.0, 20.0, 20.0, 2.0, 2.0, 0.0, 4.0, (1.5, 100.0), 100, 1e-5, 5).unwrap();
        assert!(r.beta >= 95.0, "beta {}", r.beta);
        let fwhm = 2.0 * r.alpha_x * (2f64.powf(1.0 / r.beta) - 1.0).sqrt();
        assert!((fwhm / (2.0 * 2.354_820_045) - 1.0).abs() < 0.004, "fwhm {fwhm}");
        // the old entry point on the same stamp: None or β at its own 20 cap — the two arms stay distinct
        let old = fit_moffat_2d(&pixels, 0.0, 1.0, 20.0, 20.0, 2.0, 2.0, 0.0, 100, 1e-5, 5);
        assert!(old.as_ref().map_or(true, |o| o.beta >= 19.9), "{:?}", old.as_ref().map(|o| o.beta));
    }

    #[test]
    fn bounded_fit_recovers_beta_4_and_2_5_inside_the_bounds() {
        for (beta, tol) in [(4.0, 0.05), (2.5, 0.03)] {
            let pixels = moffat_stamp(41, 20.0, 20.0, 3.0, 3.0, 0.0, beta);
            let r = fit_moffat_2d_bounded(&pixels, 0.0, 1.0, 20.0, 20.0, 1.5, 1.5, 0.0, 3.0, (1.5, 100.0), 100, 1e-5, 5).unwrap();
            assert!((r.beta - beta).abs() < tol, "beta {} vs {beta}", r.beta);
            assert!(r.iters > 0);
        }
    }

    #[test]
    fn ln_beta_jacobian_matches_central_differences() {
        for beta in [2.0f64, 10.0, 60.0] {
            let (base, a) = (1.7f64, 0.8f64);             // 1 + Q at one pixel, the amplitude
            let analytic = a * base.powf(-beta) * (-base.ln()) * beta;            // ∂M/∂(ln β)
            let h = 1e-6_f64;
            let numeric = (a * base.powf(-beta * (h).exp()) - a * base.powf(-beta * (-h).exp())) / (2.0 * h);
            assert!(((analytic - numeric) / numeric).abs() < 1e-6, "beta {beta}: {analytic} vs {numeric}");
        }
    }

    /// The β column the solver builds (`BetaParam::column` exactly as `lm_solve_moffat_with` calls it) against a
    /// central difference of the module's own model in its `params[7]`, pixel by pixel on an elliptical stamp: in ln β
    /// for `FreeLn`, in β for `FreeLinear`.
    #[test]
    fn the_solvers_beta_columns_are_the_central_differences_of_the_model() {
        let h = 1e-6;
        let ln = BetaParam::FreeLn { ln_min: 1.5f64.ln(), ln_max: 100f64.ln() };
        let mut checked = 0_usize;
        for &beta in &[1.8_f64, 4.0, 37.0] {
            let params = [0.02, 0.9, 4.3, 3.8, 2.6, 3.4, 0.5, beta];
            let (cos_t, sin_t) = (params[6].cos(), params[6].sin());
            for y in 0..9 {
                for x in 0..9 {
                    let (x, y) = (x as f64, y as f64);
                    let (dx, dy) = (x - params[2], y - params[3]);
                    let u = dx * cos_t + dy * sin_t;
                    let v = -dx * sin_t + dy * cos_t;
                    let base = 1.0 + (u / params[4]).powi(2) + (v / params[5]).powi(2);
                    let a_power = params[1] * base.powf(-beta);
                    let (mut up, mut down) = (params, params);
                    up[7] = (beta.ln() + h).exp();
                    down[7] = (beta.ln() - h).exp();
                    let fd_ln = (evaluate_moffat_2d(&up, x, y) - evaluate_moffat_2d(&down, x, y)) / (2.0 * h);
                    let col = ln.column(a_power, base, beta);
                    assert!((col - fd_ln).abs() <= 1e-7 * fd_ln.abs().max(1e-3), "ln column {col:e} vs {fd_ln:e}");
                    let (mut up, mut down) = (params, params);
                    up[7] += h;
                    down[7] -= h;
                    let fd = (evaluate_moffat_2d(&up, x, y) - evaluate_moffat_2d(&down, x, y)) / (2.0 * h);
                    let col = BetaParam::FreeLinear.column(a_power, base, beta);
                    assert!((col - fd).abs() <= 1e-7 * fd.abs().max(1e-3), "linear column {col:e} vs {fd:e}");
                    checked += 1;
                }
            }
        }
        assert_eq!(checked, 243);
    }

    /// A β at a bound is returned, exactly the bound: a star steeper than the lower bound allows (true β 1.2) ends at
    /// 1.5, a Gaussian under a cap of 30 at 30; the bounds are the caller's.
    #[test]
    fn a_fit_at_a_bound_returns_the_bound_itself() {
        let steep = moffat_stamp(41, 20.0, 20.0, 2.0, 2.0, 0.0, 1.2);
        let r = fit_moffat_2d_bounded(&steep, 0.0, 1.0, 20.0, 20.0, 1.0, 1.0, 0.0, 3.0, (1.5, 100.0), 100, 1e-5, 5).unwrap();
        assert_eq!(r.beta, 1.5);
        let gauss = gauss_stamp(41, 20.0, 20.0, 2.0);
        let r = fit_moffat_2d_bounded(&gauss, 0.0, 1.0, 20.0, 20.0, 2.0, 2.0, 0.0, 4.0, (1.5, 30.0), 100, 1e-5, 5).unwrap();
        assert_eq!(r.beta, 30.0);
        assert!(fit_moffat_2d_bounded(&gauss, 0.0, 1.0, 20.0, 20.0, 2.0, 2.0, 0.0, 4.0, (0.0, 30.0), 100, 1e-5, 5).is_none());
        assert!(fit_moffat_2d_bounded(&gauss, 0.0, 1.0, 20.0, 20.0, 2.0, 2.0, 0.0, 4.0, (40.0, 30.0), 100, 1e-5, 5).is_none());
    }

    /// The bounded fit carries the corrected θ column: the elliptical stamp the legacy free fit stalls on
    /// (`the_legacy_entry_points_are_output_identical_to_1_3_0`: θ 48.4°, β 1.87) is recovered at θ 30°, β 4.
    #[test]
    fn the_bounded_fit_recovers_the_elliptical_stamp_the_legacy_free_fit_stalls_on() {
        let pixels = moffat_stamp(31, 15.0, 15.0, 2.5, 4.0, 30f64.to_radians(), 4.0);
        let r = fit_moffat_2d_bounded(&pixels, 0.0, 1.0, 15.0, 15.0, 2.5 * 0.6, 4.0 * 0.6, 35f64.to_radians(), 3.0, (1.5, 100.0), 100, 1e-5, 5).unwrap();
        assert!((r.theta.to_degrees() - 30.0).abs() < 1.0, "theta {}", r.theta.to_degrees());
        assert!((r.beta - 4.0).abs() < 0.05, "beta {}", r.beta);
    }

    /// `fit_moffat_2d`'s free-β output on the C3 stamp, recorded at 197fd93 (= 1.3.0 for this entry point,
    /// `the_legacy_entry_points_are_output_identical_to_1_3_0`) before `BetaParam` replaced the solver's `Option<f64>`.
    const C3_GOLDEN_FREE_BETA: (f64, f64, f64, f64, f64) =
        (2.9776194062739187, 2.9823068067601484, 0.5346392752190088, 2.156163990752592, 0.013907919637681579);

    /// A pin, not RED (ruling F28): `fit_moffat_2d` through `BetaParam::FreeLinear` returns the C3 stamp's exact
    /// parameters and iteration count recorded before this task (the golden values are read from
    /// `c3_moffat_arithmetic_delta_against_the_pre_c3_solver`'s fixture at 8321412 and pasted here).
    #[test]
    fn fit_moffat_2d_is_byte_identical_to_before() {
        let pixels = c3_fixture_stamp();
        let r = fit_moffat_2d(&pixels, 100.0, 5000.0, 15.0, 15.0, 2.0, 2.0, 0.0, CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS).unwrap();
        assert_eq!((r.alpha_x, r.alpha_y, r.theta, r.beta, r.fit_residual), C3_GOLDEN_FREE_BETA);
        // every field, bit for bit, and the solve's control flow (10 iterations, not converged), recorded alongside
        assert_eq!(
            result_bits(&r),
            [
                0x406f49199cadc6bf,
                0x408bc0276d58fd3a,
                0x402df790873d9301,
                0x402e368c97d9fe56,
                0x4007d22a1f8f0d25,
                0x4007dbc3abcd625d,
                0x3fe11bc3d3472110,
                0x40013fd2e808c307,
                0x3f8c7bc15ffe70c8,
            ]
        );
        assert!(!r.converged);
        let scale = (2.0_f64.powf(1.0 / 3.0) - 1.0).sqrt();
        let alpha0 = 2.0 * 2.3548 / (2.0 * scale);
        let mut p = [100.0, 5000.0, 15.0, 15.0, alpha0, alpha0, 0.0, 3.0];
        let lm = lm_solve_moffat_impl(
            &pixels,
            &mut p,
            BetaParam::FreeLinear,
            CALIBRATION_MAX_ITER,
            CALIBRATION_CONV_TOL,
            CALIBRATION_MAX_REJECTS,
            ThetaColumn::Legacy,
        );
        assert_eq!((lm.iters, lm.converged), (10, false));
    }
}
