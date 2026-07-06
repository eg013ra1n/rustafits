//! Diagnostic probe: eccentricity + position-angle recovery accuracy.
//!
//! Sweeps synthetic elliptical Moffat stars over (size × ecc × theta × SNR),
//! runs the production measurement path (`measure_stars` with fixed field
//! beta, the same call shape `run_analysis` uses), and prints recovery bias
//! tables. Ground truth is exact, so this separates estimator-intrinsic bias
//! from population/statistics effects seen in PixInsight comparisons.
//!
//! Ignored by default — it prints tables rather than asserting:
//!     cargo test --release --features debug-pipeline --test ecc_theta_probe -- --ignored --nocapture
#![cfg(feature = "debug-pipeline")]

use astroimage::analysis::detection::DetectedStar;
use astroimage::analysis::metrics;

const W: usize = 96;
const H: usize = 96;
const BG: f32 = 500.0;
const NOISE_SIGMA: f32 = 16.0;
const BETA: f32 = 3.0;
const REALIZATIONS: usize = 24;

/// Deterministic LCG → Box-Muller Gaussian noise (no rand dependency).
struct Lcg(u64);
impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 32) as u32
    }
    fn uniform(&mut self) -> f32 {
        (self.next_u32() as f32 + 0.5) / (u32::MAX as f32 + 1.0)
    }
    fn gaussian(&mut self) -> f32 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f32::consts::PI * u2).cos()
    }
}

/// Elliptical Moffat, 4×4 sub-pixel integrated (models real detector sampling).
fn render_star(
    data: &mut [f32],
    cx: f32,
    cy: f32,
    amp: f32,
    alpha_x: f32,
    alpha_y: f32,
    theta: f32,
) {
    let (ct, st) = (theta.cos(), theta.sin());
    for y in 0..H {
        for x in 0..W {
            let mut acc = 0.0_f32;
            for sy in 0..4 {
                for sx in 0..4 {
                    let px = x as f32 + (sx as f32 + 0.5) / 4.0 - 0.5;
                    let py = y as f32 + (sy as f32 + 0.5) / 4.0 - 0.5;
                    let dx = px - cx;
                    let dy = py - cy;
                    let u = dx * ct + dy * st;
                    let v = -dx * st + dy * ct;
                    let r2 = (u / alpha_x).powi(2) + (v / alpha_y).powi(2);
                    acc += amp * (1.0 + r2).powf(-BETA);
                }
            }
            data[y * W + x] += acc / 16.0;
        }
    }
}

/// FWHM of a Moffat along one axis for α: FWHM = 2α√(2^(1/β) − 1).
fn moffat_fwhm(alpha: f32) -> f32 {
    2.0 * alpha * (2.0_f32.powf(1.0 / BETA) - 1.0).sqrt()
}

/// Circular angle difference mod π, in degrees.
fn theta_err_deg(measured: f32, truth: f32) -> f32 {
    let mut d = (measured - truth) % std::f32::consts::PI;
    if d > std::f32::consts::FRAC_PI_2 {
        d -= std::f32::consts::PI;
    } else if d < -std::f32::consts::FRAC_PI_2 {
        d += std::f32::consts::PI;
    }
    d.to_degrees()
}

#[test]
#[ignore]
fn ecc_theta_recovery_sweep() {
    // (label, alpha_major) — alpha 1.7 → FWHM ≈ 2.6px (tight), 3.9 → 6px, 7.8 → 12px (defocus-ish)
    let sizes: &[(&str, f32)] = &[("FWHM~2.6", 1.7), ("FWHM~6.0", 3.9), ("FWHM~12", 7.8)];
    let eccs: &[f32] = &[0.0, 0.2, 0.4, 0.6, 0.8];
    let thetas_deg: &[f32] = &[0.0, 30.0, 45.0, 60.0, 90.0, 120.0, 150.0];
    // amp 20000 → per-pixel SNR >> 1 (estimator-intrinsic bias); amp 300 → faint regime
    let amps: &[(&str, f32)] = &[("bright", 20000.0), ("faint", 300.0)];

    println!();
    println!(
        "{:<10} {:<7} {:>5} {:>7} | {:>9} {:>9} | {:>9} {:>9} {:>5}",
        "size", "amp", "ecc", "θ(deg)", "ecc_meas", "ecc_bias", "θerr_mean", "θerr_rms", "n"
    );
    println!("{}", "-".repeat(88));

    for &(size_label, alpha_major) in sizes {
        for &(amp_label, amp) in amps {
            for &true_ecc in eccs {
                // axis ratio from ecc: b/a = sqrt(1 − e²); ecc from FWHM ratio
                // is identical to alpha ratio at fixed beta.
                let ratio = (1.0 - true_ecc * true_ecc).sqrt();
                let alpha_x = alpha_major;
                let alpha_y = alpha_major * ratio;
                let true_fwhm_major = moffat_fwhm(alpha_x);
                let true_fwhm_geo = moffat_fwhm(alpha_x) * ratio.sqrt();

                // Aggregate over angles × realizations; per-angle rows only
                // matter when hunting an angle-dependent bug, so we print a
                // summary row per (size, amp, ecc) plus per-angle θ stats.
                let mut ecc_sum = 0.0_f64;
                let mut fwhm_ratio_sum = 0.0_f64;
                let mut n_ok = 0usize;
                let mut theta_errs: Vec<f32> = Vec::new();

                for &th_deg in thetas_deg {
                    let true_theta = th_deg.to_radians();
                    for r in 0..REALIZATIONS {
                        let mut rng = Lcg(
                            0x9E3779B97F4A7C15
                                ^ ((th_deg as u64) << 32)
                                ^ ((true_ecc * 100.0) as u64) << 16
                                ^ ((alpha_major * 10.0) as u64) << 8
                                ^ (amp as u64)
                                ^ r as u64,
                        );
                        let cx = 48.0 + rng.uniform() - 0.5;
                        let cy = 48.0 + rng.uniform() - 0.5;

                        let mut data = vec![BG; W * H];
                        render_star(&mut data, cx, cy, amp, alpha_x, alpha_y, true_theta);
                        for v in data.iter_mut() {
                            *v += NOISE_SIGMA * rng.gaussian();
                        }

                        let star = DetectedStar {
                            x: cx,
                            y: cy,
                            peak: amp,
                            flux: amp * alpha_x * alpha_y * 3.0,
                            area: (true_fwhm_geo * true_fwhm_geo * 2.0) as usize + 5,
                            theta: 0.0,
                            eccentricity: 0.0,
                        };

                        let measured = metrics::measure_stars(
                            &data, W, H,
                            std::slice::from_ref(&star),
                            BG, None, None,
                            Some(BETA as f64),
                            25, 1e-4, 5,
                            Some(true_fwhm_geo),
                        );

                        if let Some(m) = measured.first() {
                            ecc_sum += m.eccentricity as f64;
                            fwhm_ratio_sum += (m.fwhm / true_fwhm_geo) as f64;
                            n_ok += 1;
                            if true_ecc >= 0.15 {
                                theta_errs.push(theta_err_deg(m.theta, true_theta));
                            }
                        }
                    }
                }

                let total = thetas_deg.len() * REALIZATIONS;
                let ecc_mean = if n_ok > 0 { (ecc_sum / n_ok as f64) as f32 } else { f32::NAN };
                let fwhm_pct = if n_ok > 0 {
                    ((fwhm_ratio_sum / n_ok as f64 - 1.0) * 100.0) as f32
                } else {
                    f32::NAN
                };
                let (te_mean, te_rms) = if theta_errs.is_empty() {
                    (f32::NAN, f32::NAN)
                } else {
                    let mean = theta_errs.iter().sum::<f32>() / theta_errs.len() as f32;
                    let rms = (theta_errs.iter().map(|e| e * e).sum::<f32>()
                        / theta_errs.len() as f32)
                        .sqrt();
                    (mean, rms)
                };

                println!(
                    "{:<10} {:<7} {:>5.2} {:>7} | {:>9.3} {:>+9.3} | {:>+8.2}% | {:>9.2} {:>9.2} {:>3}/{:<3}",
                    size_label, amp_label, true_ecc, "all",
                    ecc_mean, ecc_mean - true_ecc, fwhm_pct,
                    te_mean, te_rms, n_ok, total,
                );
                let _ = true_fwhm_major;
            }
        }
        println!("{}", "-".repeat(88));
    }
}
