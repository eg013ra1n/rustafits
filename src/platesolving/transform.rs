use nalgebra::{DMatrix, DVector};

use crate::platesolving::projection::GnomonicProjection;
use crate::platesolving::types::{FitModel, SipCoefficients};
use crate::platesolving::wcs::WcsSolution;

/// Fits coordinate transforms from matched star pairs.
pub struct TransformFitter;

/// Result of a pixel-to-pixel transform fit (for registration, no sky coords).
#[derive(Clone, Debug)]
pub struct PixelTransform {
    pub model: FitModel,
    pub params: Vec<f64>,
}

impl PixelTransform {
    /// Transform a source pixel to target pixel coordinates.
    pub fn transform(&self, x: f64, y: f64) -> (f64, f64) {
        match self.model {
            FitModel::Affine => {
                // params: [a, b, tx, c, d, ty]
                let xp = self.params[0] * x + self.params[1] * y + self.params[2];
                let yp = self.params[3] * x + self.params[4] * y + self.params[5];
                (xp, yp)
            }
            FitModel::Projective => {
                // params: [h11, h12, h13, h21, h22, h23, h31, h32] (h33=1)
                let w = self.params[6] * x + self.params[7] * y + 1.0;
                let xp = (self.params[0] * x + self.params[1] * y + self.params[2]) / w;
                let yp = (self.params[3] * x + self.params[4] * y + self.params[5]) / w;
                (xp, yp)
            }
            FitModel::Sip { .. } => {
                // SIP only applies to WCS, not pixel transforms
                (x, y)
            }
        }
    }
}

impl TransformFitter {
    /// Fit a WCS solution from matched pixel↔sky pairs.
    ///
    /// `pairs`: (pixel_x, pixel_y) → (ra_deg, dec_deg)
    /// `model`: which transform model to fit
    /// `image_center`: (cx, cy) pixel coords of image center (used as CRPIX)
    pub fn fit_wcs(
        pairs: &[((f64, f64), (f64, f64))],
        model: FitModel,
        image_center: (f64, f64),
    ) -> Result<WcsSolution, String> {
        if pairs.len() < 3 {
            return Err("Need at least 3 matched pairs for WCS fit".to_string());
        }

        // Estimate field center from the median of sky positions
        let mut ras: Vec<f64> = pairs.iter().map(|(_, (ra, _))| *ra).collect();
        let mut decs: Vec<f64> = pairs.iter().map(|(_, (_, dec))| *dec).collect();
        ras.sort_by(|a, b| a.partial_cmp(b).unwrap());
        decs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let ra0 = ras[ras.len() / 2];
        let dec0 = decs[decs.len() / 2];

        // Project all sky positions onto tangent plane centered at (ra0, dec0)
        let projected: Vec<((f64, f64), (f64, f64))> = pairs
            .iter()
            .map(|((px, py), (ra, dec))| {
                let (xi, eta) = GnomonicProjection::sky_to_tangent(*ra, *dec, ra0, dec0);
                // Convert tangent plane from radians to degrees for CD matrix
                ((*px, *py), (xi.to_degrees(), eta.to_degrees()))
            })
            .collect();

        match model {
            FitModel::Affine | FitModel::Projective => {
                fit_affine_wcs(&projected, image_center, ra0, dec0)
            }
            FitModel::Sip { order } => {
                fit_sip_wcs(&projected, image_center, ra0, dec0, order)
            }
        }
    }

    /// Fit a pixel-to-pixel transform from matched pairs (for registration).
    pub fn fit_pixel(
        pairs: &[((f64, f64), (f64, f64))],
        model: FitModel,
    ) -> Result<PixelTransform, String> {
        match model {
            FitModel::Affine => fit_affine_pixel(pairs),
            FitModel::Projective => fit_projective_pixel(pairs),
            FitModel::Sip { .. } => {
                Err("SIP model only applies to WCS fitting, not pixel transforms".to_string())
            }
        }
    }
}

/// Fit an affine WCS: solve for CD matrix using least squares.
///
/// Model: xi_deg  = CD1_1 * u + CD1_2 * v
///        eta_deg = CD2_1 * u + CD2_2 * v
/// where u = x - crpix_x, v = y - crpix_y
fn fit_affine_wcs(
    pairs: &[((f64, f64), (f64, f64))],
    image_center: (f64, f64),
    ra0: f64,
    dec0: f64,
) -> Result<WcsSolution, String> {
    let n = pairs.len();

    // Build system: A * [CD1_1, CD1_2, CD2_1, CD2_2]^T = b
    // Two equations per pair:
    //   u * CD1_1 + v * CD1_2 = xi_deg
    //   u * CD2_1 + v * CD2_2 = eta_deg
    let mut a_mat = DMatrix::zeros(2 * n, 4);
    let mut b_vec = DVector::zeros(2 * n);

    for (i, ((px, py), (xi_deg, eta_deg))) in pairs.iter().enumerate() {
        let u = px - image_center.0;
        let v = py - image_center.1;

        a_mat[(2 * i, 0)] = u;
        a_mat[(2 * i, 1)] = v;
        b_vec[2 * i] = *xi_deg;

        a_mat[(2 * i + 1, 2)] = u;
        a_mat[(2 * i + 1, 3)] = v;
        b_vec[2 * i + 1] = *eta_deg;
    }

    // Solve via normal equations: A^T A x = A^T b
    let ata = a_mat.transpose() * &a_mat;
    let atb = a_mat.transpose() * &b_vec;

    let decomp = ata.lu();
    let x = decomp
        .solve(&atb)
        .ok_or("Singular matrix in WCS fit")?;

    Ok(WcsSolution {
        crpix: image_center,
        crval: (ra0, dec0),
        cd: [[x[0], x[1]], [x[2], x[3]]],
        sip_forward: None,
        sip_reverse: None,
    })
}

/// Fit a WCS with SIP distortion correction.
///
/// 1. Fit the linear CD matrix (affine)
/// 2. Compute residuals in pixel space (difference between where CD predicts and actual)
/// 3. Fit SIP polynomial A/B coefficients to these residuals
/// 4. Compute reverse SIP polynomials AP/BP
fn fit_sip_wcs(
    pairs: &[((f64, f64), (f64, f64))],
    image_center: (f64, f64),
    ra0: f64,
    dec0: f64,
    order: u8,
) -> Result<WcsSolution, String> {
    let order = order.clamp(2, 5) as usize;

    // Step 1: Fit the linear CD matrix
    let linear = fit_affine_wcs(pairs, image_center, ra0, dec0)?;
    let cd = linear.cd;

    // Invert CD matrix for sky→pixel direction
    let det = cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0];
    if det.abs() < 1e-30 {
        return Err("Singular CD matrix, cannot fit SIP".to_string());
    }
    let cd_inv = [
        [cd[1][1] / det, -cd[0][1] / det],
        [-cd[1][0] / det, cd[0][0] / det],
    ];

    // Step 2: SIP convention: xi_deg = CD * (u + A(u,v), v + B(u,v))
    // where u = x - crpix (detected/distorted pixel offset).
    // So: (u + A(u,v)) = CD^-1 * xi_deg = u_undistorted
    // Therefore: A(u,v) = u_undistorted - u
    //
    // The polynomial is a function of the detected pixel offsets (u, v).

    let n = pairs.len();
    let mut u_detected = Vec::with_capacity(n);
    let mut v_detected = Vec::with_capacity(n);
    let mut du = Vec::with_capacity(n);
    let mut dv = Vec::with_capacity(n);

    for ((px, py), (xi_deg, eta_deg)) in pairs {
        let u = px - image_center.0;
        let v = py - image_center.1;

        let u_undist = cd_inv[0][0] * xi_deg + cd_inv[0][1] * eta_deg;
        let v_undist = cd_inv[1][0] * xi_deg + cd_inv[1][1] * eta_deg;

        u_detected.push(u);
        v_detected.push(v);
        du.push(u_undist - u);
        dv.push(v_undist - v);
    }

    // Step 3: Fit SIP polynomials A and B
    let num_terms = count_sip_terms(order);

    if n < num_terms + 4 {
        return Ok(linear);
    }

    let sip_a = fit_sip_polynomial(&u_detected, &v_detected, &du, order)?;
    let sip_b = fit_sip_polynomial(&u_detected, &v_detected, &dv, order)?;

    // Step 4: Compute reverse SIP polynomials AP/BP
    // For sky_to_pixel: CD^-1 gives u_undistorted, then we need to find u_detected.
    // u_detected = u_undistorted + AP(u_undistorted, v_undistorted)
    // So: AP(u_undist, v_undist) = u_detected - u_undistorted = -A(u_detected, v_detected)
    // We approximate by fitting a polynomial on the undistorted coords.
    let mut u_undist_all = Vec::with_capacity(n);
    let mut v_undist_all = Vec::with_capacity(n);
    let mut du_rev = Vec::with_capacity(n);
    let mut dv_rev = Vec::with_capacity(n);

    for i in 0..n {
        let u_undist = u_detected[i] + du[i];
        let v_undist = v_detected[i] + dv[i];
        u_undist_all.push(u_undist);
        v_undist_all.push(v_undist);
        du_rev.push(u_detected[i] - u_undist);
        dv_rev.push(v_detected[i] - v_undist);
    }

    let sip_ap = fit_sip_polynomial(&u_undist_all, &v_undist_all, &du_rev, order)
        .unwrap_or_else(|_| SipCoefficients::new(order as u8));
    let sip_bp = fit_sip_polynomial(&u_undist_all, &v_undist_all, &dv_rev, order)
        .unwrap_or_else(|_| SipCoefficients::new(order as u8));

    Ok(WcsSolution {
        crpix: image_center,
        crval: (ra0, dec0),
        cd,
        sip_forward: Some((sip_a, sip_b)),
        sip_reverse: Some((sip_ap, sip_bp)),
    })
}

/// Count the number of SIP polynomial terms for a given order.
/// Terms are those where 2 <= i+j <= order.
fn count_sip_terms(order: usize) -> usize {
    let mut count = 0;
    for i in 0..=order {
        for j in 0..=order {
            if i + j >= 2 && i + j <= order {
                count += 1;
            }
        }
    }
    count
}

/// Fit a single SIP polynomial (A or B) via least squares.
///
/// Model: correction(u, v) = sum_{2 <= i+j <= order} c_ij * u^i * v^j
fn fit_sip_polynomial(
    u: &[f64],
    v: &[f64],
    corrections: &[f64],
    order: usize,
) -> Result<SipCoefficients, String> {
    let n = u.len();
    let num_terms = count_sip_terms(order);

    if n < num_terms {
        return Err(format!(
            "Need at least {num_terms} points for SIP order {order}, have {n}"
        ));
    }

    // Build mapping from flat index to (i, j) pairs
    let mut term_indices: Vec<(usize, usize)> = Vec::with_capacity(num_terms);
    for i in 0..=order {
        for j in 0..=order {
            if i + j >= 2 && i + j <= order {
                term_indices.push((i, j));
            }
        }
    }

    // Build design matrix
    let mut a_mat = DMatrix::zeros(n, num_terms);
    let mut b_vec = DVector::zeros(n);

    for row in 0..n {
        for (col, &(i, j)) in term_indices.iter().enumerate() {
            a_mat[(row, col)] = u[row].powi(i as i32) * v[row].powi(j as i32);
        }
        b_vec[row] = corrections[row];
    }

    let ata = a_mat.transpose() * &a_mat;
    let atb = a_mat.transpose() * &b_vec;

    let decomp = ata.lu();
    let x = decomp
        .solve(&atb)
        .ok_or("Singular matrix in SIP polynomial fit")?;

    let mut coeffs = SipCoefficients::new(order as u8);
    for (col, &(i, j)) in term_indices.iter().enumerate() {
        coeffs.coeffs[i][j] = x[col];
    }

    Ok(coeffs)
}

/// Fit an affine pixel-to-pixel transform.
fn fit_affine_pixel(
    pairs: &[((f64, f64), (f64, f64))],
) -> Result<PixelTransform, String> {
    if pairs.len() < 3 {
        return Err("Need at least 3 pairs for affine fit".to_string());
    }

    let n = pairs.len();
    let mut a_mat = DMatrix::zeros(2 * n, 6);
    let mut b_vec = DVector::zeros(2 * n);

    for (i, ((sx, sy), (tx, ty))) in pairs.iter().enumerate() {
        a_mat[(2 * i, 0)] = *sx;
        a_mat[(2 * i, 1)] = *sy;
        a_mat[(2 * i, 2)] = 1.0;
        b_vec[2 * i] = *tx;

        a_mat[(2 * i + 1, 3)] = *sx;
        a_mat[(2 * i + 1, 4)] = *sy;
        a_mat[(2 * i + 1, 5)] = 1.0;
        b_vec[2 * i + 1] = *ty;
    }

    let ata = a_mat.transpose() * &a_mat;
    let atb = a_mat.transpose() * &b_vec;

    let decomp = ata.lu();
    let x = decomp
        .solve(&atb)
        .ok_or("Singular matrix in affine pixel fit")?;

    Ok(PixelTransform {
        model: FitModel::Affine,
        params: x.iter().copied().collect(),
    })
}

/// Fit a projective (homography) pixel-to-pixel transform.
fn fit_projective_pixel(
    pairs: &[((f64, f64), (f64, f64))],
) -> Result<PixelTransform, String> {
    if pairs.len() < 4 {
        return Err("Need at least 4 pairs for projective fit".to_string());
    }

    let n = pairs.len();
    // Each pair gives 2 equations in 8 unknowns (h33 = 1):
    // tx = (h11*sx + h12*sy + h13) / (h31*sx + h32*sy + 1)
    // ty = (h21*sx + h22*sy + h23) / (h31*sx + h32*sy + 1)
    //
    // Rearranged:
    // h11*sx + h12*sy + h13 - h31*sx*tx - h32*sy*tx = tx
    // h21*sx + h22*sy + h23 - h31*sx*ty - h32*sy*ty = ty
    let mut a_mat = DMatrix::zeros(2 * n, 8);
    let mut b_vec = DVector::zeros(2 * n);

    for (i, ((sx, sy), (tx, ty))) in pairs.iter().enumerate() {
        a_mat[(2 * i, 0)] = *sx;
        a_mat[(2 * i, 1)] = *sy;
        a_mat[(2 * i, 2)] = 1.0;
        a_mat[(2 * i, 6)] = -sx * tx;
        a_mat[(2 * i, 7)] = -sy * tx;
        b_vec[2 * i] = *tx;

        a_mat[(2 * i + 1, 3)] = *sx;
        a_mat[(2 * i + 1, 4)] = *sy;
        a_mat[(2 * i + 1, 5)] = 1.0;
        a_mat[(2 * i + 1, 6)] = -sx * ty;
        a_mat[(2 * i + 1, 7)] = -sy * ty;
        b_vec[2 * i + 1] = *ty;
    }

    let ata = a_mat.transpose() * &a_mat;
    let atb = a_mat.transpose() * &b_vec;

    let decomp = ata.lu();
    let x = decomp
        .solve(&atb)
        .ok_or("Singular matrix in projective fit")?;

    Ok(PixelTransform {
        model: FitModel::Projective,
        params: x.iter().copied().collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn affine_wcs_fit_simple() {
        // Generate synthetic star pairs with a known WCS
        let ra0 = 180.0;
        let dec0 = 45.0;
        let pixel_scale_deg = 1.0 / 3600.0; // 1 arcsec/px
        let image_center = (512.5, 512.5);

        // Ground truth CD matrix (no rotation)
        let cd = [[pixel_scale_deg, 0.0], [0.0, pixel_scale_deg]];

        // Generate pairs
        let mut pairs = Vec::new();
        let offsets = [
            (-200.0, -200.0),
            (200.0, -200.0),
            (-200.0, 200.0),
            (200.0, 200.0),
            (0.0, 0.0),
            (100.0, -100.0),
            (-150.0, 50.0),
        ];

        for (dx, dy) in &offsets {
            let px = image_center.0 + dx;
            let py = image_center.1 + dy;
            let u = px - image_center.0;
            let v = py - image_center.1;
            let xi_deg: f64 = cd[0][0] * u + cd[0][1] * v;
            let eta_deg: f64 = cd[1][0] * u + cd[1][1] * v;
            let (ra, dec) = GnomonicProjection::tangent_to_sky(
                xi_deg.to_radians(),
                eta_deg.to_radians(),
                ra0,
                dec0,
            );
            pairs.push(((px, py), (ra, dec)));
        }

        let wcs =
            TransformFitter::fit_wcs(&pairs, FitModel::Affine, image_center).unwrap();

        // Check recovered CD matrix
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (wcs.cd[i][j] - cd[i][j]).abs() < 1e-12,
                    "CD[{i}][{j}]: expected {}, got {}",
                    cd[i][j],
                    wcs.cd[i][j]
                );
            }
        }

        // Check pixel_to_sky roundtrip
        let (ra, dec) = wcs.pixel_to_sky(600.0, 400.0);
        let (x_back, y_back) = wcs.sky_to_pixel(ra, dec);
        assert!((x_back - 600.0).abs() < 0.01);
        assert!((y_back - 400.0).abs() < 0.01);
    }

    #[test]
    fn affine_wcs_fit_rotated() {
        let ra0 = 83.6;
        let dec0 = -5.4;
        let pixel_scale_deg = 1.5 / 3600.0;
        let image_center = (2048.0, 2048.0);

        // 30 degree rotation
        let angle = 30.0f64.to_radians();
        let cd = [
            [
                pixel_scale_deg * angle.cos(),
                -pixel_scale_deg * angle.sin(),
            ],
            [
                pixel_scale_deg * angle.sin(),
                pixel_scale_deg * angle.cos(),
            ],
        ];

        let mut pairs = Vec::new();
        for dx in [-500.0, 0.0, 500.0] {
            for dy in [-500.0, 0.0, 500.0] {
                let px = image_center.0 + dx;
                let py = image_center.1 + dy;
                let u = px - image_center.0;
                let v = py - image_center.1;
                let xi_deg = cd[0][0] * u + cd[0][1] * v;
                let eta_deg = cd[1][0] * u + cd[1][1] * v;
                let (ra, dec) = GnomonicProjection::tangent_to_sky(
                    xi_deg.to_radians(),
                    eta_deg.to_radians(),
                    ra0,
                    dec0,
                );
                pairs.push(((px, py), (ra, dec)));
            }
        }

        let wcs =
            TransformFitter::fit_wcs(&pairs, FitModel::Affine, image_center).unwrap();

        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (wcs.cd[i][j] - cd[i][j]).abs() < 1e-10,
                    "CD[{i}][{j}]: expected {}, got {}",
                    cd[i][j],
                    wcs.cd[i][j]
                );
            }
        }
    }

    #[test]
    fn affine_pixel_fit() {
        // Known affine: x' = 2x + 0.1y + 10, y' = -0.1x + 2y + 20
        let pairs: Vec<((f64, f64), (f64, f64))> = vec![
            ((0.0, 0.0), (10.0, 20.0)),
            ((100.0, 0.0), (210.0, 10.0)),
            ((0.0, 100.0), (20.0, 220.0)),
            ((100.0, 100.0), (220.0, 210.0)),
            ((50.0, 50.0), (115.0, 115.0)),
        ];

        let t = TransformFitter::fit_pixel(&pairs, FitModel::Affine).unwrap();

        let (xp, yp) = t.transform(75.0, 25.0);
        let expected_x = 2.0 * 75.0 + 0.1 * 25.0 + 10.0;
        let expected_y = -0.1 * 75.0 + 2.0 * 25.0 + 20.0;
        assert!((xp - expected_x).abs() < 0.01);
        assert!((yp - expected_y).abs() < 0.01);
    }

    #[test]
    fn sip_reduces_distortion_residuals() {
        // Synthetic field with realistic barrel distortion.
        // Setup: true sky positions → linear CD → undistorted pixels → apply distortion → detected pixels.
        // The fitter receives (detected_pixel, sky_position) pairs and should recover the distortion.
        let ra0 = 180.0;
        let dec0 = 45.0;
        let pixel_scale_deg = 2.0 / 3600.0; // 2 arcsec/px
        let image_center = (2048.0, 2048.0);
        let cd = [[pixel_scale_deg, 0.0], [0.0, pixel_scale_deg]];

        // Realistic barrel distortion: ~0.5% at corners of a 4096x4096 sensor
        let k = 3e-10;

        // Generate a grid of sky positions and their distorted pixel positions
        let mut pairs = Vec::new();
        for dx_step in -5..=5 {
            for dy_step in -5..=5 {
                // Undistorted pixel offsets
                let u_true = dx_step as f64 * 350.0;
                let v_true = dy_step as f64 * 350.0;

                // Sky position from undistorted pixels via CD (ground truth)
                let xi_deg = cd[0][0] * u_true + cd[0][1] * v_true;
                let eta_deg = cd[1][0] * u_true + cd[1][1] * v_true;
                let (ra, dec) = GnomonicProjection::tangent_to_sky(
                    xi_deg.to_radians(),
                    eta_deg.to_radians(),
                    ra0, dec0,
                );

                // Apply barrel distortion to get detected pixel positions
                let r2 = u_true * u_true + v_true * v_true;
                let u_detected = u_true * (1.0 + k * r2);
                let v_detected = v_true * (1.0 + k * r2);
                let px = image_center.0 + u_detected;
                let py = image_center.1 + v_detected;

                pairs.push(((px, py), (ra, dec)));
            }
        }

        let wcs_affine =
            TransformFitter::fit_wcs(&pairs, FitModel::Affine, image_center).unwrap();
        let wcs_sip =
            TransformFitter::fit_wcs(&pairs, FitModel::Sip { order: 3 }, image_center).unwrap();

        // Evaluate at corners: compute pixel_to_sky error vs ground truth
        let test_offsets = [(-1800.0, -1800.0), (1800.0, -1800.0), (-1800.0, 1800.0), (1800.0, 1800.0)];

        let mut rms_affine = 0.0;
        let mut rms_sip = 0.0;

        for &(u_true, v_true) in &test_offsets {
            let xi_deg = cd[0][0] * u_true + cd[0][1] * v_true;
            let eta_deg = cd[1][0] * u_true + cd[1][1] * v_true;
            let (ra_true, dec_true) = GnomonicProjection::tangent_to_sky(
                xi_deg.to_radians(), eta_deg.to_radians(), ra0, dec0,
            );

            let r2 = u_true * u_true + v_true * v_true;
            let px = image_center.0 + u_true * (1.0 + k * r2);
            let py = image_center.1 + v_true * (1.0 + k * r2);

            let (ra_a, dec_a) = wcs_affine.pixel_to_sky(px, py);
            let err_a = ((ra_a - ra_true).powi(2) + (dec_a - dec_true).powi(2)).sqrt() * 3600.0;
            rms_affine += err_a * err_a;

            let (ra_s, dec_s) = wcs_sip.pixel_to_sky(px, py);
            let err_s = ((ra_s - ra_true).powi(2) + (dec_s - dec_true).powi(2)).sqrt() * 3600.0;
            rms_sip += err_s * err_s;
        }

        rms_affine = (rms_affine / 4.0).sqrt();
        rms_sip = (rms_sip / 4.0).sqrt();

        assert!(
            wcs_sip.sip_forward.is_some(),
            "SIP forward polynomials should be present"
        );
        assert!(
            rms_sip < rms_affine,
            "SIP should reduce corner residuals: affine={rms_affine:.4}\" vs sip={rms_sip:.4}\""
        );
    }
}
