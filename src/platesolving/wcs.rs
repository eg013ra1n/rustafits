use crate::platesolving::projection::GnomonicProjection;
use crate::platesolving::types::SipCoefficients;

/// A complete WCS (World Coordinate System) solution mapping pixels to sky coordinates.
#[derive(Clone, Debug)]
pub struct WcsSolution {
    /// Reference pixel (typically image center).
    pub crpix: (f64, f64),
    /// Sky coordinates at reference pixel (RA, Dec in degrees).
    pub crval: (f64, f64),
    /// CD matrix: linear part of the pixel-to-sky transform.
    /// cd[0][0]=CD1_1, cd[0][1]=CD1_2, cd[1][0]=CD2_1, cd[1][1]=CD2_2
    pub cd: [[f64; 2]; 2],
    /// SIP forward distortion polynomials (A, B). Applied to pixel coords before CD.
    pub sip_forward: Option<(SipCoefficients, SipCoefficients)>,
    /// SIP reverse polynomials (AP, BP). Applied to intermediate coords for sky-to-pixel.
    pub sip_reverse: Option<(SipCoefficients, SipCoefficients)>,
}

impl WcsSolution {
    /// Convert pixel coordinates to sky coordinates (RA, Dec in degrees).
    pub fn pixel_to_sky(&self, x: f64, y: f64) -> (f64, f64) {
        let u = x - self.crpix.0;
        let v = y - self.crpix.1;

        // Apply SIP forward distortion if present
        let (u_corr, v_corr) = if let Some((ref a, ref b)) = self.sip_forward {
            (u + a.evaluate(u, v), v + b.evaluate(u, v))
        } else {
            (u, v)
        };

        // Apply CD matrix to get tangent plane coordinates (degrees)
        let xi_deg = self.cd[0][0] * u_corr + self.cd[0][1] * v_corr;
        let eta_deg = self.cd[1][0] * u_corr + self.cd[1][1] * v_corr;

        // Convert to radians for gnomonic inverse
        let xi_rad = xi_deg.to_radians();
        let eta_rad = eta_deg.to_radians();

        GnomonicProjection::tangent_to_sky(xi_rad, eta_rad, self.crval.0, self.crval.1)
    }

    /// Convert sky coordinates (RA, Dec in degrees) to pixel coordinates.
    pub fn sky_to_pixel(&self, ra: f64, dec: f64) -> (f64, f64) {
        // Project sky to tangent plane
        let (xi_rad, eta_rad) =
            GnomonicProjection::sky_to_tangent(ra, dec, self.crval.0, self.crval.1);
        let xi_deg = xi_rad.to_degrees();
        let eta_deg = eta_rad.to_degrees();

        // Invert CD matrix
        let det = self.cd[0][0] * self.cd[1][1] - self.cd[0][1] * self.cd[1][0];
        if det.abs() < 1e-30 {
            return (f64::NAN, f64::NAN);
        }
        let inv_det = 1.0 / det;
        let u = inv_det * (self.cd[1][1] * xi_deg - self.cd[0][1] * eta_deg);
        let v = inv_det * (-self.cd[1][0] * xi_deg + self.cd[0][0] * eta_deg);

        // Apply SIP reverse correction if present
        let (u_corr, v_corr) = if let Some((ref ap, ref bp)) = self.sip_reverse {
            (u + ap.evaluate(u, v), v + bp.evaluate(u, v))
        } else {
            (u, v)
        };

        (u_corr + self.crpix.0, v_corr + self.crpix.1)
    }

    /// Pixel scale in arcseconds per pixel (geometric mean of axes).
    pub fn pixel_scale_arcsec(&self) -> f64 {
        let scale_x = (self.cd[0][0].powi(2) + self.cd[1][0].powi(2)).sqrt();
        let scale_y = (self.cd[0][1].powi(2) + self.cd[1][1].powi(2)).sqrt();
        (scale_x * scale_y).sqrt() * 3600.0
    }

    /// Field rotation angle in degrees (position angle of the Y axis, N through E).
    pub fn field_rotation_deg(&self) -> f64 {
        // Rotation = atan2(-CD1_2, CD2_2) for standard orientation
        let rot_rad = (-self.cd[0][1]).atan2(self.cd[1][1]);
        let mut deg = rot_rad.to_degrees();
        if deg < 0.0 {
            deg += 360.0;
        }
        deg
    }

    /// Generate FITS WCS header keywords from this solution.
    pub fn to_fits_headers(&self) -> Vec<(String, String)> {
        let mut headers = vec![
            ("WCSAXES".to_string(), "2".to_string()),
            ("CTYPE1".to_string(), "'RA---TAN-SIP'".to_string()),
            ("CTYPE2".to_string(), "'DEC--TAN-SIP'".to_string()),
            (
                "CRPIX1".to_string(),
                format!("{:.6}", self.crpix.0),
            ),
            (
                "CRPIX2".to_string(),
                format!("{:.6}", self.crpix.1),
            ),
            (
                "CRVAL1".to_string(),
                format!("{:.10}", self.crval.0),
            ),
            (
                "CRVAL2".to_string(),
                format!("{:.10}", self.crval.1),
            ),
            (
                "CD1_1".to_string(),
                format!("{:.15E}", self.cd[0][0]),
            ),
            (
                "CD1_2".to_string(),
                format!("{:.15E}", self.cd[0][1]),
            ),
            (
                "CD2_1".to_string(),
                format!("{:.15E}", self.cd[1][0]),
            ),
            (
                "CD2_2".to_string(),
                format!("{:.15E}", self.cd[1][1]),
            ),
        ];

        // Use TAN if no SIP
        if self.sip_forward.is_none() {
            headers[1].1 = "'RA---TAN'".to_string();
            headers[2].1 = "'DEC--TAN'".to_string();
        }

        if let Some((ref a, ref b)) = self.sip_forward {
            headers.push(("A_ORDER".to_string(), format!("{}", a.order)));
            headers.push(("B_ORDER".to_string(), format!("{}", b.order)));

            for i in 0..=(a.order as usize) {
                for j in 0..=(a.order as usize) {
                    if i + j >= 2 && i + j <= a.order as usize && a.coeffs[i][j].abs() > 1e-20 {
                        headers.push((
                            format!("A_{}_{}", i, j),
                            format!("{:.15E}", a.coeffs[i][j]),
                        ));
                    }
                }
            }

            for i in 0..=(b.order as usize) {
                for j in 0..=(b.order as usize) {
                    if i + j >= 2 && i + j <= b.order as usize && b.coeffs[i][j].abs() > 1e-20 {
                        headers.push((
                            format!("B_{}_{}", i, j),
                            format!("{:.15E}", b.coeffs[i][j]),
                        ));
                    }
                }
            }
        }

        if let Some((ref ap, ref bp)) = self.sip_reverse {
            headers.push(("AP_ORDER".to_string(), format!("{}", ap.order)));
            headers.push(("BP_ORDER".to_string(), format!("{}", bp.order)));

            for i in 0..=(ap.order as usize) {
                for j in 0..=(ap.order as usize) {
                    if i + j >= 1 && i + j <= ap.order as usize && ap.coeffs[i][j].abs() > 1e-20 {
                        headers.push((
                            format!("AP_{}_{}", i, j),
                            format!("{:.15E}", ap.coeffs[i][j]),
                        ));
                    }
                }
            }

            for i in 0..=(bp.order as usize) {
                for j in 0..=(bp.order as usize) {
                    if i + j >= 1 && i + j <= bp.order as usize && bp.coeffs[i][j].abs() > 1e-20 {
                        headers.push((
                            format!("BP_{}_{}", i, j),
                            format!("{:.15E}", bp.coeffs[i][j]),
                        ));
                    }
                }
            }
        }

        headers
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_wcs() -> WcsSolution {
        // 1 arcsec/px, no rotation, centered at RA=180 Dec=45
        let scale_deg = 1.0 / 3600.0;
        WcsSolution {
            crpix: (512.5, 512.5),
            crval: (180.0, 45.0),
            cd: [[scale_deg, 0.0], [0.0, scale_deg]],
            sip_forward: None,
            sip_reverse: None,
        }
    }

    #[test]
    fn pixel_to_sky_roundtrip() {
        let wcs = make_simple_wcs();

        let test_pixels = [
            (512.5, 512.5), // center
            (600.0, 512.5),
            (512.5, 400.0),
            (300.0, 700.0),
        ];

        for (x, y) in &test_pixels {
            let (ra, dec) = wcs.pixel_to_sky(*x, *y);
            let (x_back, y_back) = wcs.sky_to_pixel(ra, dec);
            assert!(
                (x_back - x).abs() < 0.01,
                "X roundtrip: {x} → {x_back}"
            );
            assert!(
                (y_back - y).abs() < 0.01,
                "Y roundtrip: {y} → {y_back}"
            );
        }
    }

    #[test]
    fn center_pixel_maps_to_crval() {
        let wcs = make_simple_wcs();
        let (ra, dec) = wcs.pixel_to_sky(512.5, 512.5);
        assert!((ra - 180.0).abs() < 1e-10);
        assert!((dec - 45.0).abs() < 1e-10);
    }

    #[test]
    fn pixel_scale() {
        let wcs = make_simple_wcs();
        let scale = wcs.pixel_scale_arcsec();
        assert!(
            (scale - 1.0).abs() < 0.01,
            "Expected ~1 arcsec/px, got {scale}"
        );
    }

    #[test]
    fn fits_headers_produced() {
        let wcs = make_simple_wcs();
        let headers = wcs.to_fits_headers();

        let keys: Vec<&str> = headers.iter().map(|(k, _)| k.as_str()).collect();
        assert!(keys.contains(&"CRPIX1"));
        assert!(keys.contains(&"CRVAL1"));
        assert!(keys.contains(&"CD1_1"));
        assert!(keys.contains(&"CTYPE1"));
    }
}
