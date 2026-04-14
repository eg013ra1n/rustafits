/// Gnomonic (TAN) projection: converts between sky coordinates and a tangent plane.
///
/// The tangent point (ra0, dec0) is the projection center. All inputs/outputs
/// are in degrees for RA/Dec and radians for tangent plane coordinates (xi, eta).
pub struct GnomonicProjection;

impl GnomonicProjection {
    /// Project a sky position onto the tangent plane centered at (ra0, dec0).
    ///
    /// Returns (xi, eta) in radians on the tangent plane.
    /// Panics if the point is on the opposite side of the sphere (cos_c <= 0).
    pub fn sky_to_tangent(ra: f64, dec: f64, ra0: f64, dec0: f64) -> (f64, f64) {
        let ra_rad = ra.to_radians();
        let dec_rad = dec.to_radians();
        let ra0_rad = ra0.to_radians();
        let dec0_rad = dec0.to_radians();

        let delta_ra = ra_rad - ra0_rad;
        let sin_dec = dec_rad.sin();
        let cos_dec = dec_rad.cos();
        let sin_dec0 = dec0_rad.sin();
        let cos_dec0 = dec0_rad.cos();

        let cos_c = sin_dec0 * sin_dec + cos_dec0 * cos_dec * delta_ra.cos();
        assert!(
            cos_c > 0.0,
            "Point is on the opposite hemisphere from tangent point"
        );

        let xi = cos_dec * delta_ra.sin() / cos_c;
        let eta = (cos_dec0 * sin_dec - sin_dec0 * cos_dec * delta_ra.cos()) / cos_c;

        (xi, eta)
    }

    /// Inverse gnomonic projection: tangent plane (xi, eta) back to sky (ra, dec).
    ///
    /// (xi, eta) in radians, (ra0, dec0) in degrees. Returns (ra, dec) in degrees.
    pub fn tangent_to_sky(xi: f64, eta: f64, ra0: f64, dec0: f64) -> (f64, f64) {
        let ra0_rad = ra0.to_radians();
        let dec0_rad = dec0.to_radians();

        let sin_dec0 = dec0_rad.sin();
        let cos_dec0 = dec0_rad.cos();

        let rho = (xi * xi + eta * eta).sqrt();

        if rho < 1e-15 {
            return (ra0, dec0);
        }

        let c = rho.atan();
        let sin_c = c.sin();
        let cos_c = c.cos();

        let dec_rad = (cos_c * sin_dec0 + eta * sin_c * cos_dec0 / rho).asin();
        let ra_rad =
            ra0_rad + (xi * sin_c).atan2(rho * cos_dec0 * cos_c - eta * sin_dec0 * sin_c);

        let mut ra_deg = ra_rad.to_degrees();
        if ra_deg < 0.0 {
            ra_deg += 360.0;
        }
        if ra_deg >= 360.0 {
            ra_deg -= 360.0;
        }

        (ra_deg, dec_rad.to_degrees())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_identity() {
        let ra0 = 180.0;
        let dec0 = 45.0;

        let test_points = [
            (180.5, 45.3),
            (179.5, 44.7),
            (181.0, 46.0),
            (178.0, 43.0),
            (180.0, 45.0), // tangent point itself
        ];

        for (ra, dec) in &test_points {
            let (xi, eta) = GnomonicProjection::sky_to_tangent(*ra, *dec, ra0, dec0);
            let (ra_back, dec_back) = GnomonicProjection::tangent_to_sky(xi, eta, ra0, dec0);

            let ra_err = (ra_back - ra).abs();
            let dec_err = (dec_back - dec).abs();
            assert!(
                ra_err < 1e-12,
                "RA roundtrip error {ra_err} for ({ra}, {dec})"
            );
            assert!(
                dec_err < 1e-12,
                "Dec roundtrip error {dec_err} for ({ra}, {dec})"
            );
        }
    }

    #[test]
    fn tangent_point_maps_to_origin() {
        let (xi, eta) = GnomonicProjection::sky_to_tangent(123.456, -45.678, 123.456, -45.678);
        assert!(xi.abs() < 1e-15);
        assert!(eta.abs() < 1e-15);
    }

    #[test]
    fn roundtrip_near_pole() {
        let ra0 = 0.0;
        let dec0 = 89.0;

        let (xi, eta) = GnomonicProjection::sky_to_tangent(10.0, 89.5, ra0, dec0);
        let (ra_back, dec_back) = GnomonicProjection::tangent_to_sky(xi, eta, ra0, dec0);
        assert!((ra_back - 10.0).abs() < 1e-10);
        assert!((dec_back - 89.5).abs() < 1e-10);
    }

    #[test]
    fn roundtrip_ra_wraparound() {
        let ra0 = 359.5;
        let dec0 = 0.0;

        let (xi, eta) = GnomonicProjection::sky_to_tangent(0.5, 0.0, ra0, dec0);
        let (ra_back, dec_back) = GnomonicProjection::tangent_to_sky(xi, eta, ra0, dec0);
        assert!((ra_back - 0.5).abs() < 1e-10);
        assert!((dec_back - 0.0).abs() < 1e-10);
    }
}
