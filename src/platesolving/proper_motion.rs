/// Proper motion correction: propagate star positions across epochs.
pub struct ProperMotionCorrector;

impl ProperMotionCorrector {
    /// Propagate a star position from one epoch to another using linear proper motion.
    ///
    /// - `ra`, `dec`: position in degrees at `epoch_from`
    /// - `pmra_mas_yr`: proper motion in RA in milliarcseconds/year (includes cos(dec) factor)
    /// - `pmdec_mas_yr`: proper motion in Dec in milliarcseconds/year
    /// - `epoch_from`, `epoch_to`: Julian epoch (e.g. 2000.0, 2016.0, 2025.5)
    ///
    /// Returns (ra', dec') in degrees at `epoch_to`.
    pub fn propagate(
        ra: f64,
        dec: f64,
        pmra_mas_yr: f64,
        pmdec_mas_yr: f64,
        epoch_from: f64,
        epoch_to: f64,
    ) -> (f64, f64) {
        let dt = epoch_to - epoch_from;
        if dt.abs() < 1e-10 {
            return (ra, dec);
        }

        // pmra already includes cos(dec) factor per IAU convention,
        // so we divide by cos(dec) to get the actual RA shift.
        let cos_dec = dec.to_radians().cos();
        let ra_shift_deg = if cos_dec.abs() > 1e-10 {
            pmra_mas_yr * dt / (cos_dec * 3_600_000.0)
        } else {
            0.0
        };
        let dec_shift_deg = pmdec_mas_yr * dt / 3_600_000.0;

        let mut new_ra = ra + ra_shift_deg;
        let new_dec = (dec + dec_shift_deg).clamp(-90.0, 90.0);

        if new_ra < 0.0 {
            new_ra += 360.0;
        }
        if new_ra >= 360.0 {
            new_ra -= 360.0;
        }

        (new_ra, new_dec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_proper_motion() {
        let (ra, dec) = ProperMotionCorrector::propagate(180.0, 45.0, 0.0, 0.0, 2000.0, 2025.0);
        assert!((ra - 180.0).abs() < 1e-15);
        assert!((dec - 45.0).abs() < 1e-15);
    }

    #[test]
    fn zero_time_delta() {
        let (ra, dec) =
            ProperMotionCorrector::propagate(180.0, 45.0, 100.0, -50.0, 2000.0, 2000.0);
        assert!((ra - 180.0).abs() < 1e-15);
        assert!((dec - 45.0).abs() < 1e-15);
    }

    #[test]
    fn barnards_star() {
        // Barnard's star: RA ~269.452°, Dec ~4.694° (J2000)
        // pmra = -798.58 mas/yr, pmdec = 10328.12 mas/yr
        // Over 25 years (J2000 → J2025):
        //   Dec shift = 10328.12 * 25 / 3_600_000 = 0.07172° ≈ 258.2"
        let (ra, dec) = ProperMotionCorrector::propagate(
            269.452,
            4.694,
            -798.58,
            10328.12,
            2000.0,
            2025.0,
        );

        let dec_shift_arcsec = (dec - 4.694) * 3600.0;
        assert!(
            (dec_shift_arcsec - 258.2).abs() < 1.0,
            "Barnard's star Dec shift: expected ~258.2\", got {dec_shift_arcsec}\""
        );

        // RA shift: -798.58 * 25 / (cos(4.694°) * 3_600_000) ≈ -0.00557°
        let ra_shift = ra - 269.452;
        assert!(
            ra_shift < 0.0,
            "Barnard's star RA should decrease, got shift {ra_shift}"
        );
        assert!(
            ra_shift.abs() < 0.01,
            "Barnard's star RA shift should be small, got {ra_shift}"
        );
    }

    #[test]
    fn ra_wraparound() {
        // Star near RA=0 moving westward with large proper motion
        let (ra, _dec) =
            ProperMotionCorrector::propagate(0.001, 0.0, -10000.0, 0.0, 2000.0, 2025.0);
        // shift = -10000 * 25 / 3600000 = -0.0694° → 0.001 - 0.0694 < 0 → wraps
        assert!(ra > 359.0, "RA should wrap around, got {ra}");
    }
}
