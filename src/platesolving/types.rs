/// A star from a reference catalog, with epoch-corrected coordinates.
#[derive(Clone, Debug)]
pub struct CatalogStar {
    /// Right ascension in degrees (epoch-corrected).
    pub ra: f64,
    /// Declination in degrees (epoch-corrected).
    pub dec: f64,
    /// Apparent magnitude.
    pub mag: f32,
}

/// A catalog star projected onto a gnomonic tangent plane.
#[derive(Clone, Debug)]
pub struct ProjectedStar {
    /// Tangent plane coordinate (radians).
    pub xi: f64,
    /// Tangent plane coordinate (radians).
    pub eta: f64,
    /// Apparent magnitude.
    pub mag: f32,
    /// Original RA in degrees.
    pub ra: f64,
    /// Original Dec in degrees.
    pub dec: f64,
}

/// A matched pair: image star index to catalog star index.
#[derive(Clone, Debug)]
pub struct StarMatch {
    pub image_idx: usize,
    pub catalog_idx: usize,
    pub residual_px: f64,
}

/// Hints for the plate solver derived from FITS headers or user input.
#[derive(Clone, Debug, Default)]
pub struct SolveHints {
    /// Approximate RA of field center in degrees.
    pub ra: Option<f64>,
    /// Approximate Dec of field center in degrees.
    pub dec: Option<f64>,
    /// Approximate field of view in degrees.
    pub fov_deg: Option<f64>,
    /// Approximate rotation angle in degrees.
    pub rotation: Option<f64>,
    /// Expected pixel scale in arcseconds per pixel.
    pub pixel_scale_arcsec: Option<f64>,
}

/// Configuration for the pattern matcher.
#[derive(Clone, Debug)]
pub struct PatternMatcherConfig {
    /// Maximum number of brightest stars to use.
    pub max_stars: usize,
    /// Quantization bin size for hash keys.
    pub hash_tolerance: f64,
    /// Expected arcsec/px — reject quads whose scale differs >20%.
    pub scale_hint: Option<f64>,
    /// Check adjacent hash bins (±1 in each dimension).
    pub multi_probe: bool,
}

impl Default for PatternMatcherConfig {
    fn default() -> Self {
        Self {
            max_stars: 100,
            hash_tolerance: 0.01,
            scale_hint: None,
            multi_probe: true,
        }
    }
}

/// Configuration for RANSAC outlier rejection.
#[derive(Clone, Debug)]
pub struct RansacConfig {
    /// Inlier distance threshold in pixels.
    pub threshold_px: f64,
    /// Maximum iterations.
    pub max_iterations: u32,
    /// Minimum inlier count to accept.
    pub min_inliers: usize,
}

impl Default for RansacConfig {
    fn default() -> Self {
        Self {
            threshold_px: 2.5,
            max_iterations: 100,
            min_inliers: 6,
        }
    }
}

/// Transform model to fit.
#[derive(Clone, Debug)]
pub enum FitModel {
    /// 6-parameter affine transform.
    Affine,
    /// 8-parameter projective (homography) transform.
    Projective,
    /// SIP polynomial distortion correction.
    Sip { order: u8 },
}

/// SIP polynomial coefficients (max order 5).
#[derive(Clone, Debug)]
pub struct SipCoefficients {
    pub order: u8,
    pub coeffs: [[f64; 6]; 6],
}

impl SipCoefficients {
    pub fn new(order: u8) -> Self {
        assert!(order <= 5, "SIP order must be <= 5");
        Self {
            order,
            coeffs: [[0.0; 6]; 6],
        }
    }

    pub fn evaluate(&self, u: f64, v: f64) -> f64 {
        let mut result = 0.0;
        for i in 0..=(self.order as usize) {
            for j in 0..=(self.order as usize) {
                if i + j >= 2 && i + j <= self.order as usize {
                    result += self.coeffs[i][j] * u.powi(i as i32) * v.powi(j as i32);
                }
            }
        }
        result
    }
}

/// A detected star in image pixel coordinates, used as input to the pattern matcher.
/// This mirrors the fields needed from the analysis::DetectedStar struct
/// but uses f64 for numerical precision in plate solving math.
#[derive(Clone, Debug)]
pub struct ImageStar {
    /// Centroid X in pixels (subpixel).
    pub x: f64,
    /// Centroid Y in pixels (subpixel).
    pub y: f64,
    /// Total flux (for brightness sorting).
    pub flux: f64,
}
