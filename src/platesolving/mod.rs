pub mod pattern_matcher;
pub mod projection;
pub mod proper_motion;
pub mod ransac;
pub mod transform;
pub mod types;
pub mod wcs;

pub use pattern_matcher::{
    build_quads, fit_affine_from_centers, match_quads, AffineTransform, Quad, QuadMatch,
};
pub use projection::GnomonicProjection;
pub use proper_motion::ProperMotionCorrector;
pub use ransac::RansacFilter;
pub use transform::{PixelTransform, TransformFitter};
pub use types::{
    CatalogStar, FitModel, ImageStar, PatternMatcherConfig, ProjectedStar, RansacConfig,
    SipCoefficients, SolveHints, StarMatch,
};
pub use wcs::WcsSolution;
