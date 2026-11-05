#[cfg(feature = "debug-pipeline")]
pub mod analysis;
#[cfg(not(feature = "debug-pipeline"))]
mod analysis;
mod annotate;
mod converter;
mod geometry;
#[cfg(feature = "debug-pipeline")]
pub mod formats;
#[cfg(not(feature = "debug-pipeline"))]
mod formats;
mod output;
mod pipeline;
pub mod platesolving;
pub mod processing;
pub mod star_measure;
pub mod stats;
#[cfg(test)]
mod test_support;
mod types;

pub use analysis::{
    AnalysisConfig, AnalysisResult, DetectionLevels, FastAnalysisResult, FastDetectTiming,
    FastStar, FitMethod, ImageAnalyzer, StageTiming, StarMetrics,
};
pub use annotate::{
    annotate_image, annotate_image_with, compute_annotations, compute_overlay,
    create_annotation_layer, create_overlay_layer, AnnotationConfig, ColorScheme, OverlayStar,
    StarAnnotation, StarOverlay,
};
pub use converter::ImageConverter;
/// In-memory JPEG encoding (baseline, 4:2:0). Accepts RGB (3 ch) or RGBA (4 ch)
/// — alpha is discarded by the encoder, so callers must not pre-strip it.
pub use output::encode_jpeg;
pub use output::{encode_png, OutputFormat, PngCompression};
pub use rayon::{ThreadPool, ThreadPoolBuilder};
/// The Moffat LM solve counters (feature `lm-counters`; zero numeric effect), reachable without `debug-pipeline`.
#[cfg(feature = "lm-counters")]
pub use analysis::fitting::lm_counters;
pub use types::{
    BayerPattern, DataType, FileFormat, HeaderKeyword, ImageHeader, ImageMetadata, PixelData,
    ProcessedImage,
};
