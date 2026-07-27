#[cfg(feature = "debug-pipeline")]
pub mod analysis;
#[cfg(not(feature = "debug-pipeline"))]
mod analysis;
mod annotate;
mod converter;
#[cfg(feature = "debug-pipeline")]
pub mod formats;
#[cfg(not(feature = "debug-pipeline"))]
mod formats;
mod output;
mod pipeline;
pub mod platesolving;
mod processing;
mod types;

pub use analysis::{
    AnalysisConfig, AnalysisResult, FastAnalysisResult, FastDetectTiming, FastStar, FitMethod,
    ImageAnalyzer, StageTiming, StarMetrics,
};
pub use annotate::{
    annotate_image, compute_annotations, create_annotation_layer, AnnotationConfig, ColorScheme,
    StarAnnotation,
};
pub use converter::ImageConverter;
/// In-memory JPEG encoding (baseline, 4:2:0). Accepts RGB (3 ch) or RGBA (4 ch)
/// — alpha is discarded by the encoder, so callers must not pre-strip it.
pub use output::encode_jpeg;
pub use rayon::{ThreadPool, ThreadPoolBuilder};
pub use types::{BayerPattern, DataType, ImageMetadata, PixelData, ProcessedImage};
