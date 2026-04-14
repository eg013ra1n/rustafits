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
pub use rayon::{ThreadPool, ThreadPoolBuilder};
pub use types::{BayerPattern, ImageMetadata, PixelData, ProcessedImage};
