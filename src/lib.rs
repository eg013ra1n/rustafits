mod analysis;
mod annotate;
mod converter;
mod formats;
mod output;
mod pipeline;
mod processing;
mod types;

pub use analysis::{AnalysisConfig, AnalysisResult, ImageAnalyzer, StarMetrics};
pub use annotate::{
    annotate_image, compute_annotations, create_annotation_layer, AnnotationConfig, ColorScheme,
    StarAnnotation,
};
pub use converter::ImageConverter;
pub use rayon::{ThreadPool, ThreadPoolBuilder};
pub use types::ProcessedImage;
