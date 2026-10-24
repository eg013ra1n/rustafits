use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};

use crate::formats;
use crate::output::{self, PngCompression};
use crate::pipeline;
use crate::types::{ImageHeader, ImageMetadata, PixelData, ProcessConfig, ProcessedImage};

pub struct ImageConverter {
    downscale: usize,
    quality: u8,
    apply_debayer: bool,
    preview_mode: bool,
    rgba_output: bool,
    vng_debayer: bool,
    png_compression: PngCompression,
    thread_pool: Option<Arc<rayon::ThreadPool>>,
}

impl ImageConverter {
    pub fn new() -> Self {
        ImageConverter {
            downscale: 1,
            quality: 95,
            apply_debayer: true,
            preview_mode: false,
            rgba_output: false,
            vng_debayer: false,
            png_compression: PngCompression::Fast,
            thread_pool: None,
        }
    }

    pub fn with_downscale(mut self, factor: usize) -> Self {
        // Factor 0 would divide-by-zero in the downscale kernels.
        self.downscale = factor.max(1);
        self
    }

    pub fn with_quality(mut self, quality: u8) -> Self {
        self.quality = quality.clamp(1, 100);
        self
    }

    pub fn without_debayer(mut self) -> Self {
        self.apply_debayer = false;
        self
    }

    pub fn with_preview_mode(mut self) -> Self {
        self.preview_mode = true;
        self
    }

    /// Output RGBA (4 bytes/pixel) instead of RGB, suitable for HTML Canvas `ImageData`.
    pub fn with_rgba_output(mut self) -> Self {
        self.rgba_output = true;
        self
    }

    /// Debayer CFA input at native resolution with the gradient method
    /// instead of folding 2x2 tiles into single pixels. See
    /// [`ProcessConfig::vng_debayer`].
    pub fn with_vng_debayer(mut self) -> Self {
        self.vng_debayer = true;
        self
    }

    /// PNG compression level used by [`convert`](Self::convert) when the
    /// output is `.png` (default [`PngCompression::Fast`]; ignored for JPEG).
    pub fn with_png_compression(mut self, level: PngCompression) -> Self {
        self.png_compression = level;
        self
    }

    pub fn with_thread_pool(mut self, pool: Arc<rayon::ThreadPool>) -> Self {
        self.thread_pool = Some(pool);
        self
    }

    /// Read raw pixel data from a FITS/XISF file without any processing.
    pub fn read_raw<P: AsRef<Path>>(path: P) -> Result<(ImageMetadata, PixelData)> {
        formats::read_image(path.as_ref()).context("Failed to read image")
    }

    /// Read a FITS/XISF header without its pixel data: dimensions, sample
    /// type, Bayer pattern, orientation and every keyword. Fast — a few KB of
    /// I/O — so suitable for listing large sets of files.
    pub fn read_header<P: AsRef<Path>>(path: P) -> Result<ImageHeader> {
        formats::read_header(path.as_ref()).context("Failed to read header")
    }

    /// Process pre-read image data (skips file I/O).
    pub fn process_data(
        &self,
        meta: ImageMetadata,
        pixels: PixelData,
    ) -> Result<ProcessedImage> {
        let config = self.build_config();
        match &self.thread_pool {
            Some(pool) => pool.install(|| pipeline::process_image_data(meta, pixels, &config)),
            None => pipeline::process_image_data(meta, pixels, &config),
        }
        .context("Image processing failed")
    }

    /// Process a FITS/XISF image and return raw pixel data without writing to disk.
    ///
    /// Returns a `ProcessedImage` containing interleaved RGB u8 bytes,
    /// suitable for display in a GUI, web backend, or further processing.
    pub fn process<P: AsRef<Path>>(&self, input_path: P) -> Result<ProcessedImage> {
        let config = self.build_config();
        let path = input_path.as_ref();
        match &self.thread_pool {
            Some(pool) => pool.install(|| pipeline::process_image(path, &config)),
            None => pipeline::process_image(path, &config),
        }
        .context("Image processing failed")
    }

    fn build_config(&self) -> ProcessConfig {
        ProcessConfig {
            downscale_factor: self.downscale,
            jpeg_quality: self.quality,
            apply_debayer: self.apply_debayer,
            preview_mode: self.preview_mode,
            auto_stretch: true,
            rgba_output: self.rgba_output,
            vng_debayer: self.vng_debayer,
        }
    }

    /// Save a `ProcessedImage` to disk as JPEG or PNG.
    pub fn save_processed<P: AsRef<Path>>(
        image: &ProcessedImage,
        output_path: P,
        quality: u8,
    ) -> Result<()> {
        output::save_image(image, output_path.as_ref(), quality)
            .context("Image save failed")
    }

    /// [`save_processed`](Self::save_processed) with an explicit PNG
    /// compression level (ignored for JPEG).
    pub fn save_processed_with<P: AsRef<Path>>(
        image: &ProcessedImage,
        output_path: P,
        quality: u8,
        png_compression: PngCompression,
    ) -> Result<()> {
        output::save_image_with(image, output_path.as_ref(), quality, png_compression)
            .context("Image save failed")
    }

    /// Process a FITS/XISF image and save the result as JPEG or PNG.
    pub fn convert<P: AsRef<Path>, Q: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: Q,
    ) -> Result<()> {
        // Reject an unsupported output extension before the (expensive) processing.
        output::OutputFormat::from_path(output_path.as_ref())?;
        let image = self.process(&input_path)?;

        output::save_image_with(&image, output_path.as_ref(), self.quality, self.png_compression)
            .context("Image save failed")?;

        Ok(())
    }
}

impl Default for ImageConverter {
    fn default() -> Self {
        Self::new()
    }
}
