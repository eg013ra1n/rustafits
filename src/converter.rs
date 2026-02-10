use std::path::Path;

use anyhow::{Context, Result};

use crate::output;
use crate::pipeline;
use crate::types::{ProcessConfig, ProcessedImage};

pub struct FitsConverter {
    downscale: usize,
    quality: u8,
    apply_debayer: bool,
    preview_mode: bool,
}

impl FitsConverter {
    pub fn new() -> Self {
        FitsConverter {
            downscale: 1,
            quality: 95,
            apply_debayer: true,
            preview_mode: false,
        }
    }

    pub fn with_downscale(mut self, factor: usize) -> Self {
        self.downscale = factor;
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

    /// Process a FITS/XISF image and return raw pixel data without writing to disk.
    ///
    /// Returns a `ProcessedImage` containing interleaved RGB u8 bytes,
    /// suitable for display in a GUI, web backend, or further processing.
    pub fn process<P: AsRef<Path>>(&self, input_path: P) -> Result<ProcessedImage> {
        let config = ProcessConfig {
            downscale_factor: self.downscale,
            jpeg_quality: self.quality,
            apply_debayer: self.apply_debayer,
            preview_mode: self.preview_mode,
            auto_stretch: true,
        };

        pipeline::process_image(input_path.as_ref(), &config)
            .context("Image processing failed")
    }

    /// Process a FITS/XISF image and save the result as JPEG or PNG.
    pub fn convert<P: AsRef<Path>, Q: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: Q,
    ) -> Result<()> {
        let image = self.process(&input_path)?;

        output::save_image(&image, output_path.as_ref(), self.quality)
            .context("Image save failed")?;

        Ok(())
    }
}

impl Default for FitsConverter {
    fn default() -> Self {
        Self::new()
    }
}
