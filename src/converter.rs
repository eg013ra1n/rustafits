use std::path::Path;

use anyhow::{Context, Result};

use crate::output;
use crate::pipeline;
use crate::types::ProcessConfig;

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

    pub fn convert<P: AsRef<Path>, Q: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: Q,
    ) -> Result<()> {
        let input = input_path.as_ref();
        let output = output_path.as_ref();

        let config = ProcessConfig {
            downscale_factor: self.downscale,
            jpeg_quality: self.quality,
            apply_debayer: self.apply_debayer,
            preview_mode: self.preview_mode,
            auto_stretch: true,
        };

        let image = pipeline::process_image(input, &config)
            .context("Image processing failed")?;

        output::save_image(&image, output, self.quality)
            .context("Image save failed")?;

        Ok(())
    }
}

impl Default for FitsConverter {
    fn default() -> Self {
        Self::new()
    }
}
