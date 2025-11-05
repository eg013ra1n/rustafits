use anyhow::Result;
use ndarray::{Array2, Array3};
use crate::fits::{FitsImage, FitsImageData, RowOrder};
use crate::downscale::{downscale_mono, downscale_color};
use crate::debayer::super_pixel_debayer;
use crate::stretch::{stretch_mono, stretch_color, StretchParams};

pub struct ProcessorConfig {
    /// Downscale factor (1 = no downscaling)
    pub downscale_factor: usize,

    /// Manual stretch parameters (None = auto-stretch)
    pub manual_stretch: Option<Vec<StretchParams>>,

    /// Whether to apply debayering if Bayer pattern is detected
    pub apply_debayer: bool,

    /// Skip stretching and return raw 16-bit data (for TIFF output)
    pub skip_stretch: bool,
}

impl Default for ProcessorConfig {
    fn default() -> Self {
        ProcessorConfig {
            downscale_factor: 1,
            manual_stretch: None,
            apply_debayer: true,
            skip_stretch: false,
        }
    }
}

pub enum ProcessedImage {
    Mono(Array2<u8>),
    RGB(Array3<u8>),
    RGB16(Array3<f32>),  // Raw 16-bit float data (before stretching)
}

impl ProcessedImage {
    pub fn width(&self) -> usize {
        match self {
            ProcessedImage::Mono(arr) => arr.dim().1,
            ProcessedImage::RGB(arr) => arr.dim().2,
            ProcessedImage::RGB16(arr) => arr.dim().2,
        }
    }

    pub fn height(&self) -> usize {
        match self {
            ProcessedImage::Mono(arr) => arr.dim().0,
            ProcessedImage::RGB(arr) => arr.dim().1,
            ProcessedImage::RGB16(arr) => arr.dim().1,
        }
    }

    pub fn is_color(&self) -> bool {
        matches!(self, ProcessedImage::RGB(_) | ProcessedImage::RGB16(_))
    }

    /// Convert to interleaved RGB bytes (for JPEG encoding)
    pub fn to_rgb_bytes(&self) -> Vec<u8> {
        match self {
            ProcessedImage::Mono(arr) => {
                // Convert grayscale to RGB by duplicating the value
                let (height, width) = arr.dim();
                let mut rgb = Vec::with_capacity(height * width * 3);
                for &pixel in arr.iter() {
                    rgb.push(pixel);
                    rgb.push(pixel);
                    rgb.push(pixel);
                }
                rgb
            }
            ProcessedImage::RGB(arr) => {
                // Convert from planar (CHW) to interleaved (HWC)
                let (_, height, width) = arr.dim();
                let mut rgb = Vec::with_capacity(height * width * 3);

                for y in 0..height {
                    for x in 0..width {
                        rgb.push(arr[[0, y, x]]);  // R
                        rgb.push(arr[[1, y, x]]);  // G
                        rgb.push(arr[[2, y, x]]);  // B
                    }
                }
                rgb
            }
            ProcessedImage::RGB16(_) => {
                panic!("Cannot convert RGB16 to 8-bit bytes. Use save_tiff_16bit instead.");
            }
        }
    }

    /// Flip image vertically (for ROWORDER handling)
    pub fn flip_vertical(&mut self) {
        match self {
            ProcessedImage::Mono(arr) => {
                arr.invert_axis(ndarray::Axis(0));
            }
            ProcessedImage::RGB(arr) => {
                arr.invert_axis(ndarray::Axis(1));
            }
            ProcessedImage::RGB16(arr) => {
                arr.invert_axis(ndarray::Axis(1));
            }
        }
    }
}

pub fn process_fits(fits: FitsImage, config: ProcessorConfig) -> Result<ProcessedImage> {
    let row_order = fits.metadata.row_order;
    let has_bayer = fits.has_bayer_pattern();
    let bayer_pattern = fits.metadata.bayer_pattern;

    let mut result = match fits.data {
        FitsImageData::U16(data) => {
            process_mono_u16(data, has_bayer, bayer_pattern, &config)?
        }
        FitsImageData::F32(data) => {
            process_mono_f32(data, has_bayer, bayer_pattern, &config)?
        }
        FitsImageData::U16Color(data) => {
            process_color_u16(data, &config)?
        }
        FitsImageData::F32Color(data) => {
            process_color_f32(data, &config)?
        }
    };

    // Handle row order (flip if needed)
    if matches!(row_order, RowOrder::TopDown) {
        result.flip_vertical();
    }

    Ok(result)
}

fn process_mono_u16(
    mut data: Array2<u16>,
    has_bayer: bool,
    bayer_pattern: Option<crate::fits::BayerPattern>,
    config: &ProcessorConfig,
) -> Result<ProcessedImage> {
    // Step 1: Downscale if needed
    if config.downscale_factor > 1 {
        data = downscale_mono(&data, config.downscale_factor);
    }

    // Step 2: Check for Bayer pattern and debayer
    if config.apply_debayer && has_bayer {
        let pattern = bayer_pattern.unwrap();
        let rgb_f32 = super_pixel_debayer(&data, pattern);

        // Step 3: Stretch or return raw
        if config.skip_stretch {
            Ok(ProcessedImage::RGB16(rgb_f32))
        } else {
            let stretched = stretch_color(&rgb_f32, config.manual_stretch.clone());
            Ok(ProcessedImage::RGB(stretched))
        }
    } else {
        // Step 3: Stretch mono image (no RGB16 output for mono)
        let stretched = stretch_mono(&data, config.manual_stretch.as_ref().map(|v| v[0]));
        Ok(ProcessedImage::Mono(stretched))
    }
}

fn process_mono_f32(
    mut data: Array2<f32>,
    has_bayer: bool,
    bayer_pattern: Option<crate::fits::BayerPattern>,
    config: &ProcessorConfig,
) -> Result<ProcessedImage> {
    // Step 1: Downscale if needed
    if config.downscale_factor > 1 {
        data = downscale_mono(&data, config.downscale_factor);
    }

    // Step 2: Check for Bayer pattern and debayer
    if config.apply_debayer && has_bayer {
        let pattern = bayer_pattern.unwrap();
        let rgb_f32 = super_pixel_debayer(&data, pattern);

        // Step 3: Stretch or return raw
        if config.skip_stretch {
            Ok(ProcessedImage::RGB16(rgb_f32))
        } else {
            let stretched = stretch_color(&rgb_f32, config.manual_stretch.clone());
            Ok(ProcessedImage::RGB(stretched))
        }
    } else {
        // Step 3: Stretch mono image (no RGB16 output for mono)
        let stretched = stretch_mono(&data, config.manual_stretch.as_ref().map(|v| v[0]));
        Ok(ProcessedImage::Mono(stretched))
    }
}

fn process_color_u16(
    mut data: Array3<u16>,
    config: &ProcessorConfig,
) -> Result<ProcessedImage> {
    // Step 1: Downscale if needed
    if config.downscale_factor > 1 {
        data = downscale_color(&data, config.downscale_factor);
    }

    // Step 2: Convert to f32 for stretching
    let data_f32 = data.mapv(|v| v as f32);

    // Step 3: Stretch
    let stretched = stretch_color(&data_f32, config.manual_stretch.clone());

    Ok(ProcessedImage::RGB(stretched))
}

fn process_color_f32(
    mut data: Array3<f32>,
    config: &ProcessorConfig,
) -> Result<ProcessedImage> {
    // Step 1: Downscale if needed
    if config.downscale_factor > 1 {
        data = downscale_color(&data, config.downscale_factor);
    }

    // Step 2: Stretch
    let stretched = stretch_color(&data, config.manual_stretch.clone());

    Ok(ProcessedImage::RGB(stretched))
}
