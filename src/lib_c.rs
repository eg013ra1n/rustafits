use std::ffi::CString;
use std::os::raw::{c_char, c_int};
use std::path::Path;
use anyhow::{Context, Result};

// C structures and functions (from fits_processor.h)
#[repr(C)]
struct ProcessConfig {
    downscale_factor: c_int,
    jpeg_quality: c_int,
    apply_debayer: c_int,
    preview_mode: c_int,
    auto_stretch: c_int,
    manual_shadows: f32,
    manual_highlights: f32,
    manual_midtones: f32,
}

#[repr(C)]
struct ProcessedImage {
    data: *mut u8,
    width: usize,
    height: usize,
    is_color: c_int,
}

extern "C" {
    fn process_fits_file(
        fits_path: *const c_char,
        config: *const ProcessConfig,
        out_image: *mut ProcessedImage,
    ) -> c_int;

    fn save_jpeg(
        image: *const ProcessedImage,
        output_path: *const c_char,
        quality: c_int,
    ) -> c_int;

    fn free_processed_image(image: *mut ProcessedImage);

    fn get_last_error() -> *const c_char;
}

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
        let input_cstr = CString::new(input_path.as_ref().to_str().unwrap())
            .context("Invalid input path")?;
        let output_cstr = CString::new(output_path.as_ref().to_str().unwrap())
            .context("Invalid output path")?;

        let config = ProcessConfig {
            downscale_factor: self.downscale as c_int,
            jpeg_quality: self.quality as c_int,
            apply_debayer: if self.apply_debayer { 1 } else { 0 },
            preview_mode: if self.preview_mode { 1 } else { 0 },
            auto_stretch: 1,
            manual_shadows: 0.0,
            manual_highlights: 1.0,
            manual_midtones: 0.5,
        };

        let mut image = ProcessedImage {
            data: std::ptr::null_mut(),
            width: 0,
            height: 0,
            is_color: 0,
        };

        unsafe {
            let result = process_fits_file(
                input_cstr.as_ptr(),
                &config as *const ProcessConfig,
                &mut image as *mut ProcessedImage,
            );

            if result != 0 {
                let err_ptr = get_last_error();
                let err_msg = if !err_ptr.is_null() {
                    std::ffi::CStr::from_ptr(err_ptr)
                        .to_string_lossy()
                        .into_owned()
                } else {
                    "Unknown error".to_string()
                };
                anyhow::bail!("FITS processing failed: {}", err_msg);
            }

            let result = save_jpeg(
                &image as *const ProcessedImage,
                output_cstr.as_ptr(),
                self.quality as c_int,
            );

            free_processed_image(&mut image as *mut ProcessedImage);

            if result != 0 {
                anyhow::bail!("JPEG save failed");
            }
        }

        Ok(())
    }

}

impl Default for FitsConverter {
    fn default() -> Self {
        Self::new()
    }
}
