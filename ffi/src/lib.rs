use fits_converter::FitsConverter;
use libc::{c_char, c_int, size_t};
use std::ffi::CStr;
use std::ptr;

// C debayer function bindings
extern "C" {
    pub fn super_pixel_rggb_u16(buf: *const u16, newbuf: *mut f32, width: c_int, height: c_int);
    pub fn super_pixel_rggb_f32(buf: *const f32, newbuf: *mut f32, width: c_int, height: c_int);
    pub fn super_pixel_bggr_u16(buf: *const u16, newbuf: *mut f32, width: c_int, height: c_int);
    pub fn super_pixel_bggr_f32(buf: *const f32, newbuf: *mut f32, width: c_int, height: c_int);
    pub fn super_pixel_gbrg_u16(buf: *const u16, newbuf: *mut f32, width: c_int, height: c_int);
    pub fn super_pixel_gbrg_f32(buf: *const f32, newbuf: *mut f32, width: c_int, height: c_int);
    pub fn super_pixel_grbg_u16(buf: *const u16, newbuf: *mut f32, width: c_int, height: c_int);
    pub fn super_pixel_grbg_f32(buf: *const f32, newbuf: *mut f32, width: c_int, height: c_int);
}

/// Opaque handle to FitsConverter
pub struct FitsConverterHandle {
    converter: FitsConverter,
}

/// Error codes
pub const FITS_OK: c_int = 0;
pub const FITS_ERROR_NULL_POINTER: c_int = -1;
pub const FITS_ERROR_INVALID_UTF8: c_int = -2;
pub const FITS_ERROR_CONVERSION: c_int = -3;

/// Create a new FITS converter with default settings
/// Returns: Opaque pointer to FitsConverter (NULL on error)
#[no_mangle]
pub extern "C" fn fits_converter_create() -> *mut FitsConverterHandle {
    let converter = FitsConverter::new();
    let handle = Box::new(FitsConverterHandle { converter });
    Box::into_raw(handle)
}

/// Set downscale factor
/// handle: Converter handle
/// factor: Downscale factor (1 = no downscaling)
/// Returns: FITS_OK on success
#[no_mangle]
pub extern "C" fn fits_converter_set_downscale(
    handle: *mut FitsConverterHandle,
    factor: size_t,
) -> c_int {
    if handle.is_null() {
        return FITS_ERROR_NULL_POINTER;
    }

    unsafe {
        let handle_ref = &mut *handle;
        handle_ref.converter = std::mem::take(&mut handle_ref.converter)
            .with_downscale(factor);
    }

    FITS_OK
}

/// Set JPEG quality
/// handle: Converter handle
/// quality: JPEG quality (1-100)
/// Returns: FITS_OK on success
#[no_mangle]
pub extern "C" fn fits_converter_set_quality(
    handle: *mut FitsConverterHandle,
    quality: c_int,
) -> c_int {
    if handle.is_null() {
        return FITS_ERROR_NULL_POINTER;
    }

    unsafe {
        let handle_ref = &mut *handle;
        handle_ref.converter = std::mem::take(&mut handle_ref.converter)
            .with_quality(quality as u8);
    }

    FITS_OK
}

/// Disable automatic debayering
/// handle: Converter handle
/// Returns: FITS_OK on success
#[no_mangle]
pub extern "C" fn fits_converter_disable_debayer(
    handle: *mut FitsConverterHandle,
) -> c_int {
    if handle.is_null() {
        return FITS_ERROR_NULL_POINTER;
    }

    unsafe {
        let handle_ref = &mut *handle;
        handle_ref.converter = std::mem::take(&mut handle_ref.converter)
            .without_debayer();
    }

    FITS_OK
}

/// Convert FITS file to JPEG
/// handle: Converter handle
/// input_path: Path to input FITS file (null-terminated C string)
/// output_path: Path to output JPEG file (null-terminated C string)
/// Returns: FITS_OK on success, error code on failure
#[no_mangle]
pub extern "C" fn fits_converter_convert(
    handle: *mut FitsConverterHandle,
    input_path: *const c_char,
    output_path: *const c_char,
) -> c_int {
    if handle.is_null() || input_path.is_null() || output_path.is_null() {
        return FITS_ERROR_NULL_POINTER;
    }

    unsafe {
        let handle_ref = &*handle;

        // Convert C strings to Rust strings
        let input_cstr = match CStr::from_ptr(input_path).to_str() {
            Ok(s) => s,
            Err(_) => return FITS_ERROR_INVALID_UTF8,
        };

        let output_cstr = match CStr::from_ptr(output_path).to_str() {
            Ok(s) => s,
            Err(_) => return FITS_ERROR_INVALID_UTF8,
        };

        // Perform conversion
        match handle_ref.converter.convert(input_cstr, output_cstr) {
            Ok(_) => FITS_OK,
            Err(_) => FITS_ERROR_CONVERSION,
        }
    }
}

/// Get last error message (not implemented in this simple version)
/// Returns: Pointer to error message string (may be NULL)
#[no_mangle]
pub extern "C" fn fits_converter_get_error() -> *const c_char {
    ptr::null()
}

/// Destroy the FITS converter and free memory
/// handle: Converter handle to destroy
#[no_mangle]
pub extern "C" fn fits_converter_destroy(handle: *mut FitsConverterHandle) {
    if !handle.is_null() {
        unsafe {
            let _ = Box::from_raw(handle);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_destroy() {
        let handle = fits_converter_create();
        assert!(!handle.is_null());
        fits_converter_destroy(handle);
    }
}
