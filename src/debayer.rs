use ndarray::{Array2, Array3};
use crate::fits::BayerPattern;
use std::os::raw::c_int;

// C debayer function bindings
extern "C" {
    fn super_pixel_rggb_u16(buf: *const u16, newbuf: *mut f32, width: c_int, height: c_int);
    fn super_pixel_rggb_f32(buf: *const f32, newbuf: *mut f32, width: c_int, height: c_int);
    fn super_pixel_bggr_u16(buf: *const u16, newbuf: *mut f32, width: c_int, height: c_int);
    fn super_pixel_bggr_f32(buf: *const f32, newbuf: *mut f32, width: c_int, height: c_int);
    fn super_pixel_gbrg_u16(buf: *const u16, newbuf: *mut f32, width: c_int, height: c_int);
    fn super_pixel_gbrg_f32(buf: *const f32, newbuf: *mut f32, width: c_int, height: c_int);
    fn super_pixel_grbg_u16(buf: *const u16, newbuf: *mut f32, width: c_int, height: c_int);
    fn super_pixel_grbg_f32(buf: *const f32, newbuf: *mut f32, width: c_int, height: c_int);
}

/// Super-pixel debayering algorithm using C implementation
/// Reduces resolution by 50% in each dimension but avoids interpolation artifacts
/// Each 2x2 Bayer block becomes 1 RGB pixel
pub fn super_pixel_debayer<T>(data: &Array2<T>, pattern: BayerPattern) -> Array3<f32>
where
    T: Copy + Into<f32>,
{
    let (height, width) = data.dim();

    // Ensure even dimensions (crop if needed)
    let new_height = (height / 2) * 2;
    let new_width = (width / 2) * 2;

    let out_height = new_height / 2;
    let out_width = new_width / 2;
    let out_plane_size = out_height * out_width;

    // Allocate output buffer (3 planes for RGB)
    let mut output = vec![0.0f32; out_plane_size * 3];

    // Check if input is u16 or f32 and call appropriate C function
    unsafe {
        // Try to get u16 pointer first
        if let Some(u16_data) = try_as_u16_slice(data) {
            match pattern {
                BayerPattern::RGGB => {
                    super_pixel_rggb_u16(u16_data.as_ptr(), output.as_mut_ptr(), new_width as c_int, new_height as c_int);
                }
                BayerPattern::BGGR => {
                    super_pixel_bggr_u16(u16_data.as_ptr(), output.as_mut_ptr(), new_width as c_int, new_height as c_int);
                }
                BayerPattern::GBRG => {
                    super_pixel_gbrg_u16(u16_data.as_ptr(), output.as_mut_ptr(), new_width as c_int, new_height as c_int);
                }
                BayerPattern::GRBG => {
                    super_pixel_grbg_u16(u16_data.as_ptr(), output.as_mut_ptr(), new_width as c_int, new_height as c_int);
                }
            }
        } else {
            // Convert to f32 if not u16
            let f32_data: Vec<f32> = data.iter().map(|&x| x.into()).collect();
            match pattern {
                BayerPattern::RGGB => {
                    super_pixel_rggb_f32(f32_data.as_ptr(), output.as_mut_ptr(), new_width as c_int, new_height as c_int);
                }
                BayerPattern::BGGR => {
                    super_pixel_bggr_f32(f32_data.as_ptr(), output.as_mut_ptr(), new_width as c_int, new_height as c_int);
                }
                BayerPattern::GBRG => {
                    super_pixel_gbrg_f32(f32_data.as_ptr(), output.as_mut_ptr(), new_width as c_int, new_height as c_int);
                }
                BayerPattern::GRBG => {
                    super_pixel_grbg_f32(f32_data.as_ptr(), output.as_mut_ptr(), new_width as c_int, new_height as c_int);
                }
            }
        }
    }

    // Convert output to Array3 in planar format (C, H, W)
    let rgb = Array3::from_shape_vec((3, out_height, out_width), output)
        .expect("Failed to create output array");

    rgb
}

// Helper function to try casting Array2<T> to &[u16]
fn try_as_u16_slice<T: Copy>(data: &Array2<T>) -> Option<Vec<u16>> {
    // This is a bit of a hack - check if T is u16 by trying to cast
    if std::mem::size_of::<T>() == std::mem::size_of::<u16>() {
        // Safety: We've checked the size matches
        let slice = unsafe {
            std::slice::from_raw_parts(
                data.as_ptr() as *const u16,
                data.len()
            )
        };
        Some(slice.to_vec())
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_super_pixel_rggb() {
        // Create a simple 4x4 Bayer pattern
        // R G R G
        // G B G B
        // R G R G
        // G B G B
        let data = Array2::from_shape_fn((4, 4), |(y, x)| {
            let is_red_row = y % 2 == 0;
            let is_red_col = x % 2 == 0;

            if is_red_row && is_red_col {
                100u16  // R
            } else if !is_red_row && !is_red_col {
                50u16   // B
            } else {
                75u16   // G
            }
        });

        let rgb = super_pixel_debayer(&data, BayerPattern::RGGB);

        // Should be 2x2 output
        assert_eq!(rgb.dim(), (3, 2, 2));

        // Check first pixel
        assert_eq!(rgb[[0, 0, 0]], 100.0);  // R
        assert_eq!(rgb[[1, 0, 0]], 75.0);   // G (average of two 75s)
        assert_eq!(rgb[[2, 0, 0]], 50.0);   // B
    }

    #[test]
    fn test_super_pixel_bggr() {
        let data = Array2::from_shape_fn((4, 4), |(y, x)| {
            let is_red_row = y % 2 == 1;
            let is_red_col = x % 2 == 1;

            if is_red_row && is_red_col {
                100u16  // R
            } else if !is_red_row && !is_red_col {
                50u16   // B
            } else {
                75u16   // G
            }
        });

        let rgb = super_pixel_debayer(&data, BayerPattern::BGGR);

        // Check first pixel
        assert_eq!(rgb[[0, 0, 0]], 100.0);  // R
        assert_eq!(rgb[[1, 0, 0]], 75.0);   // G
        assert_eq!(rgb[[2, 0, 0]], 50.0);   // B
    }
}
