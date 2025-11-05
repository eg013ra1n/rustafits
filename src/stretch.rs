use ndarray::{Array2, Array3, Axis};
use rayon::prelude::*;
use std::os::raw::c_int;

const MAX_SAMPLES: usize = 500_000;

// C stretch function bindings
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CStretchParams {
    pub shadows: f32,
    pub highlights: f32,
    pub midtones: f32,
    pub max_input: c_int,
}

extern "C" {
    fn compute_stretch_params(data: *const f32, width: c_int, height: c_int,
                             input_range: c_int, params: *mut CStretchParams);
    fn apply_stretch(data: *mut f32, width: c_int, height: c_int,
                    params: *const CStretchParams);
}

/// Parameters for stretching a single channel
#[derive(Debug, Clone, Copy)]
pub struct StretchParams {
    pub shadows: f32,
    pub highlights: f32,
    pub midtones: f32,
}

impl From<CStretchParams> for StretchParams {
    fn from(c_params: CStretchParams) -> Self {
        StretchParams {
            shadows: c_params.shadows,
            highlights: c_params.highlights,
            midtones: c_params.midtones,
        }
    }
}

impl StretchParams {
    /// Compute stretch parameters automatically using median-based statistics
    /// This is the QuickFits/PixInsight algorithm (using C implementation)
    pub fn compute_auto(data: &[f32], max_input: f32) -> Self {
        if data.is_empty() {
            return Self::default();
        }

        // Use C implementation for exact QuickFits algorithm
        let mut c_params = CStretchParams {
            shadows: 0.0,
            highlights: 1.0,
            midtones: 0.5,
            max_input: max_input as c_int,
        };

        // Assume data is a flat array representing width x height
        let width = data.len() as c_int;
        let height = 1;
        let input_range = max_input as c_int;

        unsafe {
            compute_stretch_params(
                data.as_ptr(),
                width,
                height,
                input_range,
                &mut c_params as *mut CStretchParams
            );
        }

        c_params.into()
    }

    /// Compute midtones parameter from median and shadows/highlights
    fn compute_midtones(normalized_median: f32, shadows: f32, highlights: f32) -> f32 {
        if highlights == shadows {
            return 0.5;
        }

        // Midtones transformation function (from PixInsight)
        let x = (normalized_median - shadows) / (highlights - shadows);
        if x <= 0.0 {
            0.0
        } else if x >= 1.0 {
            1.0
        } else {
            // Calculate midtones such that f(normalized_median) ≈ 0.5
            // Using the inverse of the midtones transfer function
            ((x - 1.0) * 0.5) / (x - 0.5)
        }
    }

    /// Default parameters (no stretching)
    pub fn default() -> Self {
        StretchParams {
            shadows: 0.0,
            highlights: 1.0,
            midtones: 0.5,
        }
    }

    /// Apply stretch transformation to a single value
    pub fn apply(&self, normalized: f32) -> f32 {
        if normalized <= self.shadows {
            0.0
        } else if normalized >= self.highlights {
            1.0
        } else {
            // Midtones transformation
            let x = (normalized - self.shadows) / (self.highlights - self.shadows);

            if self.midtones == 0.5 {
                // Linear transformation
                x
            } else {
                // Non-linear midtones transformation
                let k1 = self.highlights - self.shadows;
                let k2 = (2.0 * self.midtones - 1.0) / k1;

                (x * k1) / (x * k2 - self.midtones + 1.0)
            }
        }
    }
}

/// Detect the maximum input value based on data type
pub fn calculate_max_input<T>(data: &[T]) -> f32
where
    T: Copy + Into<f32>,
{
    // Try to detect range
    let max_val = data.iter()
        .map(|&v| v.into())
        .fold(f32::NEG_INFINITY, f32::max);

    if max_val <= 255.0 {
        255.0
    } else if max_val <= 65535.0 {
        65535.0
    } else {
        max_val
    }
}

/// Apply auto-stretch to a 2D mono image
pub fn stretch_mono<T>(data: &Array2<T>, params: Option<StretchParams>) -> Array2<u8>
where
    T: Copy + Into<f32> + Send + Sync,
{
    let flat_data: Vec<f32> = data.iter().map(|&v| v.into()).collect();
    let max_input = calculate_max_input(&flat_data);

    let params = params.unwrap_or_else(|| StretchParams::compute_auto(&flat_data, max_input));

    let (height, width) = data.dim();
    let mut result = Array2::zeros((height, width));

    // Apply stretch in parallel
    result.iter_mut()
        .zip(data.iter())
        .for_each(|(out, &val)| {
            let normalized = val.into() / max_input;
            let stretched = params.apply(normalized);
            *out = (stretched * 255.0).round().clamp(0.0, 255.0) as u8;
        });

    result
}

/// Apply auto-stretch to a 3D color image (channels, height, width)
/// Each channel is stretched independently in parallel
pub fn stretch_color<T>(data: &Array3<T>, params: Option<Vec<StretchParams>>) -> Array3<u8>
where
    T: Copy + Into<f32> + Send + Sync,
{
    let (channels, height, width) = data.dim();
    let mut result = Array3::zeros((channels, height, width));

    // Compute parameters for each channel if not provided
    let all_params: Vec<StretchParams> = if let Some(p) = params {
        p
    } else {
        // Compute parameters for each channel
        (0..channels)
            .into_par_iter()
            .map(|c| {
                let channel_data: Vec<f32> = data.index_axis(Axis(0), c)
                    .iter()
                    .map(|&v| v.into())
                    .collect();
                let max_input = calculate_max_input(&channel_data);
                StretchParams::compute_auto(&channel_data, max_input)
            })
            .collect()
    };

    // Apply stretch to each channel using C implementation
    for c in 0..channels {
        let mut channel_data: Vec<f32> = data.index_axis(Axis(0), c)
            .iter()
            .map(|&v| v.into())
            .collect();
        let max_input = calculate_max_input(&channel_data);
        let params = all_params[c];

        // Debug: print stretch parameters
        eprintln!("Channel {}: shadows={:.6}, highlights={:.6}, midtones={:.6}, max_input={:.1}",
                  c, params.shadows, params.highlights, params.midtones, max_input);

        // Convert to C params and apply stretch
        let c_params = CStretchParams {
            shadows: params.shadows,
            highlights: params.highlights,
            midtones: params.midtones,
            max_input: max_input as c_int,
        };

        unsafe {
            apply_stretch(
                channel_data.as_mut_ptr(),
                width as c_int,
                height as c_int,
                &c_params as *const CStretchParams
            );
        }

        // Copy stretched data to output (now in 0-255 range)
        let mut channel_out = result.index_axis_mut(Axis(0), c);
        channel_out.iter_mut()
            .zip(channel_data.iter())
            .for_each(|(out, &val)| {
                *out = val.round().clamp(0.0, 255.0) as u8;
            });
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stretch_params_compute() {
        let data: Vec<f32> = (0..1000).map(|x| x as f32).collect();
        let params = StretchParams::compute_auto(&data, 1000.0);

        assert!(params.shadows >= 0.0 && params.shadows < 1.0);
        assert!(params.highlights > 0.0 && params.highlights <= 1.0);
        assert!(params.highlights > params.shadows);
    }

    #[test]
    fn test_stretch_mono() {
        let data = Array2::from_shape_fn((100, 100), |(y, x)| ((y * 100 + x) * 65) as u16);
        let stretched = stretch_mono(&data, None);

        assert_eq!(stretched.dim(), (100, 100));
        // Check that values are in valid range
        assert!(stretched.iter().all(|&v| v <= 255));
    }

    #[test]
    fn test_stretch_color() {
        let data = Array3::from_shape_fn((3, 100, 100), |(c, y, x)|
            ((c * 10000 + y * 100 + x) * 20) as u16);
        let stretched = stretch_color(&data, None);

        assert_eq!(stretched.dim(), (3, 100, 100));
        assert!(stretched.iter().all(|&v| v <= 255));
    }
}
