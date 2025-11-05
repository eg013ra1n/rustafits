use ndarray::{Array2, Array3};
use num_traits::Zero;

/// Downscale a 2D array using nearest-neighbor sampling
/// This is fast but sacrifices quality - matches QuickLook implementation
pub fn downscale_mono<T: Copy + Clone + Zero>(data: &Array2<T>, factor: usize) -> Array2<T> {
    if factor <= 1 {
        return data.clone();
    }

    let (height, width) = data.dim();
    let new_height = height / factor;
    let new_width = width / factor;

    let mut result = Array2::zeros((new_height, new_width));

    for y in 0..new_height {
        for x in 0..new_width {
            result[[y, x]] = data[[y * factor, x * factor]];
        }
    }

    result
}

/// Downscale a 3D array (channels, height, width) using nearest-neighbor sampling
pub fn downscale_color<T: Copy + Clone + Zero>(data: &Array3<T>, factor: usize) -> Array3<T> {
    if factor <= 1 {
        return data.clone();
    }

    let (channels, height, width) = data.dim();
    let new_height = height / factor;
    let new_width = width / factor;

    let mut result = Array3::zeros((channels, new_height, new_width));

    for c in 0..channels {
        for y in 0..new_height {
            for x in 0..new_width {
                result[[c, y, x]] = data[[c, y * factor, x * factor]];
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_downscale_mono() {
        let data = Array2::from_shape_fn((100, 100), |(y, x)| (y * 100 + x) as u16);
        let downscaled = downscale_mono(&data, 2);
        assert_eq!(downscaled.dim(), (50, 50));
        assert_eq!(downscaled[[0, 0]], data[[0, 0]]);
        assert_eq!(downscaled[[1, 1]], data[[2, 2]]);
    }

    #[test]
    fn test_downscale_color() {
        let data = Array3::from_shape_fn((3, 100, 100), |(c, y, x)| (c * 10000 + y * 100 + x) as u16);
        let downscaled = downscale_color(&data, 2);
        assert_eq!(downscaled.dim(), (3, 50, 50));
        assert_eq!(downscaled[[0, 0, 0]], data[[0, 0, 0]]);
        assert_eq!(downscaled[[1, 1, 1]], data[[1, 2, 2]]);
    }
}
