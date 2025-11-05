use anyhow::{Context, Result};
use fitsio::FitsFile;
use fitsio::hdu::HduInfo;
use ndarray::{Array2, Array3};
use std::path::Path;

#[derive(Debug, Clone, Copy)]
pub enum BayerPattern {
    RGGB,
    BGGR,
    GBRG,
    GRBG,
}

impl BayerPattern {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_uppercase().as_str() {
            "RGGB" => Some(BayerPattern::RGGB),
            "BGGR" => Some(BayerPattern::BGGR),
            "GBRG" => Some(BayerPattern::GBRG),
            "GRBG" => Some(BayerPattern::GRBG),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum RowOrder {
    TopDown,
    BottomUp,
}

pub struct FitsMetadata {
    pub width: usize,
    pub height: usize,
    pub channels: usize,
    pub bayer_pattern: Option<BayerPattern>,
    pub row_order: RowOrder,
}

pub enum FitsImageData {
    U16(Array2<u16>),
    F32(Array2<f32>),
    U16Color(Array3<u16>),
    F32Color(Array3<f32>),
}

pub struct FitsImage {
    pub metadata: FitsMetadata,
    pub data: FitsImageData,
}

impl FitsImage {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut fptr = FitsFile::open(path.as_ref())
            .context("Failed to open FITS file")?;

        let hdu = fptr.primary_hdu()
            .context("Failed to get primary HDU")?;

        // Get HDU info to determine dimensions and type
        let (width, height, channels) = match &hdu.info {
            HduInfo::ImageInfo { shape, .. } => {
                match shape.len() {
                    2 => (shape[1], shape[0], 1),
                    3 => (shape[2], shape[1], shape[0]),
                    _ => anyhow::bail!("Unsupported FITS image dimensions: {:?}", shape),
                }
            }
            _ => anyhow::bail!("Primary HDU is not an image"),
        };

        // Read header keywords
        let bayer_pattern = hdu.read_key::<String>(&mut fptr, "BAYERPAT")
            .ok()
            .and_then(|s| BayerPattern::from_str(&s));

        let row_order = hdu.read_key::<String>(&mut fptr, "ROWORDER")
            .ok()
            .map(|s| match s.to_uppercase().as_str() {
                "TOP-DOWN" => RowOrder::TopDown,
                _ => RowOrder::BottomUp,
            })
            .unwrap_or(RowOrder::BottomUp);

        let metadata = FitsMetadata {
            width,
            height,
            channels,
            bayer_pattern,
            row_order,
        };

        // Read image data based on BITPIX
        let data = if channels == 1 {
            // Try reading as different types
            if let Ok(raw_data) = hdu.read_image::<Vec<u16>>(&mut fptr) {
                let array = Array2::from_shape_vec((height, width), raw_data)
                    .context("Failed to reshape u16 data")?;
                FitsImageData::U16(array)
            } else if let Ok(raw_data) = hdu.read_image::<Vec<f32>>(&mut fptr) {
                let array = Array2::from_shape_vec((height, width), raw_data)
                    .context("Failed to reshape f32 data")?;
                FitsImageData::F32(array)
            } else {
                anyhow::bail!("Unsupported FITS data type");
            }
        } else {
            // Color image (3 channels)
            if let Ok(raw_data) = hdu.read_image::<Vec<u16>>(&mut fptr) {
                let array = Array3::from_shape_vec((channels, height, width), raw_data)
                    .context("Failed to reshape u16 color data")?;
                FitsImageData::U16Color(array)
            } else if let Ok(raw_data) = hdu.read_image::<Vec<f32>>(&mut fptr) {
                let array = Array3::from_shape_vec((channels, height, width), raw_data)
                    .context("Failed to reshape f32 color data")?;
                FitsImageData::F32Color(array)
            } else {
                anyhow::bail!("Unsupported FITS color data type");
            }
        };

        Ok(FitsImage { metadata, data })
    }

    pub fn width(&self) -> usize {
        self.metadata.width
    }

    pub fn height(&self) -> usize {
        self.metadata.height
    }

    pub fn channels(&self) -> usize {
        self.metadata.channels
    }

    pub fn has_bayer_pattern(&self) -> bool {
        self.metadata.bayer_pattern.is_some()
    }
}
