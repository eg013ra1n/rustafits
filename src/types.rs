use std::collections::HashMap;

#[allow(dead_code)]
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum BayerPattern {
    None,
    Rggb,
    Bggr,
    Gbrg,
    Grbg,
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum DataType {
    Uint16,
    Float32,
}

pub enum PixelData {
    Uint16(Vec<u16>),
    Float32(Vec<f32>),
}

#[derive(Clone, Default)]
pub struct FitsMetadata {
    pub date_obs: Option<String>,
    pub time_obs: Option<String>,
    pub exposure: Option<f64>,
    pub iso: Option<u32>,
    pub gain: Option<f64>,
    pub aperture: Option<f64>,
    pub shutter: Option<f64>,
    pub ccd_temp: Option<f64>,
    pub instrument: Option<String>,
    pub telescope: Option<String>,
    pub lens: Option<String>,
    pub filter: Option<String>,
    pub object: Option<String>,
    pub ra: Option<String>,
    pub dec: Option<String>,
    pub alt: Option<f64>,
    pub az: Option<f64>,
    pub observer: Option<String>,
    pub site: Option<String>,
    pub notes: Option<String>,
    pub bayerpat: Option<String>,
    pub all_keywords: HashMap<String, String>,
}

#[derive(Clone)]
#[allow(dead_code)]
pub struct ImageMetadata {
    pub width: usize,
    pub height: usize,
    pub channels: usize,
    pub dtype: DataType,
    pub bayer_pattern: BayerPattern,
    pub flip_vertical: bool,
    pub observational: FitsMetadata,
}

#[allow(dead_code)]
pub struct ProcessedImage {
    pub data: Vec<u8>,
    pub width: usize,
    pub height: usize,
    pub is_color: bool,
    pub channels: u8,
    pub flip_vertical: bool,
    pub observational: FitsMetadata,
}

#[allow(dead_code)]
pub struct ProcessConfig {
    pub downscale_factor: usize,
    pub jpeg_quality: u8,
    pub apply_debayer: bool,
    pub preview_mode: bool,
    pub auto_stretch: bool,
    pub rgba_output: bool,
}
