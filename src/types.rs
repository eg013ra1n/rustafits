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

/// Container format of an input file.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum FileFormat {
    Fits,
    Xisf,
}

/// One header keyword: a FITS card, an XISF `<FITSKeyword>`, or an XISF
/// `<Property>` (key = its id, e.g. `Instrument:ExposureTime`). String values
/// are unquoted; commentary cards (COMMENT, HISTORY) carry their text as value.
#[derive(Clone, Debug, PartialEq)]
pub struct HeaderKeyword {
    pub key: String,
    pub value: String,
    pub comment: String,
}

/// What a file's header says, read without touching the pixel data.
/// `width`, `height`, `channels`, `bayer_pattern` and `flip_vertical` are the
/// values `ImageConverter::read_raw` would report for the same file.
#[derive(Clone, Debug)]
pub struct ImageHeader {
    pub format: FileFormat,
    pub width: usize,
    pub height: usize,
    pub channels: usize,
    /// Stored sample type: `u8`, `u16`, `i16`, `u32`, `i32`, `f32` or `f64`.
    pub sample_type: &'static str,
    pub bayer_pattern: BayerPattern,
    pub flip_vertical: bool,
    /// Every keyword in file order.
    pub keywords: Vec<HeaderKeyword>,
}

impl ImageHeader {
    /// Value of the first keyword named `key` (case-insensitive).
    pub fn get(&self, key: &str) -> Option<&str> {
        self.keywords
            .iter()
            .find(|k| k.key.eq_ignore_ascii_case(key))
            .map(|k| k.value.as_str())
    }
}

pub enum PixelData {
    Uint16(Vec<u16>),
    Float32(Vec<f32>),
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
}

#[allow(dead_code)]
pub struct ProcessedImage {
    pub data: Vec<u8>,
    pub width: usize,
    pub height: usize,
    pub is_color: bool,
    /// Number of channels in `data`: 3 = RGB, 4 = RGBA.
    pub channels: u8,
    /// Whether the image was vertically flipped during processing.
    pub flip_vertical: bool,
}

#[allow(dead_code)]
pub struct ProcessConfig {
    pub downscale_factor: usize,
    pub jpeg_quality: u8,
    pub apply_debayer: bool,
    pub preview_mode: bool,
    pub auto_stretch: bool,
    /// Output RGBA (4 bytes/pixel) instead of RGB (3 bytes/pixel).
    pub rgba_output: bool,
    /// Debayer a CFA frame with the full-resolution gradient method
    /// ([`crate::processing::vng`]) instead of the super-pixel one. Keeps the
    /// native pixel grid — four times the pixels, and roughly a hundred times
    /// the debayer cost — so it is opt-in and meant for one frame a user is
    /// actually looking at, not for bulk conversion.
    pub vng_debayer: bool,
}
