/// Replicate grayscale u8 to interleaved RGB u8.
pub fn replicate_gray_to_rgb(gray: &[u8]) -> Vec<u8> {
    let mut rgb = vec![0u8; gray.len() * 3];
    for (i, &val) in gray.iter().enumerate() {
        rgb[i * 3] = val;
        rgb[i * 3 + 1] = val;
        rgb[i * 3 + 2] = val;
    }
    rgb
}

/// Vertical flip of an RGB image (3 bytes per pixel).
pub fn vertical_flip_rgb(data: &mut [u8], width: usize, height: usize) {
    let row_bytes = width * 3;
    let mut temp = vec![0u8; row_bytes];
    for y in 0..height / 2 {
        let top = y * row_bytes;
        let bot = (height - 1 - y) * row_bytes;
        temp.copy_from_slice(&data[top..top + row_bytes]);
        data.copy_within(bot..bot + row_bytes, top);
        data[bot..bot + row_bytes].copy_from_slice(&temp);
    }
}

/// Convert u16 slice to f32.
pub fn u16_to_f32(data: &[u16]) -> Vec<f32> {
    data.iter().map(|&v| v as f32).collect()
}
