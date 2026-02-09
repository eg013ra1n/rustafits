/// 2x2 averaging for mono preview mode.
pub fn bin_2x2_float(input: &[f32], width: usize, height: usize) -> (Vec<f32>, usize, usize) {
    let out_w = width / 2;
    let out_h = height / 2;
    let mut output = vec![0f32; out_w * out_h];

    for y in 0..out_h {
        let in_y = y * 2;
        let row0 = in_y * width;
        let row1 = (in_y + 1) * width;

        for x in 0..out_w {
            let in_x = x * 2;
            let p00 = input[row0 + in_x];
            let p01 = input[row0 + in_x + 1];
            let p10 = input[row1 + in_x];
            let p11 = input[row1 + in_x + 1];
            output[y * out_w + x] = (p00 + p01 + p10 + p11) * 0.25;
        }
    }

    (output, out_w, out_h)
}
