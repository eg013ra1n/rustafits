// Minimal reproduction: trailing-MCU padding divergence at 4:2:0.
//
// Cargo.toml:
//   [dependencies]
//   libjpeg-turbo-rs = "0.6"
//   turbojpeg        = "1.5"   # C libjpeg-turbo, the reference
//
// Run with `cargo run --release`.

fn synth(w: usize, h: usize) -> Vec<u8> {
    let mut v = vec![0u8; w * h * 3];
    let mut s: u32 = 0x0123_4567;
    for p in v.chunks_exact_mut(3) {
        s ^= s << 13;
        s ^= s >> 17;
        s ^= s << 5;
        let n = (s % 200) as u8 + 20;
        p[0] = n;
        p[1] = n.wrapping_add(30);
        p[2] = n.wrapping_add(60);
    }
    v
}

/// Expand to `wp` columns by replicating the last real column (libjpeg's
/// `expand_right_edge` semantics, applied in RGB space).
fn pad_replicate(d: &[u8], w: usize, h: usize, wp: usize) -> Vec<u8> {
    let mut o = vec![0u8; wp * h * 3];
    for y in 0..h {
        for x in 0..wp {
            let sx = x.min(w - 1);
            o[(y * wp + x) * 3..][..3].copy_from_slice(&d[(y * w + sx) * 3..][..3]);
        }
    }
    o
}

fn scan(j: &[u8]) -> &[u8] {
    let i = j.windows(2).position(|w| w == [0xFF, 0xDA]).unwrap();
    &j[i..]
}

fn c_enc(d: &[u8], w: usize, h: usize, s: turbojpeg::Subsamp) -> Vec<u8> {
    let img = turbojpeg::Image { pixels: d, width: w, pitch: w * 3, height: h,
                                 format: turbojpeg::PixelFormat::RGB };
    turbojpeg::compress(img, 90, s).unwrap().to_vec()
}

fn rs_enc(d: &[u8], w: usize, h: usize, s: libjpeg_turbo_rs::Subsampling) -> Vec<u8> {
    libjpeg_turbo_rs::compress(d, w, h, libjpeg_turbo_rs::PixelFormat::Rgb, 90, s).unwrap()
}

fn main() {
    const H: usize = 320;

    println!("width sweep, height {H} — '=' means byte-identical to C libjpeg-turbo\n");
    for (name, cs, rs) in [
        ("4:4:4", turbojpeg::Subsamp::None,   libjpeg_turbo_rs::Subsampling::S444),
        ("4:2:2", turbojpeg::Subsamp::Sub2x1, libjpeg_turbo_rs::Subsampling::S422),
        ("4:2:0", turbojpeg::Subsamp::Sub2x2, libjpeg_turbo_rs::Subsampling::S420),
    ] {
        let marks: String = (320..=336)
            .map(|w| {
                let d = synth(w, H);
                if c_enc(&d, w, H, cs) == rs_enc(&d, w, H, rs) { '=' } else { '#' }
            })
            .collect();
        println!("  {name}  W=320..336  {marks}");
    }
    println!("        w%16       01234567890123450\n");

    println!("the divergent output equals the C encoder fed a pre-padded image:\n");
    for w in [322usize, 324, 326, 328] {
        let wp = w.div_ceil(16) * 16;
        let d = synth(w, H);
        let padded = pad_replicate(&d, w, H, wp);
        let rs_native = rs_enc(&d, w, H, libjpeg_turbo_rs::Subsampling::S420);
        let c_padded = c_enc(&padded, wp, H, turbojpeg::Subsamp::Sub2x2);
        let c_native = c_enc(&d, w, H, turbojpeg::Subsamp::Sub2x2);
        println!(
            "  W={w} (w%16={:>2})  rust(W).scan == C(replicate-padded to {wp}).scan : {}   \
             [{} B vs C native {} B, +{:.3}%]",
            w % 16,
            scan(&rs_native) == scan(&c_padded),
            rs_native.len(),
            c_native.len(),
            100.0 * (rs_native.len() as f64 - c_native.len() as f64) / c_native.len() as f64
        );
    }
}
