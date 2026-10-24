//! Single-file timing of the background mesh and the MRS noise on a float32 FITS plane.
//! usage: bg_bench <file.fits> [reps]   (run with RAYON_NUM_THREADS=1 for CPU seconds)
use astroimage::analysis::background::{estimate_background_mesh, estimate_noise_mrs};
use std::time::Instant;

fn read_f32_fits(path: &str) -> (Vec<f32>, usize, usize) {
    let b = std::fs::read(path).unwrap();
    let (mut w, mut h, mut bitpix) = (0usize, 0usize, 0i64);
    let mut off = 0;
    'outer: loop {
        for c in 0..36 {
            let card = std::str::from_utf8(&b[off + c * 80..off + c * 80 + 80]).unwrap();
            let v = |k: &str| card.starts_with(k).then(|| card[10..].split('/').next().unwrap().trim().parse::<usize>().unwrap());
            if let Some(x) = v("NAXIS1  ") { w = x; }
            if let Some(x) = v("NAXIS2  ") { h = x; }
            if card.starts_with("BITPIX  ") { bitpix = card[10..].split('/').next().unwrap().trim().parse().unwrap(); }
            if card.starts_with("END ") { off += 2880; break 'outer; }
        }
        off += 2880;
    }
    assert_eq!(bitpix, -32, "bg_bench reads float32 FITS only (BITPIX = -32); {path} has BITPIX = {bitpix}");
    let d = (0..w * h).map(|i| f32::from_be_bytes(b[off + i * 4..off + i * 4 + 4].try_into().unwrap())).collect();
    (d, w, h)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (mut d, w, h) = read_f32_fits(&a[1]);
    for v in d.iter_mut() { *v *= 65535.0; }
    let reps: usize = a.get(2).map(|s| s.parse().unwrap()).unwrap_or(3);
    for _ in 0..reps {
        let t = Instant::now();
        let r = estimate_background_mesh(&d, w, h, 128, false);
        let tm = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let n = estimate_noise_mrs(&d, w, h, 4);
        let tn = t.elapsed().as_secs_f64();
        let m = r.background_map.as_ref().unwrap();
        let mut hsh = 0u64;
        for v in m { hsh = hsh.wrapping_mul(1099511628211).wrapping_add(v.to_bits() as u64); }
        println!("{w}x{h} mesh {tm:.3}s mrs {tn:.3}s  bg={} noise={} mrs={} maphash={hsh:016x}", r.background, r.noise, n);
    }
}
