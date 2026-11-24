//! Decode-hash harness: a change to a reader must not move one decoded bit.
//!
//! ```text
//! cargo run --release --example decode_hashes -- [--synthetic <dir>] FILES_OR_DIRS...
//! ```
//!
//! For every file one line `<file>\tread_raw\t<sha256>` hashes the layout (width, height, channels, sample type,
//! Bayer pattern, flip) and every pixel's bytes as `ImageConverter::read_raw` returns them. Then one line per
//! converter mode, `<file>\t<mode>\t<sha256>`, hashes the processed image (size, colour flag, channel count, flip and
//! every output byte) for `plain`, `downscale2`, `downscale3`, `preview`, `no-debayer`, `vng`, `rgb` and `rgba`.
//! `rgb` runs the default converter over pre-read data (`process_data`), `plain` over the path, so the two cover the
//! file read and the processing apart. A read or a mode that fails prints `error` and a hash of its message.
//!
//! `--synthetic <dir>` first writes 20 synthetic FITS files into `<dir>` and puts them in front of `FILES`:
//! every BITPIX with BZERO 0 and 32768, ROWORDER BOTTOM-UP and TOP-DOWN, widths with `w % 8` = 0 to 7 (the odd ones
//! are OSC), an OSC file with Bayer offsets, f32 with NaN and f32 RGB. A directory among `FILES` stands for its
//! `.fit`, `.fits`, `.fts` and `.xisf` files, sorted by name.
//!
//! The same file builds against an older checkout: it uses the public API only.
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use astroimage::{ImageConverter, PixelData};

// ----- SHA-256 (FIPS 180-4), kept here so the harness adds no dependency -----

const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

struct Sha256 {
    h: [u32; 8],
    buf: [u8; 64],
    fill: usize,
    len: u64,
}

impl Sha256 {
    fn new() -> Self {
        Sha256 {
            h: [
                0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
                0x5be0cd19,
            ],
            buf: [0; 64],
            fill: 0,
            len: 0,
        }
    }

    fn block(&mut self, b: &[u8; 64]) {
        let mut w = [0u32; 64];
        for (i, c) in b.chunks_exact(4).enumerate() {
            w[i] = u32::from_be_bytes([c[0], c[1], c[2], c[3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let [mut a, mut b_, mut c, mut d, mut e, mut f, mut g, mut h] = self.h;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ (!e & g);
            let t1 = h
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b_) ^ (a & c) ^ (b_ & c);
            let t2 = s0.wrapping_add(maj);
            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b_;
            b_ = a;
            a = t1.wrapping_add(t2);
        }
        for (s, v) in self.h.iter_mut().zip([a, b_, c, d, e, f, g, h]) {
            *s = s.wrapping_add(v);
        }
    }

    fn update(&mut self, mut data: &[u8]) {
        self.len += data.len() as u64;
        if self.fill > 0 {
            let take = (64 - self.fill).min(data.len());
            self.buf[self.fill..self.fill + take].copy_from_slice(&data[..take]);
            self.fill += take;
            data = &data[take..];
            if self.fill == 64 {
                let b = self.buf;
                self.block(&b);
                self.fill = 0;
            }
        }
        while data.len() >= 64 {
            let (head, rest) = data.split_at(64);
            self.block(head.try_into().unwrap());
            data = rest;
        }
        self.buf[..data.len()].copy_from_slice(data);
        self.fill += data.len();
    }

    fn finish(mut self) -> String {
        let bits = self.len * 8;
        let mut pad = vec![0x80u8];
        while (self.fill + pad.len()) % 64 != 56 {
            pad.push(0);
        }
        pad.extend_from_slice(&bits.to_be_bytes());
        let len = self.len;
        self.update(&pad);
        self.len = len;
        self.h.iter().fold(String::new(), |mut s, v| {
            let _ = write!(s, "{v:08x}");
            s
        })
    }
}

fn sha_of(parts: &[&[u8]]) -> String {
    let mut s = Sha256::new();
    for p in parts {
        s.update(p);
    }
    s.finish()
}

// ----- hashing a decode -----

fn pixel_bytes(p: &PixelData) -> Vec<u8> {
    match p {
        PixelData::Uint16(v) => v.iter().flat_map(|x| x.to_le_bytes()).collect(),
        PixelData::Float32(v) => v.iter().flat_map(|x| x.to_bits().to_le_bytes()).collect(),
    }
}

fn read_raw_hash(path: &Path) -> Result<String, String> {
    let (meta, pixels) = ImageConverter::read_raw(path).map_err(|e| format!("{e:#}"))?;
    let layout = format!(
        "{}x{}x{} {:?} {:?} flip={}",
        meta.width, meta.height, meta.channels, meta.dtype, meta.bayer_pattern, meta.flip_vertical
    );
    Ok(sha_of(&[layout.as_bytes(), &pixel_bytes(&pixels)]))
}

fn processed_hash(img: &astroimage::ProcessedImage) -> String {
    let head = format!(
        "{}x{} color={} ch={} flip={}",
        img.width, img.height, img.is_color, img.channels, img.flip_vertical
    );
    sha_of(&[head.as_bytes(), &img.data])
}

/// The eight converter modes, in output order.
const MODES: [&str; 8] = [
    "plain",
    "downscale2",
    "downscale3",
    "preview",
    "no-debayer",
    "vng",
    "rgb",
    "rgba",
];

fn mode_hash(mode: &str, path: &Path) -> Result<String, String> {
    let e = |e: anyhow::Error| format!("{e:#}");
    let img = match mode {
        "plain" => ImageConverter::new().process(path).map_err(e)?,
        "downscale2" => ImageConverter::new()
            .with_downscale(2)
            .process(path)
            .map_err(e)?,
        "downscale3" => ImageConverter::new()
            .with_downscale(3)
            .process(path)
            .map_err(e)?,
        "preview" => ImageConverter::new()
            .with_preview_mode()
            .process(path)
            .map_err(e)?,
        "no-debayer" => ImageConverter::new()
            .without_debayer()
            .process(path)
            .map_err(e)?,
        "vng" => ImageConverter::new()
            .with_vng_debayer()
            .process(path)
            .map_err(e)?,
        "rgb" => {
            let (meta, pixels) = ImageConverter::read_raw(path).map_err(e)?;
            ImageConverter::new().process_data(meta, pixels).map_err(e)?
        }
        "rgba" => ImageConverter::new()
            .with_rgba_output()
            .process(path)
            .map_err(e)?,
        other => return Err(format!("unknown mode {other}")),
    };
    Ok(processed_hash(&img))
}

fn line(file: &str, what: &str, r: Result<String, String>) {
    match r {
        Ok(h) => println!("{file}\t{what}\t{h}"),
        Err(m) => {
            let first = m.lines().next().unwrap_or("");
            println!("{file}\t{what}\terror {}", sha_of(&[first.as_bytes()]));
        }
    }
}

// ----- the synthetic files -----

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    /// A float in [0, 1).
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A smooth background plus noise plus a few bright blobs, in [0, 1): the frame an auto-stretch can work on.
fn scene(w: usize, h: usize, planes: usize, seed: u64) -> Vec<f64> {
    let mut rng = Rng(seed | 1);
    let blobs: Vec<(f64, f64, f64)> = (0..12)
        .map(|_| {
            (
                rng.unit() * w as f64,
                rng.unit() * h as f64,
                0.4 + rng.unit() * 0.6,
            )
        })
        .collect();
    let mut out = Vec::with_capacity(w * h * planes);
    for p in 0..planes {
        for y in 0..h {
            for x in 0..w {
                let mut v = 0.08 + 0.05 * (x as f64 / w as f64) + 0.02 * p as f64 + 0.01 * rng.unit();
                for &(bx, by, a) in &blobs {
                    let d2 = (x as f64 - bx).powi(2) + (y as f64 - by).powi(2);
                    v += a * (-d2 / 3.0).exp();
                }
                out.push(v.min(0.999));
            }
        }
    }
    out
}

fn card(key: &str, value: &str) -> [u8; 80] {
    let mut c = [b' '; 80];
    let text = format!("{key:<8}= {value:>20}");
    c[..text.len()].copy_from_slice(text.as_bytes());
    c
}

struct Spec {
    name: String,
    bitpix: i32,
    bzero: Option<f64>,
    w: usize,
    h: usize,
    planes: usize,
    bayer: Option<&'static str>,
    offsets: Option<(i32, i32)>,
    roworder: Option<&'static str>,
    nan: bool,
}

impl Spec {
    fn mono(name: &str, bitpix: i32, bzero: Option<f64>, w: usize, h: usize) -> Self {
        Spec {
            name: name.to_string(),
            bitpix,
            bzero,
            w,
            h,
            planes: 1,
            bayer: None,
            offsets: None,
            roworder: None,
            nan: false,
        }
    }
}

fn write_fits(dir: &Path, s: &Spec, seed: u64) -> std::io::Result<PathBuf> {
    let mut cards: Vec<[u8; 80]> = vec![
        card("SIMPLE", "T"),
        card("BITPIX", &s.bitpix.to_string()),
        card("NAXIS", &(if s.planes > 1 { 3 } else { 2 }).to_string()),
        card("NAXIS1", &s.w.to_string()),
        card("NAXIS2", &s.h.to_string()),
    ];
    if s.planes > 1 {
        cards.push(card("NAXIS3", &s.planes.to_string()));
    }
    if let Some(z) = s.bzero {
        cards.push(card("BZERO", &format!("{z:.1}")));
        cards.push(card("BSCALE", "1.0"));
    }
    if let Some(b) = s.bayer {
        cards.push(card("BAYERPAT", &format!("'{b}'")));
    }
    if let Some((x, y)) = s.offsets {
        cards.push(card("XBAYROFF", &x.to_string()));
        cards.push(card("YBAYROFF", &y.to_string()));
    }
    if let Some(r) = s.roworder {
        cards.push(card("ROWORDER", &format!("'{r}'")));
    }
    let mut end = [b' '; 80];
    end[..3].copy_from_slice(b"END");
    cards.push(end);
    let mut bytes: Vec<u8> = cards.concat();
    while bytes.len() % 2880 != 0 {
        bytes.push(b' ');
    }

    let mut vals = scene(s.w, s.h, s.planes, seed);
    if s.nan {
        for (i, v) in vals.iter_mut().enumerate() {
            if i % 37 == 5 {
                *v = f64::NAN;
            }
        }
    }
    let bzero = s.bzero.unwrap_or(0.0);
    for v in &vals {
        match s.bitpix {
            8 => bytes.push((v * 255.0).round() as u8),
            16 => {
                // Physical value v * 65535 stored as value - BZERO (signed 16-bit).
                let phys = (v * 65535.0).round();
                let stored = if bzero == 0.0 { (phys / 2.0).round() } else { phys - bzero };
                bytes.extend_from_slice(&(stored as i16).to_be_bytes());
            }
            32 => {
                let phys = (v * 65535.0).round();
                bytes.extend_from_slice(&((phys - bzero) as i32).to_be_bytes());
            }
            -32 => bytes.extend_from_slice(&((*v as f32) + bzero as f32).to_bits().to_be_bytes()),
            -64 => bytes.extend_from_slice(&(*v + bzero).to_bits().to_be_bytes()),
            other => panic!("bitpix {other}"),
        }
    }
    while bytes.len() % 2880 != 0 {
        bytes.push(0);
    }
    let path = dir.join(format!("{}.fits", s.name));
    std::fs::write(&path, bytes)?;
    Ok(path)
}

fn synthetic_specs() -> Vec<Spec> {
    let mut v = Vec::new();
    // Every BITPIX with BZERO 0 (written without the card) and with BZERO 32768.
    for (bp, tag) in [(8, "u8"), (16, "i16"), (32, "i32"), (-32, "f32"), (-64, "f64")] {
        v.push(Spec::mono(&format!("bitpix_{tag}_bz0"), bp, None, 47, 31));
    }
    for (bp, tag) in [(16, "i16"), (32, "i32")] {
        v.push(Spec::mono(&format!("bitpix_{tag}_bz32768"), bp, Some(32768.0), 47, 31));
    }
    // Both row orders, as OSC so the flip meets the debayer.
    for (name, order) in [("order_bottom_up", "BOTTOM-UP"), ("order_top_down", "TOP-DOWN")] {
        let mut s = Spec::mono(name, 16, Some(32768.0), 34, 18);
        s.bayer = Some("RGGB");
        s.roworder = Some(order);
        v.push(s);
    }
    // Widths with w % 8 = 0..7; the odd ones are OSC.
    for k in 0..8usize {
        let w = 48 + k;
        let mut s = Spec::mono(&format!("width_{w}"), 16, Some(32768.0), w, 20);
        if w % 2 == 1 {
            s.bayer = Some(["RGGB", "BGGR", "GRBG", "GBRG"][(k / 2) % 4]);
            s.roworder = Some("TOP-DOWN");
        }
        v.push(s);
    }
    // OSC with Bayer offsets.
    let mut s = Spec::mono("osc_offsets", 16, Some(32768.0), 40, 22);
    s.bayer = Some("RGGB");
    s.offsets = Some((1, 1));
    s.roworder = Some("TOP-DOWN");
    v.push(s);
    // f32 with NaN, and f32 RGB.
    let mut s = Spec::mono("f32_nan", -32, None, 45, 29);
    s.nan = true;
    v.push(s);
    let mut s = Spec::mono("f32_rgb", -32, None, 38, 26);
    s.planes = 3;
    v.push(s);
    v
}

fn write_synthetic(dir: &Path) -> std::io::Result<Vec<PathBuf>> {
    std::fs::create_dir_all(dir)?;
    synthetic_specs()
        .iter()
        .enumerate()
        .map(|(i, s)| write_fits(dir, s, 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(i as u64 + 1)))
        .collect()
}

fn expand(arg: &str) -> Vec<PathBuf> {
    let p = PathBuf::from(arg);
    if !p.is_dir() {
        return vec![p];
    }
    let mut files: Vec<PathBuf> = std::fs::read_dir(&p)
        .unwrap_or_else(|e| panic!("{arg}: {e}"))
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|f| {
            f.extension()
                .and_then(|e| e.to_str())
                .map(|e| {
                    matches!(
                        e.to_ascii_lowercase().as_str(),
                        "fit" | "fits" | "fts" | "xisf"
                    )
                })
                .unwrap_or(false)
        })
        .collect();
    files.sort();
    files
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut files: Vec<PathBuf> = Vec::new();
    if args.first().map(String::as_str) == Some("--synthetic") {
        if args.len() < 2 {
            eprintln!("usage: decode_hashes [--synthetic <dir>] FILES...");
            std::process::exit(2);
        }
        let dir = PathBuf::from(&args[1]);
        files = write_synthetic(&dir).unwrap_or_else(|e| panic!("{}: {e}", dir.display()));
        args.drain(..2);
    }
    for a in &args {
        files.extend(expand(a));
    }
    for f in &files {
        let name = f.display().to_string();
        line(&name, "read_raw", read_raw_hash(f));
        for mode in MODES {
            line(&name, mode, mode_hash(mode, f));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_known_answers() {
        assert_eq!(
            sha_of(&[b"abc"]),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(
            sha_of(&[b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"]),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
        // Split input across updates equals one update.
        let data = vec![0x5au8; 1000];
        assert_eq!(
            sha_of(&[&data[..7], &data[7..333], &data[333..]]),
            "8fe15844cfeedd35f5dc30a9fa5ed38afd849dbe4f8dcae5642d934be0afb13d"
        );
    }
}
