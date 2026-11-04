//! Port digest: the same measurements through one `star_measure`, written out in full (`Debug`, timings zeroed), so
//! two copies of the module can be compared byte for byte (port plan phases 2 and 3). Args: <out file> <dir with the
//! rustafits test frames>.
//!
//! Not built by cargo here. To use it, copy it into an `examples/` directory and replace `SM_PATH`:
//! `astroimage::star_measure` in rustafits, `athenaeum_core::star_measure` in `athenaeum-core`; build both with the
//! same toolchain, run both on `tests/`, and `cmp` the two outputs.
use std::io::Write;
use std::path::Path;

use astroimage::{BayerPattern, ImageConverter, PixelData};
use SM_PATH as sm;
use sm::{
    Domain, FitStarsOptions, MeasureOptions, Mode, PixelNoise, PlaneMeasurement, PsfModel, Seed, SeedSpec,
};

struct Rng(u64);
impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    }
    fn gauss(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// Elliptical Moffat stars (β 3) on sky 0.05 with white noise, unit range; peaks above 1 clip at 1.
fn synth(w: usize, h: usize, n: usize, fwhm: f64, ell: f64, theta: f64, max_amp: f64, seed: u64) -> Vec<f32> {
    let mut rng = Rng(seed);
    let mut d = vec![0.05f64; w * h];
    let beta: f64 = 3.0;
    let alpha = fwhm / (2.0 * (2f64.powf(1.0 / beta) - 1.0).sqrt());
    let (ax, ay) = (alpha, alpha * (1.0 - ell));
    let (st, ct) = theta.sin_cos();
    for _ in 0..n {
        let (x, y) = (20.0 + rng.next_f64() * (w as f64 - 40.0), 20.0 + rng.next_f64() * (h as f64 - 40.0));
        let amp = 0.01 + rng.next_f64().powi(3) * max_amp;
        let r = (10.0 * ax).ceil() as i64;
        for yy in (y as i64 - r).max(0)..=(y as i64 + r).min(h as i64 - 1) {
            for xx in (x as i64 - r).max(0)..=(x as i64 + r).min(w as i64 - 1) {
                let (dx, dy) = (xx as f64 - x, yy as f64 - y);
                let (u, v) = (dx * ct + dy * st, -dx * st + dy * ct);
                d[yy as usize * w + xx as usize] += amp * (1.0 + (u / ax).powi(2) + (v / ay).powi(2)).powf(-beta);
            }
        }
    }
    d.iter().map(|&v| (v + 0.002 * rng.gauss()).min(1.0) as f32).collect()
}

fn crop(p: &[f32], w: usize, h: usize, cw: usize, ch: usize) -> (Vec<f32>, usize, usize) {
    let (cw, ch) = (cw.min(w), ch.min(h));
    let (x0, y0) = ((w - cw) / 2, (h - ch) / 2);
    let mut out = Vec::with_capacity(cw * ch);
    for y in y0..y0 + ch {
        out.extend_from_slice(&p[y * w + x0..y * w + x0 + cw]);
    }
    (out, cw, ch)
}

fn real(path: &Path) -> (Vec<f32>, usize, usize) {
    let (meta, px) = ImageConverter::read_raw(path).unwrap();
    let data: Vec<f32> = match px {
        PixelData::Uint16(v) => v.iter().map(|&x| x as f32).collect(),
        PixelData::Float32(v) => v,
    };
    let (w, h) = (meta.width, meta.height);
    let (plane, w, h) = if meta.channels == 3 {
        (data[w * h..2 * w * h].to_vec(), w, h)
    } else if meta.bayer_pattern != BayerPattern::None {
        astroimage::processing::cfa::green_sum(&data, w, h, meta.bayer_pattern)
    } else {
        (data, w, h)
    };
    crop(&plane, w, h, 2048, 1536)
}

fn norm(m: &PlaneMeasurement) -> String {
    let mut m = m.clone();
    m.timings = Default::default();
    m.sub_steps = Default::default();
    format!("{m:?}")
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut out = std::io::BufWriter::new(std::fs::File::create(&args[1]).unwrap());
    let dir = Path::new(&args[2]);
    let unit = Domain::Calibrated { gain_e_per_unit: None };
    let adu = Domain::RawAdu { clip: 65535.0, full_scale: 65535.0, gain_e_per_adu: None };
    let mut planes: Vec<(String, Vec<f32>, usize, usize, Domain)> = vec![
        ("synth-round".into(), synth(1024, 768, 250, 3.0, 0.0, 0.0, 0.6, 7), 1024, 768, unit),
        ("synth-elliptic".into(), synth(1024, 768, 250, 3.4, 0.35, 0.5, 0.6, 8), 1024, 768, unit),
        ("synth-crowded-saturated".into(), synth(1024, 768, 900, 2.6, 0.1, 1.2, 1.4, 9), 1024, 768, unit),
    ];
    for f in ["mono.fits", "osc.fits", "test.xisf"] {
        let (p, w, h) = real(&dir.join(f));
        planes.push((f.to_string(), p, w, h, adu));
    }
    for (name, p, w, h, domain) in &planes {
        let (w, h) = (*w, *h);
        let all = MeasureOptions {
            domain: *domain,
            psf_model: PsfModel::Auto,
            seeds: SeedSpec::MatchedFilter { threshold: 6.0 },
            ..Default::default()
        };
        let m_all = sm::measure_plane(p, w, h, &all, None);
        writeln!(out, "== {name} all_auto\n{}", norm(&m_all)).unwrap();
        let budget = MeasureOptions {
            psf_model: PsfModel::Moffat4,
            mode: Mode::Budget { budget: 500 },
            pixel_noise: PixelNoise::White,
            ..all.clone()
        };
        writeln!(out, "== {name} budget_m4\n{}", norm(&sm::measure_plane(p, w, h, &budget, None))).unwrap();
        let pop = MeasureOptions { mode: Mode::Population { statistic: 200, list: 300 }, ..all.clone() };
        writeln!(out, "== {name} population\n{}", norm(&sm::measure_plane(p, w, h, &pop, None))).unwrap();
        let prep = sm::prepare_plane(p, w, h, &all, None);
        let beta = sm::plane_beta(p, w, h, &all, &prep, None);
        writeln!(out, "== {name} plane_beta {beta:?}").unwrap();
        let fin = sm::finish_plane(p, w, h, &all, prep, None, None, None);
        writeln!(out, "== {name} prepare_finish\n{}", norm(&fin)).unwrap();
        let seeds: Vec<Seed> = m_all.detected.iter().take(300).copied().collect();
        let fopts = FitStarsOptions { trail_verdict: true, defocus_verdict: true, ..Default::default() };
        let mut outcome = sm::fit_stars_at(p, w, h, &seeds, 4.0, &fopts, None);
        writeln!(out, "== {name} fit_stars_at\n{outcome:?}").unwrap();
        let ring = m_all.frame.fwhm_px.unwrap_or(3.0);
        sm::correct_at_ring(p, w, h, &seeds, &mut outcome, ring, &fopts, None);
        writeln!(out, "== {name} correct_at_ring\n{outcome:?}").unwrap();
        eprintln!("{name}: {} fits (all), {} seeds", m_all.stars.len(), m_all.seeds);
    }
}
