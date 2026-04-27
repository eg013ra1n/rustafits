/// rustafits-debug: diagnostic CLI for the analysis pipeline.
///
/// Built only with `--features debug-pipeline`. Zero cost to the production binary.

use std::path::PathBuf;
use std::time::Instant;

use anyhow::{bail, Context, Result};

use astroimage::analysis::background;
use astroimage::analysis::convolution;
use astroimage::analysis::detection::{self, DetectedStar, DetectionParams};
use astroimage::analysis::fitting::{self, PixelSample,
    CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS};
use astroimage::analysis::metrics;
use astroimage::analysis::render;
use astroimage::analysis::snr;
use astroimage::analysis::{estimate_fwhm_from_stars, prepare_luminance, sigma_clipped_median};
use astroimage::ImageAnalyzer;
use astroimage::formats;

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

// ── Options ──────────────────────────────────────────────────────────────────

struct Opts {
    subcommand: String,
    input: PathBuf,
    sigma: f32,
    top: usize,
    fwhm_override: Option<f32>,
    save_dir: Option<PathBuf>,
    star_pos: Option<(f32, f32)>,
    stamp_radius: usize,
    near: Option<(f32, f32, f32)>,
    no_debayer: bool,
    mrs_layers: usize,
    measure_cap: usize,
    fit_max_iter: usize,
    fit_tolerance: f64,
    fit_max_rejects: usize,
}

fn print_usage() {
    eprintln!("rustafits-debug: diagnostic CLI for the analysis pipeline");
    eprintln!();
    eprintln!("Usage: rustafits-debug <subcommand> <input> [OPTIONS]");
    eprintln!();
    eprintln!("Subcommands:");
    eprintln!("  background   Background + noise estimation");
    eprintln!("  detect       Star detection (show all candidates)");
    eprintln!("  measure      PSF measurement (FWHM, ecc, HFR per star)");
    eprintln!("  fit          Single-star PSF fit (Gaussian + Moffat side by side)");
    eprintln!("  query        Find nearest detected star and show full metrics + fits");
    eprintln!("  pipeline     Full pipeline with step-by-step verbose output + timing");
    eprintln!("  dump         TSV output of all measured stars");
    eprintln!("  compare      Compare against PixInsight SubframeSelector CSV");
    eprintln!();
    eprintln!("Common options:");
    eprintln!("  --sigma <F>         Detection sigma (default: 5.0)");
    eprintln!("  --top <N>           Show top N stars (default: 20)");
    eprintln!("  --fwhm <F>          Override initial kernel FWHM");
    eprintln!("  --save <dir>        Save debug images to directory");
    eprintln!("  --star <x>,<y>      Star position for 'fit'/'query' subcommand");
    eprintln!("  --radius <N>        Stamp radius for 'fit'/'query' (default: 15)");
    eprintln!("  --near <x>,<y>,<r>  Filter stars within radius of position");
    eprintln!("  --no-debayer        Skip green interpolation for OSC");
    eprintln!("  --mrs <N>           MRS wavelet noise layers (default: 4)");
    eprintln!("  --measure-cap <N>   Max stars to measure (default: 2000, 0 = all)");
    eprintln!("  --fit-iter <N>      LM max iterations (default: 25)");
    eprintln!("  --fit-tol <F>       LM convergence tolerance (default: 1e-4)");
    eprintln!("  --fit-rejects <N>   LM consecutive reject bailout (default: 5)");
}

fn parse_args() -> Result<Opts> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        print_usage();
        bail!("Missing required arguments");
    }

    let subcommand = args[1].clone();
    let input = PathBuf::from(&args[2]);

    let mut opts = Opts {
        subcommand,
        input,
        sigma: 5.0,
        top: 20,
        fwhm_override: None,
        save_dir: None,
        star_pos: None,
        stamp_radius: 15,
        near: None,
        no_debayer: false,
        mrs_layers: 0,
        measure_cap: 2000,
        fit_max_iter: 25,
        fit_tolerance: 1e-4,
        fit_max_rejects: 5,
    };

    let mut i = 3;
    while i < args.len() {
        match args[i].as_str() {
            "--sigma" => {
                i += 1;
                opts.sigma = args[i].parse().context("--sigma value")?;
            }
            "--top" => {
                i += 1;
                opts.top = args[i].parse().context("--top value")?;
            }
            "--fwhm" => {
                i += 1;
                opts.fwhm_override = Some(args[i].parse().context("--fwhm value")?);
            }
            "--save" => {
                i += 1;
                opts.save_dir = Some(PathBuf::from(&args[i]));
            }
            "--star" => {
                i += 1;
                let parts: Vec<&str> = args[i].split(',').collect();
                if parts.len() != 2 {
                    bail!("--star expects x,y (e.g. --star 100,200)");
                }
                let x: f32 = parts[0].parse().context("star x")?;
                let y: f32 = parts[1].parse().context("star y")?;
                opts.star_pos = Some((x, y));
            }
            "--radius" => {
                i += 1;
                opts.stamp_radius = args[i].parse().context("--radius value")?;
            }
            "--near" => {
                i += 1;
                let parts: Vec<&str> = args[i].split(',').collect();
                if parts.len() != 3 {
                    bail!("--near expects x,y,r (e.g. --near 100,200,50)");
                }
                let x: f32 = parts[0].parse().context("near x")?;
                let y: f32 = parts[1].parse().context("near y")?;
                let r: f32 = parts[2].parse().context("near radius")?;
                opts.near = Some((x, y, r));
            }
            "--no-debayer" => {
                opts.no_debayer = true;
            }
            "--mrs" => {
                i += 1;
                opts.mrs_layers = args[i].parse().context("--mrs value")?;
            }
            "--measure-cap" => {
                i += 1;
                opts.measure_cap = args[i].parse().context("--measure-cap value")?;
            }
            "--fit-iter" => {
                i += 1;
                opts.fit_max_iter = args[i].parse().context("--fit-iter value")?;
            }
            "--fit-tol" => {
                i += 1;
                opts.fit_tolerance = args[i].parse().context("--fit-tol value")?;
            }
            "--fit-rejects" => {
                i += 1;
                opts.fit_max_rejects = args[i].parse().context("--fit-rejects value")?;
            }
            other => {
                bail!("Unknown option: {}", other);
            }
        }
        i += 1;
    }

    Ok(opts)
}

// ── Main dispatch ────────────────────────────────────────────────────────────

fn run() -> Result<()> {
    let opts = parse_args()?;

    // `compare` handles its own file I/O (reads multiple files from PI CSV)
    if opts.subcommand == "compare" {
        return cmd_compare(&opts);
    }

    // Read image
    let t0 = Instant::now();
    let (meta, pixels) =
        formats::read_image(opts.input.as_ref()).context("Failed to read image")?;
    eprintln!(
        "Read {}x{} {:?} ({} ch) in {:.1}ms",
        meta.width,
        meta.height,
        opts.input.file_name().unwrap_or_default().to_string_lossy(),
        meta.channels,
        t0.elapsed().as_secs_f64() * 1000.0
    );

    // Prepare luminance
    let t1 = Instant::now();
    let apply_debayer = !opts.no_debayer;
    let (lum, width, height, channels, green_mask) =
        prepare_luminance(&meta, &pixels, apply_debayer);
    eprintln!(
        "Luminance extraction: {:.1}ms ({}x{}, {} ch, debayer={})",
        t1.elapsed().as_secs_f64() * 1000.0,
        width,
        height,
        channels,
        apply_debayer && meta.bayer_pattern != astroimage::BayerPattern::None
    );

    if let Some(ref dir) = opts.save_dir {
        std::fs::create_dir_all(dir).context("Failed to create save directory")?;
    }

    match opts.subcommand.as_str() {
        "background" => cmd_background(&lum, width, height, &opts),
        "detect" => cmd_detect(&lum, width, height, &opts),
        "measure" => cmd_measure(&lum, width, height, green_mask.as_deref(), &opts),
        "fit" => cmd_fit(&lum, width, height, &opts),
        "query" => cmd_query(&lum, width, height, green_mask.as_deref(), &opts),
        "pipeline" => cmd_pipeline(&lum, width, height, channels, green_mask.as_deref(), &opts),
        "dump" => cmd_dump(&lum, width, height, green_mask.as_deref(), &opts),
        other => bail!("Unknown subcommand: {}", other),
    }
}

// ── Background ───────────────────────────────────────────────────────────────

fn cmd_background(lum: &[f32], w: usize, h: usize, opts: &Opts) -> Result<()> {
    let t = Instant::now();
    let cell_size = background::auto_cell_size(w, h);
    let mut bg = background::estimate_background_mesh(lum, w, h, cell_size);
    let mad_noise = bg.noise;
    if opts.mrs_layers > 0 {
        bg.noise = background::estimate_noise_mrs(lum, w, h, opts.mrs_layers.max(1)).max(0.001);
    }
    let elapsed = t.elapsed().as_secs_f64() * 1000.0;

    eprintln!();
    eprintln!("Background estimation ({:.1}ms):", elapsed);
    eprintln!("  Cell size:         {} (auto)", cell_size);
    eprintln!("  Background level:  {:.2}", bg.background);
    eprintln!("  MRS noise:         {:.2} ({}L)", bg.noise, opts.mrs_layers);
    eprintln!("  MAD noise:         {:.2} (ratio {:.2}x)", mad_noise, mad_noise / bg.noise.max(0.001));
    eprintln!("  SNR (bg/noise):    {:.1}", bg.background / bg.noise);

    if let Some(ref bg_map) = bg.background_map {
        let (bmin, bmax) = f32_min_max(bg_map);
        eprintln!("  Background map:    [{:.1}, {:.1}] range", bmin, bmax);
    }
    if let Some(ref noise_map) = bg.noise_map {
        let (nmin, nmax) = f32_min_max(noise_map);
        eprintln!("  Noise map:         [{:.2}, {:.2}] range", nmin, nmax);
    }

    if let Some(ref dir) = opts.save_dir {
        if let Some(ref bg_map) = bg.background_map {
            render::save_grayscale(bg_map, w, h, &dir.join("background_map.png"))?;
            eprintln!("  Saved: background_map.png");
        }
        if let Some(ref noise_map) = bg.noise_map {
            render::save_grayscale(noise_map, w, h, &dir.join("noise_map.png"))?;
            eprintln!("  Saved: noise_map.png");
        }
    }

    Ok(())
}

// ── Detect ───────────────────────────────────────────────────────────────────

fn run_detection(
    lum: &[f32],
    w: usize,
    h: usize,
    opts: &Opts,
) -> Result<(
    background::BackgroundResult,
    Vec<DetectedStar>,
    f32,  // final_fwhm
    f32,  // first-pass measured fwhm
)> {
    let cell_size = background::auto_cell_size(w, h);
    let mut bg = background::estimate_background_mesh(lum, w, h, cell_size);
    if opts.mrs_layers > 0 {
        bg.noise = background::estimate_noise_mrs(lum, w, h, opts.mrs_layers.max(1)).max(0.001);
    }

    let det_params = DetectionParams {
        detection_sigma: opts.sigma,
        min_star_area: 5,
        max_star_area: 2000,
        saturation_limit: 0.95 * 65535.0,
    };

    let bg_map_ref = bg.background_map.as_deref();
    let noise_map_ref = bg.noise_map.as_deref();

    let initial_fwhm = opts.fwhm_override.unwrap_or(3.0);

    // Pass 1
    let pass1 = detection::detect_stars(
        lum, w, h,
        bg.background, bg.noise,
        bg_map_ref, noise_map_ref, &det_params, initial_fwhm, None,
    );

    let measured_fwhm = estimate_fwhm_from_stars(
        lum, w, h, &pass1, bg.background, bg_map_ref,
    );

    let capped = measured_fwhm.min(initial_fwhm * 2.0);

    let (detected, final_fwhm) = if capped > 0.0
        && ((capped - initial_fwhm) / initial_fwhm).abs() > 0.30
    {
        let d = detection::detect_stars(
            lum, w, h,
            bg.background, bg.noise,
            bg_map_ref, noise_map_ref, &det_params, capped, None,
        );
        (d, capped)
    } else {
        (pass1, initial_fwhm)
    };

    Ok((bg, detected, final_fwhm, measured_fwhm))
}

fn cmd_detect(lum: &[f32], w: usize, h: usize, opts: &Opts) -> Result<()> {
    let t = Instant::now();
    let (bg, detected, final_fwhm, measured_fwhm) = run_detection(lum, w, h, opts)?;
    let elapsed = t.elapsed().as_secs_f64() * 1000.0;

    let threshold = opts.sigma * bg.noise;

    eprintln!();
    eprintln!("Detection ({:.1}ms):", elapsed);
    eprintln!("  Background:        {:.2}", bg.background);
    eprintln!("  Noise:             {:.2}", bg.noise);
    eprintln!("  Threshold:         {:.2} ({:.1}sigma)", threshold, opts.sigma);
    eprintln!("  Pass-1 FWHM:       {:.2}px (measured)", measured_fwhm);
    eprintln!("  Final kernel FWHM: {:.2}px", final_fwhm);
    eprintln!("  Stars detected:    {}", detected.len());
    eprintln!();

    // Show top N stars
    let n = detected.len().min(opts.top);
    eprintln!(
        "{:>5} {:>8} {:>8} {:>8} {:>10} {:>6} {:>6} {:>7}",
        "#", "X", "Y", "Peak", "Flux", "Area", "Ecc", "Theta"
    );
    eprintln!("{}", "-".repeat(67));
    for (i, s) in detected.iter().take(n).enumerate() {
        eprintln!(
            "{:>5} {:>8.2} {:>8.2} {:>8.1} {:>10.1} {:>6} {:>6.3} {:>7.3}",
            i + 1, s.x, s.y, s.peak, s.flux, s.area, s.eccentricity, s.theta
        );
    }

    if let Some(ref dir) = opts.save_dir {
        // Save convolution response map
        let sigma = final_fwhm / 2.3548;
        let (kernel_1d, _) = convolution::generate_1d_kernel(sigma);
        let mut conv = vec![0.0_f32; w * h];
        convolution::separable_convolve(lum, w, h, &kernel_1d, &mut conv);
        render::save_grayscale(&conv, w, h, &dir.join("convolution.png"))?;
        eprintln!("\n  Saved: convolution.png");

        // Save detection overlay
        let markers: Vec<(f32, f32, [u8; 3])> = detected
            .iter()
            .map(|s| (s.x, s.y, [0, 255, 0]))
            .collect();
        render::save_with_markers(lum, w, h, &markers, &dir.join("detections.png"))?;
        eprintln!("  Saved: detections.png");
    }

    Ok(())
}

// ── Measure ──────────────────────────────────────────────────────────────────

fn cmd_measure(
    lum: &[f32],
    w: usize,
    h: usize,
    _green_mask: Option<&[bool]>, // unused now but keep for API compat
    opts: &Opts,
) -> Result<()> {
    let t = Instant::now();
    let (bg, detected, _final_fwhm, _) = run_detection(lum, w, h, opts)?;

    let mut measured = metrics::measure_stars(
        lum, w, h, &detected,
        bg.background, bg.background_map.as_deref(),
        None, // no green mask — fit all pixels in interpolated image
        None, // free-beta for debug measure
        CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS,
        None,
    );
    let elapsed = t.elapsed().as_secs_f64() * 1000.0;

    if measured.is_empty() {
        eprintln!("No stars measured.");
        return Ok(());
    }

    // Compute SNR
    let fwhm_vals: Vec<f32> = measured.iter().map(|s| s.fwhm).collect();
    let median_fwhm = sigma_clipped_median(&fwhm_vals);
    snr::compute_star_snr(lum, w, h, &mut measured, median_fwhm);

    // Apply --near filter
    if let Some((nx, ny, nr)) = opts.near {
        let nr_sq = nr * nr;
        let before = measured.len();
        measured.retain(|s| {
            let dx = s.x - nx;
            let dy = s.y - ny;
            dx * dx + dy * dy <= nr_sq
        });
        eprintln!("--near ({:.0},{:.0},r={:.0}): {} / {} stars in range",
            nx, ny, nr, measured.len(), before);
    }

    // Statistics
    let ecc_vals: Vec<f32> = measured.iter().map(|s| s.eccentricity).collect();
    let hfr_vals: Vec<f32> = measured.iter().map(|s| s.hfr).collect();
    let median_ecc = sigma_clipped_median(&ecc_vals);
    let median_hfr = sigma_clipped_median(&hfr_vals);

    eprintln!();
    eprintln!("Measurement ({:.1}ms):", elapsed);
    eprintln!("  Stars measured:    {}", measured.len());
    eprintln!("  Median FWHM:       {:.3}px (sigma-clipped)", median_fwhm);
    eprintln!("  Median eccentricity: {:.3}", median_ecc);
    eprintln!("  Median HFR:        {:.3}px", median_hfr);
    {
        let beta_vals: Vec<f32> = measured.iter().filter_map(|s| s.beta).collect();
        if !beta_vals.is_empty() {
            let median_beta = sigma_clipped_median(&beta_vals);
            let moffat_count = beta_vals.len();
            eprintln!("  Median beta:       {:.2} ({} Moffat fits)", median_beta, moffat_count);
            eprintln!("  Gaussian fallback: {} stars", measured.len() - moffat_count);
        }
    }
    eprintln!();

    // Top N table
    let n = measured.len().min(opts.top);
    eprintln!(
        "{:>5} {:>8} {:>8} {:>7} {:>7} {:>7} {:>6} {:>6} {:>7} {:>6}",
        "#", "X", "Y", "FWHMx", "FWHMy", "FWHM", "Ecc", "HFR", "Theta", "Beta"
    );
    eprintln!("{}", "-".repeat(80));
    for (i, s) in measured.iter().take(n).enumerate() {
        eprintln!(
            "{:>5} {:>8.2} {:>8.2} {:>7.3} {:>7.3} {:>7.3} {:>6.3} {:>6.3} {:>7.3} {:>6}",
            i + 1, s.x, s.y, s.fwhm_x, s.fwhm_y, s.fwhm, s.eccentricity, s.hfr, s.theta,
            s.beta.map_or("-".to_string(), |b| format!("{:.2}", b)),
        );
    }

    if let Some(ref dir) = opts.save_dir {
        let ellipses: Vec<(f32, f32, f32, f32, f32, [u8; 3])> = measured
            .iter()
            .map(|s| {
                let color = if s.beta.is_some() {
                    [0, 255, 0] // green = Moffat
                } else {
                    [255, 255, 0] // yellow = Gaussian fallback
                };
                // semi-axes = fwhm/2
                (s.x, s.y, s.fwhm_x / 2.0, s.fwhm_y / 2.0, s.theta, color)
            })
            .collect();
        render::save_with_ellipses(lum, w, h, &ellipses, &dir.join("measured.png"))?;
        eprintln!("\n  Saved: measured.png");
    }

    Ok(())
}

// ── Fit ──────────────────────────────────────────────────────────────────────

fn cmd_fit(lum: &[f32], w: usize, h: usize, opts: &Opts) -> Result<()> {
    let (sx, sy) = opts.star_pos.context("--star x,y is required for 'fit' subcommand")?;
    let r = opts.stamp_radius;

    let ix = sx.round() as i32;
    let iy = sy.round() as i32;

    if ix - r as i32 <= 0 || iy - r as i32 <= 0
        || ix + r as i32 >= w as i32 - 1 || iy + r as i32 >= h as i32 - 1
    {
        bail!("Star at ({}, {}) with radius {} is too close to the image edge", sx, sy, r);
    }

    // Background estimation for local subtraction
    let cell_size = background::auto_cell_size(w, h);
    let bg = background::estimate_background_mesh(lum, w, h, cell_size);

    let x0 = (ix - r as i32) as usize;
    let y0 = (iy - r as i32) as usize;
    let x1 = (ix + r as i32) as usize;
    let y1 = (iy + r as i32) as usize;
    let stamp_w = x1 - x0 + 1;
    let stamp_h = y1 - y0 + 1;

    // Extract background-subtracted stamp
    let mut stamp = Vec::with_capacity(stamp_w * stamp_h);
    for sy_i in y0..=y1 {
        for sx_i in x0..=x1 {
            let bg_val = bg.background_map.as_ref()
                .map_or(bg.background, |m| m[sy_i * w + sx_i]);
            stamp.push(lum[sy_i * w + sx_i] - bg_val);
        }
    }

    let rel_cx = sx - x0 as f32;
    let rel_cy = sy - y0 as f32;

    // Estimate sigma from halfmax
    let robust_sigma = metrics::estimate_sigma_halfmax(&stamp, stamp_w, rel_cx, rel_cy);

    eprintln!();
    eprintln!("Stamp: {}x{} at ({}, {}), radius={}", stamp_w, stamp_h, sx, sy, r);
    eprintln!("  Robust sigma:  {:.3}px (halfmax estimate)", robust_sigma);

    // Prepare pixel samples for fitting
    let fitting_radius = 5.0_f64.max(4.0 * robust_sigma as f64);
    let fitting_radius_sq = fitting_radius * fitting_radius;
    let cx64 = rel_cx as f64;
    let cy64 = rel_cy as f64;

    let mut pixels = Vec::new();
    for sy_i in 0..stamp_h {
        for sx_i in 0..stamp_w {
            let dx = sx_i as f64 - cx64;
            let dy = sy_i as f64 - cy64;
            if dx * dx + dy * dy <= fitting_radius_sq {
                pixels.push(PixelSample {
                    x: sx_i as f64,
                    y: sy_i as f64,
                    value: stamp[sy_i * stamp_w + sx_i] as f64,
                });
            }
        }
    }

    eprintln!("  Fitting pixels: {} (radius={:.1})", pixels.len(), fitting_radius);

    // Background from annulus
    let bg_inner_sq = fitting_radius_sq;
    let bg_outer = fitting_radius + 3.0;
    let bg_outer_sq = bg_outer * bg_outer;
    let mut annulus_vals: Vec<f64> = Vec::new();
    for sy_i in 0..stamp_h {
        for sx_i in 0..stamp_w {
            let dx = sx_i as f64 - cx64;
            let dy = sy_i as f64 - cy64;
            let r_sq = dx * dx + dy * dy;
            if r_sq > bg_inner_sq && r_sq <= bg_outer_sq {
                annulus_vals.push(stamp[sy_i * stamp_w + sx_i] as f64);
            }
        }
    }
    annulus_vals.sort_by(|a, b| a.total_cmp(b));
    let init_bg = if annulus_vals.is_empty() { 0.0 } else { annulus_vals[annulus_vals.len() / 2] };

    let peak = stamp.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let init_a = peak as f64 - init_bg;

    eprintln!("  Stamp peak:    {:.1}", peak);
    eprintln!("  Annulus bg:    {:.2}", init_bg);
    eprintln!("  Init amplitude: {:.1}", init_a);

    // Gaussian fit
    eprintln!();
    eprintln!("--- Gaussian 2D Fit ---");
    let gauss = fitting::fit_gaussian_2d(
        &pixels, init_bg, init_a,
        cx64, cy64,
        robust_sigma as f64, robust_sigma as f64, 0.0,
        CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS,
    );
    if let Some(ref g) = gauss {
        let fwhm_x = 2.3548 * g.sigma_x;
        let fwhm_y = 2.3548 * g.sigma_y;
        let fwhm = (fwhm_x * fwhm_y).sqrt();
        let ecc = if g.sigma_x >= g.sigma_y {
            (1.0 - (g.sigma_y / g.sigma_x).powi(2)).max(0.0).sqrt()
        } else {
            (1.0 - (g.sigma_x / g.sigma_y).powi(2)).max(0.0).sqrt()
        };
        eprintln!("  Converged:     {}", g.converged);
        eprintln!("  B={:.2}  A={:.1}", g.b, g.a);
        eprintln!("  x0={:.3}  y0={:.3} (stamp-local)", g.x0, g.y0);
        eprintln!("  sigma_x={:.3}  sigma_y={:.3}  theta={:.4} rad", g.sigma_x, g.sigma_y, g.theta);
        eprintln!("  FWHM_x={:.3}  FWHM_y={:.3}  FWHM={:.3}", fwhm_x, fwhm_y, fwhm);
        eprintln!("  Eccentricity:  {:.3}", ecc);
        eprintln!("  Fit residual:  {:.6}", g.fit_residual);
    } else {
        eprintln!("  FAILED");
    }

    // Moffat fit
    eprintln!();
    eprintln!("--- Moffat 2D Fit ---");
    let moffat = fitting::fit_moffat_2d(
        &pixels, init_bg, init_a,
        cx64, cy64,
        robust_sigma as f64, robust_sigma as f64, 0.0,
        CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS,
    );
    if let Some(ref m) = moffat {
        let fwhm_x = m.fwhm_x();
        let fwhm_y = m.fwhm_y();
        let fwhm = (fwhm_x * fwhm_y).sqrt();
        let ecc = if fwhm_x >= fwhm_y {
            (1.0 - (fwhm_y / fwhm_x).powi(2)).max(0.0).sqrt()
        } else {
            (1.0 - (fwhm_x / fwhm_y).powi(2)).max(0.0).sqrt()
        };
        eprintln!("  Converged:     {}", m.converged);
        eprintln!("  B={:.2}  A={:.1}", m.b, m.a);
        eprintln!("  x0={:.3}  y0={:.3} (stamp-local)", m.x0, m.y0);
        eprintln!("  alpha_x={:.3}  alpha_y={:.3}  theta={:.4} rad  beta={:.3}", m.alpha_x, m.alpha_y, m.theta, m.beta);
        eprintln!("  FWHM_x={:.3}  FWHM_y={:.3}  FWHM={:.3}", fwhm_x, fwhm_y, fwhm);
        eprintln!("  Eccentricity:  {:.3}", ecc);
        eprintln!("  Fit residual:  {:.6}", m.fit_residual);
    } else {
        eprintln!("  FAILED");
    }

    // Save debug images
    if let Some(ref dir) = opts.save_dir {
        let label = format!("{}_{}", sx as i32, sy as i32);

        render::save_stamp(&stamp, stamp_w, stamp_h, rel_cx, rel_cy,
            &dir.join(format!("stamp_{}.png", label)))?;
        eprintln!("\n  Saved: stamp_{}.png", label);

        // Build model images for comparison
        if gauss.is_some() || moffat.is_some() {
            let mut model = vec![0.0_f32; stamp_w * stamp_h];
            for y in 0..stamp_h {
                for x in 0..stamp_w {
                    model[y * stamp_w + x] = if let Some(ref m) = moffat {
                        let dx = x as f64 - m.x0;
                        let dy = y as f64 - m.y0;
                        let (ct, st) = (m.theta.cos(), m.theta.sin());
                        let u = dx * ct + dy * st;
                        let v = -dx * st + dy * ct;
                        let q = u * u / (m.alpha_x * m.alpha_x)
                            + v * v / (m.alpha_y * m.alpha_y);
                        (m.b + m.a * (1.0 + q).powf(-m.beta)) as f32
                    } else if let Some(ref g) = gauss {
                        let dx = x as f64 - g.x0;
                        let dy = y as f64 - g.y0;
                        let (ct, st) = (g.theta.cos(), g.theta.sin());
                        let u = dx * ct + dy * st;
                        let v = -dx * st + dy * ct;
                        let q = u * u / (g.sigma_x * g.sigma_x)
                            + v * v / (g.sigma_y * g.sigma_y);
                        (g.b + g.a * (-0.5 * q).exp()) as f32
                    } else {
                        0.0
                    };
                }
            }

            let residual: Vec<f32> = stamp
                .iter()
                .zip(model.iter())
                .map(|(d, m)| d - m)
                .collect();

            render::save_fit_comparison(
                &stamp, &model, &residual,
                stamp_w, stamp_h,
                &dir.join(format!("fit_{}.png", label)),
            )?;
            eprintln!("  Saved: fit_{}.png", label);
        }
    }

    Ok(())
}

// ── Query ────────────────────────────────────────────────────────────────────

fn cmd_query(
    lum: &[f32],
    w: usize,
    h: usize,
    _green_mask: Option<&[bool]>, // unused now but keep for API compat
    opts: &Opts,
) -> Result<()> {
    let (sx, sy) = opts.star_pos.context("--star x,y is required for 'query' subcommand")?;
    let r = opts.stamp_radius;

    // Run full detection + measurement pipeline
    let (bg, detected, _final_fwhm, _) = run_detection(lum, w, h, opts)?;

    let mut measured = metrics::measure_stars(
        lum, w, h, &detected,
        bg.background, bg.background_map.as_deref(),
        None, // no green mask — fit all pixels
        None, // free-beta for debug query
        CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS,
        None,
    );

    if measured.is_empty() {
        bail!("No stars measured in image");
    }

    let fwhm_vals: Vec<f32> = measured.iter().map(|s| s.fwhm).collect();
    let median_fwhm = sigma_clipped_median(&fwhm_vals);
    snr::compute_star_snr(lum, w, h, &mut measured, median_fwhm);

    // Find nearest measured star
    let nearest = measured
        .iter()
        .min_by(|a, b| {
            let da = (a.x - sx) * (a.x - sx) + (a.y - sy) * (a.y - sy);
            let db = (b.x - sx) * (b.x - sx) + (b.y - sy) * (b.y - sy);
            da.total_cmp(&db)
        })
        .unwrap();

    let dist = ((nearest.x - sx).powi(2) + (nearest.y - sy).powi(2)).sqrt();

    eprintln!();
    eprintln!("Query: nearest star to ({:.1}, {:.1}) is at ({:.2}, {:.2}), dist={:.1}px",
        sx, sy, nearest.x, nearest.y, dist);
    eprintln!();

    // Print single-star metrics table
    eprintln!(
        "{:>8} {:>8} {:>7} {:>7} {:>7} {:>6} {:>6} {:>6} {:>7} {:>6}",
        "X", "Y", "FWHMx", "FWHMy", "FWHM", "Ecc", "HFR", "SNR", "Theta", "Beta"
    );
    eprintln!("{}", "-".repeat(80));
    eprintln!(
        "{:>8.2} {:>8.2} {:>7.3} {:>7.3} {:>7.3} {:>6.3} {:>6.3} {:>6.1} {:>7.3} {:>6}",
        nearest.x, nearest.y, nearest.fwhm_x, nearest.fwhm_y, nearest.fwhm,
        nearest.eccentricity, nearest.hfr, nearest.snr, nearest.theta,
        nearest.beta.map_or("-".to_string(), |b| format!("{:.2}", b)),
    );

    // Now run detailed fits at the measured centroid
    let cx = nearest.x;
    let cy = nearest.y;
    let ix = cx.round() as i32;
    let iy = cy.round() as i32;

    if ix - r as i32 <= 0 || iy - r as i32 <= 0
        || ix + r as i32 >= w as i32 - 1 || iy + r as i32 >= h as i32 - 1
    {
        eprintln!("\nStar too close to edge for stamp extraction (radius={})", r);
        return Ok(());
    }

    let x0 = (ix - r as i32) as usize;
    let y0 = (iy - r as i32) as usize;
    let x1 = (ix + r as i32) as usize;
    let y1 = (iy + r as i32) as usize;
    let stamp_w = x1 - x0 + 1;
    let stamp_h = y1 - y0 + 1;

    // Extract background-subtracted stamp
    let mut stamp = Vec::with_capacity(stamp_w * stamp_h);
    for sy_i in y0..=y1 {
        for sx_i in x0..=x1 {
            let bg_val = bg.background_map.as_ref()
                .map_or(bg.background, |m| m[sy_i * w + sx_i]);
            stamp.push(lum[sy_i * w + sx_i] - bg_val);
        }
    }

    let rel_cx = cx - x0 as f32;
    let rel_cy = cy - y0 as f32;
    let robust_sigma = metrics::estimate_sigma_halfmax(&stamp, stamp_w, rel_cx, rel_cy);

    // Prepare fitting pixels
    let fitting_radius = 5.0_f64.max(4.0 * robust_sigma as f64);
    let fitting_radius_sq = fitting_radius * fitting_radius;
    let cx64 = rel_cx as f64;
    let cy64 = rel_cy as f64;

    let mut pixels = Vec::new();
    for sy_i in 0..stamp_h {
        for sx_i in 0..stamp_w {
            let dx = sx_i as f64 - cx64;
            let dy = sy_i as f64 - cy64;
            if dx * dx + dy * dy <= fitting_radius_sq {
                pixels.push(PixelSample {
                    x: sx_i as f64,
                    y: sy_i as f64,
                    value: stamp[sy_i * stamp_w + sx_i] as f64,
                });
            }
        }
    }

    // Background from annulus
    let bg_inner_sq = fitting_radius_sq;
    let bg_outer = fitting_radius + 3.0;
    let bg_outer_sq = bg_outer * bg_outer;
    let mut annulus_vals: Vec<f64> = Vec::new();
    for sy_i in 0..stamp_h {
        for sx_i in 0..stamp_w {
            let dx = sx_i as f64 - cx64;
            let dy = sy_i as f64 - cy64;
            let r_sq = dx * dx + dy * dy;
            if r_sq > bg_inner_sq && r_sq <= bg_outer_sq {
                annulus_vals.push(stamp[sy_i * stamp_w + sx_i] as f64);
            }
        }
    }
    annulus_vals.sort_by(|a, b| a.total_cmp(b));
    let init_bg = if annulus_vals.is_empty() { 0.0 } else { annulus_vals[annulus_vals.len() / 2] };

    let peak = stamp.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let init_a = peak as f64 - init_bg;

    // Gaussian fit
    eprintln!();
    eprintln!("--- Gaussian 2D Fit ---");
    let gauss = fitting::fit_gaussian_2d(
        &pixels, init_bg, init_a,
        cx64, cy64,
        robust_sigma as f64, robust_sigma as f64, 0.0,
        CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS,
    );
    if let Some(ref g) = gauss {
        let fwhm_x = 2.3548 * g.sigma_x;
        let fwhm_y = 2.3548 * g.sigma_y;
        let fwhm = (fwhm_x * fwhm_y).sqrt();
        let ecc = if g.sigma_x >= g.sigma_y {
            (1.0 - (g.sigma_y / g.sigma_x).powi(2)).max(0.0).sqrt()
        } else {
            (1.0 - (g.sigma_x / g.sigma_y).powi(2)).max(0.0).sqrt()
        };
        eprintln!("  Converged:     {}", g.converged);
        eprintln!("  B={:.2}  A={:.1}", g.b, g.a);
        eprintln!("  x0={:.3}  y0={:.3} (stamp-local)", g.x0, g.y0);
        eprintln!("  sigma_x={:.3}  sigma_y={:.3}  theta={:.4} rad", g.sigma_x, g.sigma_y, g.theta);
        eprintln!("  FWHM_x={:.3}  FWHM_y={:.3}  FWHM={:.3}", fwhm_x, fwhm_y, fwhm);
        eprintln!("  Eccentricity:  {:.3}", ecc);
        eprintln!("  Fit residual:  {:.6}", g.fit_residual);
    } else {
        eprintln!("  FAILED");
    }

    // Moffat fit
    eprintln!();
    eprintln!("--- Moffat 2D Fit ---");
    let moffat = fitting::fit_moffat_2d(
        &pixels, init_bg, init_a,
        cx64, cy64,
        robust_sigma as f64, robust_sigma as f64, 0.0,
        CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS,
    );
    if let Some(ref m) = moffat {
        let fwhm_x = m.fwhm_x();
        let fwhm_y = m.fwhm_y();
        let fwhm = (fwhm_x * fwhm_y).sqrt();
        let ecc = if fwhm_x >= fwhm_y {
            (1.0 - (fwhm_y / fwhm_x).powi(2)).max(0.0).sqrt()
        } else {
            (1.0 - (fwhm_x / fwhm_y).powi(2)).max(0.0).sqrt()
        };
        eprintln!("  Converged:     {}", m.converged);
        eprintln!("  B={:.2}  A={:.1}", m.b, m.a);
        eprintln!("  x0={:.3}  y0={:.3} (stamp-local)", m.x0, m.y0);
        eprintln!("  alpha_x={:.3}  alpha_y={:.3}  theta={:.4} rad  beta={:.3}", m.alpha_x, m.alpha_y, m.theta, m.beta);
        eprintln!("  FWHM_x={:.3}  FWHM_y={:.3}  FWHM={:.3}", fwhm_x, fwhm_y, fwhm);
        eprintln!("  Eccentricity:  {:.3}", ecc);
        eprintln!("  Fit residual:  {:.6}", m.fit_residual);
    } else {
        eprintln!("  FAILED");
    }

    // Save debug images
    if let Some(ref dir) = opts.save_dir {
        let label = format!("query_{}_{}", cx as i32, cy as i32);

        render::save_stamp(&stamp, stamp_w, stamp_h, rel_cx, rel_cy,
            &dir.join(format!("stamp_{}.png", label)))?;
        eprintln!("\n  Saved: stamp_{}.png", label);

        if gauss.is_some() || moffat.is_some() {
            let mut model = vec![0.0_f32; stamp_w * stamp_h];
            for y in 0..stamp_h {
                for x in 0..stamp_w {
                    model[y * stamp_w + x] = if let Some(ref m) = moffat {
                        let dx = x as f64 - m.x0;
                        let dy = y as f64 - m.y0;
                        let (ct, st) = (m.theta.cos(), m.theta.sin());
                        let u = dx * ct + dy * st;
                        let v = -dx * st + dy * ct;
                        let q = u * u / (m.alpha_x * m.alpha_x)
                            + v * v / (m.alpha_y * m.alpha_y);
                        (m.b + m.a * (1.0 + q).powf(-m.beta)) as f32
                    } else if let Some(ref g) = gauss {
                        let dx = x as f64 - g.x0;
                        let dy = y as f64 - g.y0;
                        let (ct, st) = (g.theta.cos(), g.theta.sin());
                        let u = dx * ct + dy * st;
                        let v = -dx * st + dy * ct;
                        let q = u * u / (g.sigma_x * g.sigma_x)
                            + v * v / (g.sigma_y * g.sigma_y);
                        (g.b + g.a * (-0.5 * q).exp()) as f32
                    } else {
                        0.0
                    };
                }
            }

            let residual: Vec<f32> = stamp
                .iter()
                .zip(model.iter())
                .map(|(d, m)| d - m)
                .collect();

            render::save_fit_comparison(
                &stamp, &model, &residual,
                stamp_w, stamp_h,
                &dir.join(format!("fit_{}.png", label)),
            )?;
            eprintln!("  Saved: fit_{}.png", label);
        }
    }

    Ok(())
}

// ── Pipeline ─────────────────────────────────────────────────────────────────

fn cmd_pipeline(
    lum: &[f32],
    w: usize,
    h: usize,
    channels: usize,
    _green_mask: Option<&[bool]>,  // unused now but keep for API compat
    opts: &Opts,
) -> Result<()> {
    let mut analyzer = ImageAnalyzer::new()
        .with_detection_sigma(opts.sigma)
        .with_measure_cap(opts.measure_cap)
        .with_mrs_layers(opts.mrs_layers)
        .with_max_stars(opts.top)
        .with_fit_max_iter(opts.fit_max_iter)
        .with_fit_tolerance(opts.fit_tolerance)
        .with_fit_max_rejects(opts.fit_max_rejects);

    if opts.no_debayer {
        analyzer = analyzer.without_debayer();
    }

    let result = analyzer.analyze_data(lum, w, h, channels)?;
    let t = &result.stage_timing;

    eprintln!("\n=== Stage 1: Background Estimation ===");
    eprintln!("  Time: {:.1}ms", t.background_ms);
    eprintln!("  Background: {:.2}  Noise: {:.2}", result.background, result.noise);

    eprintln!("\n=== Stage 2: Star Detection ===");
    eprintln!("  Time: {:.1}ms", t.detection_pass1_ms + t.calibration_ms + t.detection_pass2_ms);
    eprintln!("  Pass 1: {} detections", result.pass1_detections);
    eprintln!("  Calibration FWHM: {:.3}px", result.calibrated_fwhm);
    eprintln!("  Final: {} detections (FWHM={:.2})", result.stars_detected, result.measured_fwhm_kernel);
    if result.possibly_trailed {
        eprintln!("  Trail R²: {:.4}  Rayleigh trailed: true", result.trail_r_squared);
    } else {
        eprintln!("  Trail R²: {:.4}  Rayleigh trailed: false", result.trail_r_squared);
    }

    eprintln!("\n=== Stage 3: PSF Measurement ===");
    eprintln!("  Time: {:.1}ms", t.measurement_ms);
    eprintln!("  Measured: {} / {} detected ({} Moffat / {} Gauss)",
        result.stars_measured, result.stars_detected,
        result.moffat_count, result.gaussian_count);

    eprintln!("\n=== Stage 4: Statistics ===");
    eprintln!("  Time: {:.1}ms", t.statistics_ms + t.snr_ms);
    if result.possibly_trailed {
        eprintln!("  Trailed: true (Rayleigh angle coherence)");
    }
    eprintln!("  Median FWHM:     {:.3}px", result.median_fwhm);
    eprintln!("  Median ecc:      {:.3}", result.median_eccentricity);
    eprintln!("  Median HFR:      {:.3}px", result.median_hfr);
    eprintln!("  SNR weight:      {:.2}", result.snr_weight);
    eprintln!("  PSF signal:      {:.2}", result.psf_signal);
    eprintln!("  Frame SNR:       {:.1}", result.frame_snr);
    if let Some(beta) = result.median_beta {
        eprintln!("  Median beta:     {:.2} ({} Moffat / {} Gauss)",
            beta, result.moffat_count, result.gaussian_count);
    }
    eprintln!("  Image: {}x{}, {} ch", result.width, result.height, result.source_channels);

    eprintln!("\n=== Total pipeline: {:.1}ms ===\n", t.total_ms);

    // Per-star table header
    if !result.stars.is_empty() {
        eprintln!("    #        X        Y   FWHMx   FWHMy    FWHM    Ecc    HFR    SNR   Theta   Beta");
        eprintln!("----------------------------------------------------------------------------------------");
        for (i, s) in result.stars.iter().enumerate() {
            let beta_str = s.beta.map_or("     -".to_string(), |b| format!("{:6.2}", b));
            eprintln!("{:5} {:8.2} {:8.2} {:7.3} {:7.3} {:7.3} {:6.3} {:6.3} {:6.1} {:8.4} {}",
                i + 1, s.x, s.y, s.fwhm_x, s.fwhm_y, s.fwhm,
                s.eccentricity, s.hfr, s.snr, s.theta, beta_str);
        }
    }

    Ok(())
}

// ── Dump ─────────────────────────────────────────────────────────────────────

fn cmd_dump(
    lum: &[f32],
    w: usize,
    h: usize,
    _green_mask: Option<&[bool]>, // unused now but keep for API compat
    opts: &Opts,
) -> Result<()> {
    let (bg, detected, _final_fwhm, _) = run_detection(lum, w, h, opts)?;

    let mut measured = metrics::measure_stars(
        lum, w, h, &detected,
        bg.background, bg.background_map.as_deref(),
        None, // no green mask — fit all pixels
        None, // free-beta for debug dump
        CALIBRATION_MAX_ITER, CALIBRATION_CONV_TOL, CALIBRATION_MAX_REJECTS,
        None,
    );

    if measured.is_empty() {
        eprintln!("No stars measured.");
        return Ok(());
    }

    let fwhm_vals: Vec<f32> = measured.iter().map(|s| s.fwhm).collect();
    let median_fwhm = sigma_clipped_median(&fwhm_vals);
    snr::compute_star_snr(lum, w, h, &mut measured, median_fwhm);

    // Apply --near filter
    if let Some((nx, ny, nr)) = opts.near {
        let nr_sq = nr * nr;
        let before = measured.len();
        measured.retain(|s| {
            let dx = s.x - nx;
            let dy = s.y - ny;
            dx * dx + dy * dy <= nr_sq
        });
        eprintln!("--near ({:.0},{:.0},r={:.0}): {} / {} stars in range",
            nx, ny, nr, measured.len(), before);
    }

    // TSV header to stdout
    println!("x\ty\tpeak\tflux\tfwhm_x\tfwhm_y\tfwhm\tecc\thfr\tsnr\ttheta\tbeta\tfit_res");
    for s in &measured {
        println!(
            "{:.3}\t{:.3}\t{:.1}\t{:.1}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.1}\t{:.4}\t{}\t{:.4}",
            s.x, s.y, s.peak, s.flux,
            s.fwhm_x, s.fwhm_y, s.fwhm, s.eccentricity, s.hfr, s.snr, s.theta,
            s.beta.map_or("-".to_string(), |b| format!("{:.3}", b)),
            s.fit_residual,
        );
    }

    // Percentile distribution summary to stderr
    let count = measured.len();
    let mut fwhm_s: Vec<f32> = measured.iter().map(|s| s.fwhm).collect();
    let mut ecc_s: Vec<f32> = measured.iter().map(|s| s.eccentricity).collect();
    let mut hfr_s: Vec<f32> = measured.iter().map(|s| s.hfr).collect();
    let mut snr_s: Vec<f32> = measured.iter().map(|s| s.snr).collect();
    fwhm_s.sort_by(|a, b| a.total_cmp(b));
    ecc_s.sort_by(|a, b| a.total_cmp(b));
    hfr_s.sort_by(|a, b| a.total_cmp(b));
    snr_s.sort_by(|a, b| a.total_cmp(b));

    let pcts = [0.10_f32, 0.25, 0.50, 0.75, 0.90];
    eprintln!();
    eprintln!("Distribution ({} stars):", count);
    eprintln!("         {:>7} {:>7} {:>7} {:>7} {:>7}", "P10", "P25", "P50", "P75", "P90");
    eprintln!("  FWHM  {:>7.3} {:>7.3} {:>7.3} {:>7.3} {:>7.3}",
        percentile(&fwhm_s, pcts[0]), percentile(&fwhm_s, pcts[1]),
        percentile(&fwhm_s, pcts[2]), percentile(&fwhm_s, pcts[3]),
        percentile(&fwhm_s, pcts[4]));
    eprintln!("  Ecc   {:>7.3} {:>7.3} {:>7.3} {:>7.3} {:>7.3}",
        percentile(&ecc_s, pcts[0]), percentile(&ecc_s, pcts[1]),
        percentile(&ecc_s, pcts[2]), percentile(&ecc_s, pcts[3]),
        percentile(&ecc_s, pcts[4]));
    eprintln!("  HFR   {:>7.3} {:>7.3} {:>7.3} {:>7.3} {:>7.3}",
        percentile(&hfr_s, pcts[0]), percentile(&hfr_s, pcts[1]),
        percentile(&hfr_s, pcts[2]), percentile(&hfr_s, pcts[3]),
        percentile(&hfr_s, pcts[4]));
    eprintln!("  SNR   {:>7.1} {:>7.1} {:>7.1} {:>7.1} {:>7.1}",
        percentile(&snr_s, pcts[0]), percentile(&snr_s, pcts[1]),
        percentile(&snr_s, pcts[2]), percentile(&snr_s, pcts[3]),
        percentile(&snr_s, pcts[4]));

    let beta_vals: Vec<f32> = measured.iter().filter_map(|s| s.beta).collect();
    if !beta_vals.is_empty() {
        let mut beta_s = beta_vals;
        beta_s.sort_by(|a, b| a.total_cmp(b));
        eprintln!("  Beta  {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2}",
            percentile(&beta_s, pcts[0]), percentile(&beta_s, pcts[1]),
            percentile(&beta_s, pcts[2]), percentile(&beta_s, pcts[3]),
            percentile(&beta_s, pcts[4]));
    }

    let mut fitres_s: Vec<f32> = measured.iter().map(|s| s.fit_residual).collect();
    fitres_s.sort_by(|a, b| a.total_cmp(b));
    eprintln!("  FitR  {:>7.4} {:>7.4} {:>7.4} {:>7.4} {:>7.4}",
        percentile(&fitres_s, pcts[0]), percentile(&fitres_s, pcts[1]),
        percentile(&fitres_s, pcts[2]), percentile(&fitres_s, pcts[3]),
        percentile(&fitres_s, pcts[4]));

    eprintln!();
    eprintln!("{} stars written to stdout", count);
    Ok(())
}

// ── Compare ─────────────────────────────────────────────────────────────────

/// Reference data from one row of a PixInsight SubframeSelector CSV.
struct PiRef {
    file: String,
    fwhm: f64,
    eccentricity: f64,
    noise_adu: f64,
    stars: usize,
    snr: f64,
    psf_signal_weight: f64,
}

/// Our measurements for one file.
struct OurResult {
    median_fwhm: f64,
    #[allow(dead_code)]
    mean_fwhm: f64,
    median_ecc: f64,
    mean_ecc: f64,
    noise: f64,
    stars: usize,
    frame_snr: f64,
    snr_weight: f64,
    median_snr: f64,
}

/// Parse a PixInsight SubframeSelector CSV.
/// Returns (settings header lines, per-file reference data).
fn parse_pi_csv(path: &std::path::Path) -> Result<(Vec<String>, Vec<PiRef>)> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read PI CSV: {:?}", path))?;

    let lines: Vec<&str> = content.lines().collect();
    // Find the header row (starts with "Index,")
    let header_idx = lines.iter().position(|l| l.starts_with("Index,"))
        .context("No header row starting with 'Index,' found in CSV")?;

    let settings: Vec<String> = lines[..header_idx].iter().map(|s| s.to_string()).collect();

    // Parse header to find column indices
    let headers: Vec<&str> = lines[header_idx].split(',').collect();
    let col = |name: &str| -> Result<usize> {
        headers.iter().position(|h| *h == name)
            .with_context(|| format!("Missing column '{}' in PI CSV", name))
    };
    let i_file = col("File")?;
    let i_fwhm = col("FWHM")?;
    let i_ecc = col("Eccentricity")?;
    let i_noise = col("Noise")?;
    let i_stars = col("Stars")?;
    let i_snr = col("SNR")?;
    let i_psfw = col("PSF Signal Weight")?;
    let max_col = i_file.max(i_fwhm).max(i_ecc).max(i_noise).max(i_stars)
        .max(i_snr).max(i_psfw);

    let mut refs = Vec::new();
    for line in &lines[(header_idx + 1)..] {
        if line.trim().is_empty() { continue; }
        // CSV with quoted file paths — split carefully
        let fields = csv_split(line);
        if fields.len() <= max_col {
            continue;
        }
        let file = fields[i_file].trim_matches('"').to_string();
        let fwhm: f64 = fields[i_fwhm].parse().unwrap_or(0.0);
        let ecc: f64 = fields[i_ecc].parse().unwrap_or(0.0);
        let noise_norm: f64 = fields[i_noise].parse().unwrap_or(0.0);
        let stars: usize = fields[i_stars].parse().unwrap_or(0);
        let snr: f64 = fields[i_snr].parse().unwrap_or(0.0);
        let psf_signal_weight: f64 = fields[i_psfw].parse().unwrap_or(0.0);

        refs.push(PiRef {
            file,
            fwhm,
            eccentricity: ecc,
            noise_adu: noise_norm * 65535.0,
            stars,
            snr,
            psf_signal_weight,
        });
    }

    Ok((settings, refs))
}

/// Split a CSV line respecting quoted fields.
fn csv_split(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;
    for ch in line.chars() {
        match ch {
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                fields.push(current.clone());
                current.clear();
            }
            _ => current.push(ch),
        }
    }
    fields.push(current);
    fields
}

/// Run analysis on a single file and return our measurements.
fn analyze_one_file(
    path: &std::path::Path,
    opts: &Opts,
) -> Result<OurResult> {
    let mut analyzer = ImageAnalyzer::new()
        .with_detection_sigma(opts.sigma)
        .with_measure_cap(opts.measure_cap)
        .with_mrs_layers(opts.mrs_layers)
        .with_fit_max_iter(opts.fit_max_iter)
        .with_fit_tolerance(opts.fit_tolerance)
        .with_fit_max_rejects(opts.fit_max_rejects);

    if opts.no_debayer {
        analyzer = analyzer.without_debayer();
    }

    let result = analyzer.analyze(path)?;

    // Compute mean values from per-star data for compare output
    let mean_fwhm = if result.stars.is_empty() { 0.0 } else {
        result.stars.iter().map(|s| s.fwhm as f64).sum::<f64>() / result.stars.len() as f64
    };
    let mean_ecc = if result.stars.is_empty() { 0.0 } else {
        result.stars.iter().map(|s| s.eccentricity as f64).sum::<f64>() / result.stars.len() as f64
    };

    Ok(OurResult {
        median_fwhm: result.median_fwhm as f64,
        mean_fwhm,
        median_ecc: result.median_eccentricity as f64,
        mean_ecc,
        noise: result.noise as f64,
        stars: result.stars_measured,
        frame_snr: result.frame_snr as f64,
        snr_weight: result.snr_weight as f64,
        median_snr: result.median_snr as f64,
    })
}

fn cmd_compare(opts: &Opts) -> Result<()> {
    let total = Instant::now();

    // Parse PI CSV
    let (settings, pi_refs) = parse_pi_csv(opts.input.as_ref())?;
    eprintln!("PixInsight SubframeSelector reference: {} files", pi_refs.len());
    for s in &settings {
        if !s.is_empty() {
            eprintln!("  {}", s);
        }
    }

    // Show our settings
    eprintln!();
    eprintln!("Our settings:");
    eprintln!("  --sigma {:.1}  --mrs {}", opts.sigma, opts.mrs_layers);
    eprintln!();

    // Per-file comparison
    struct FileResult {
        fname: String,
        pi: PiRef,
        ours: OurResult,
        fwhm_pct: f64,
        ecc_diff: f64,
        noise_ratio: f64,
        star_ratio: f64,
        snr_ratio: f64,
    }

    let n_files = pi_refs.len();
    let mut results: Vec<FileResult> = Vec::with_capacity(n_files);
    let mut errors: Vec<(String, String)> = Vec::new();

    for (i, pi) in pi_refs.into_iter().enumerate() {
        let path = std::path::Path::new(&pi.file);
        let fname = path.file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| pi.file.clone());

        eprint!("\r  [{}/{}] {}", i + 1, n_files, fname);

        match analyze_one_file(path, opts) {
            Ok(ours) => {
                let fwhm_pct = if pi.fwhm > 0.0 {
                    (ours.median_fwhm - pi.fwhm) / pi.fwhm * 100.0
                } else { 0.0 };
                let ecc_diff = ours.median_ecc - pi.eccentricity;
                let noise_ratio = if pi.noise_adu > 0.0 {
                    ours.noise / pi.noise_adu
                } else { 0.0 };
                let star_ratio = if pi.stars > 0 {
                    ours.stars as f64 / pi.stars as f64
                } else { 0.0 };
                let snr_ratio = if pi.snr > 0.0 {
                    ours.frame_snr / pi.snr
                } else { 0.0 };

                results.push(FileResult {
                    fname,
                    pi,
                    ours,
                    fwhm_pct,
                    ecc_diff,
                    noise_ratio,
                    star_ratio,
                    snr_ratio,
                });
            }
            Err(e) => {
                errors.push((fname, format!("{:#}", e)));
            }
        }
    }
    eprintln!("\r  Processed {} files in {:.1}s{}",
        results.len(),
        total.elapsed().as_secs_f64(),
        if errors.is_empty() { String::new() }
        else { format!(" ({} errors)", errors.len()) },
    );

    if results.is_empty() {
        bail!("No files successfully processed");
    }

    // ── Per-file table ──────────────────────────────────────────────────────
    eprintln!();
    eprintln!("{:<52} {:>7} {:>7} {:>7} {:>7} {:>8} {:>7} {:>8} {:>8}",
        "File", "PI_FWHM", "FWHM", "PI_Ecc", "Ecc", "Noise_r", "Star_r", "PI_SNR", "SNR");
    eprintln!("{}", "-".repeat(120));

    let mut sorted = results.iter().collect::<Vec<_>>();
    sorted.sort_by(|a, b| a.fname.cmp(&b.fname));
    for r in &sorted {
        eprintln!("{:<52} {:>7.3} {:>7.3} {:>7.4} {:>7.4} {:>7.3}x {:>6.2}x {:>8.2} {:>8.2}",
            &r.fname[..r.fname.len().min(52)],
            r.pi.fwhm, r.ours.median_fwhm,
            r.pi.eccentricity, r.ours.median_ecc,
            r.noise_ratio, r.star_ratio,
            r.pi.snr, r.ours.frame_snr);
    }

    // ── Aggregate statistics ────────────────────────────────────────────────
    let fwhm_pcts: Vec<f64> = results.iter().map(|r| r.fwhm_pct).collect();
    let ecc_diffs: Vec<f64> = results.iter().map(|r| r.ecc_diff).collect();
    let noise_rats: Vec<f64> = results.iter().map(|r| r.noise_ratio).collect();
    let star_rats: Vec<f64> = results.iter().map(|r| r.star_ratio).collect();
    let snr_rats: Vec<f64> = results.iter().map(|r| r.snr_ratio).collect();

    // Also compute mean-ecc based offsets for comparison
    let mean_ecc_diffs: Vec<f64> = results.iter()
        .map(|r| r.ours.mean_ecc - r.pi.eccentricity).collect();

    fn mean(v: &[f64]) -> f64 { v.iter().sum::<f64>() / v.len() as f64 }
    fn median_f64(v: &[f64]) -> f64 {
        let mut s = v.to_vec();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        s[s.len() / 2]
    }
    fn std_dev(v: &[f64]) -> f64 {
        let m = mean(v);
        (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
    }
    fn percentile_f64(v: &[f64], p: f64) -> f64 {
        let mut s = v.to_vec();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let idx = ((p / 100.0) * (s.len() - 1) as f64).round() as usize;
        s[idx.min(s.len() - 1)]
    }

    eprintln!();
    eprintln!("=== AGGREGATE ({} files) ===", results.len());
    eprintln!("{:<22} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "", "FWHM %", "Ecc off", "Noise ×", "Stars ×", "SNR ×");
    eprintln!("{}", "-".repeat(76));
    eprintln!("{:<22} {:>+10.3}% {:>+10.4} {:>10.3}x {:>10.2}x {:>10.3}x",
        "Mean", mean(&fwhm_pcts), mean(&ecc_diffs), mean(&noise_rats), mean(&star_rats), mean(&snr_rats));
    eprintln!("{:<22} {:>+10.3}% {:>+10.4} {:>10.3}x {:>10.2}x {:>10.3}x",
        "Median", median_f64(&fwhm_pcts), median_f64(&ecc_diffs),
        median_f64(&noise_rats), median_f64(&star_rats), median_f64(&snr_rats));
    eprintln!("{:<22} {:>10.3}% {:>10.4} {:>10.3}x {:>10.2}x {:>10.3}x",
        "Std dev", std_dev(&fwhm_pcts), std_dev(&ecc_diffs),
        std_dev(&noise_rats), std_dev(&star_rats), std_dev(&snr_rats));

    // Ecc with mean vs median
    eprintln!();
    eprintln!("Eccentricity (ours − PI):");
    eprintln!("  Using our median: mean={:+.4}, median={:+.4}",
        mean(&ecc_diffs), median_f64(&ecc_diffs));
    eprintln!("  Using our mean:   mean={:+.4}, median={:+.4}",
        mean(&mean_ecc_diffs), median_f64(&mean_ecc_diffs));

    // ── Linear regression: PI_ecc = a × our_ecc + b ────────────────────────
    fn linear_regression(xs: &[f64], ys: &[f64]) -> (f64, f64, f64) {
        let mx = mean(xs);
        let my = mean(ys);
        let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
        let sxy: f64 = xs.iter().zip(ys.iter())
            .map(|(x, y)| (x - mx) * (y - my)).sum();
        let slope = sxy / sxx;
        let intercept = my - slope * mx;
        let ss_res: f64 = xs.iter().zip(ys.iter())
            .map(|(x, y)| (y - (slope * x + intercept)).powi(2)).sum();
        let ss_tot: f64 = ys.iter().map(|y| (y - my).powi(2)).sum();
        let r_sq = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { 0.0 };
        (slope, intercept, r_sq)
    }

    let our_eccs: Vec<f64> = results.iter().map(|r| r.ours.median_ecc).collect();
    let pi_eccs: Vec<f64> = results.iter().map(|r| r.pi.eccentricity).collect();
    let (ecc_slope, ecc_intercept, ecc_r_sq) = linear_regression(&our_eccs, &pi_eccs);

    eprintln!();
    eprintln!("Eccentricity regression: PI = {:.4} × ours + {:.4}  (R²={:.4})",
        ecc_slope, ecc_intercept, ecc_r_sq);

    // Ecc percentiles
    eprintln!();
    eprintln!("Eccentricity offset percentiles:");
    eprintln!("  P10={:+.4}  P25={:+.4}  P50={:+.4}  P75={:+.4}  P90={:+.4}",
        percentile_f64(&ecc_diffs, 10.0), percentile_f64(&ecc_diffs, 25.0),
        percentile_f64(&ecc_diffs, 50.0), percentile_f64(&ecc_diffs, 75.0),
        percentile_f64(&ecc_diffs, 90.0));

    // ── SNR analysis ────────────────────────────────────────────────────────
    let our_snrs: Vec<f64> = results.iter().map(|r| r.ours.frame_snr).collect();
    let pi_snrs: Vec<f64> = results.iter().map(|r| r.pi.snr).collect();
    let (snr_slope, snr_intercept, snr_r_sq) = linear_regression(&our_snrs, &pi_snrs);

    eprintln!();
    eprintln!("=== FRAME SNR ANALYSIS ===");
    eprintln!("  Our frame_snr = bg / noise;  PI SNR = their image SNR metric");
    eprintln!();
    eprintln!("  SNR ratio (ours/PI):  mean={:.3}x  median={:.3}x  std={:.3}x",
        mean(&snr_rats), median_f64(&snr_rats), std_dev(&snr_rats));
    eprintln!("  Regression: PI_SNR = {:.4} × ours + {:.4}  (R²={:.4})",
        snr_slope, snr_intercept, snr_r_sq);
    eprintln!("  Percentiles:");
    eprintln!("    P10={:.3}x  P25={:.3}x  P50={:.3}x  P75={:.3}x  P90={:.3}x",
        percentile_f64(&snr_rats, 10.0), percentile_f64(&snr_rats, 25.0),
        percentile_f64(&snr_rats, 50.0), percentile_f64(&snr_rats, 75.0),
        percentile_f64(&snr_rats, 90.0));

    // SNR Weight comparison (our snr_weight vs PI PSF Signal Weight)
    let our_snrw: Vec<f64> = results.iter().map(|r| r.ours.snr_weight).collect();
    let pi_psfw: Vec<f64> = results.iter().map(|r| r.pi.psf_signal_weight).collect();
    let snrw_rats: Vec<f64> = our_snrw.iter().zip(pi_psfw.iter())
        .map(|(o, p)| if *p > 0.0 { o / p } else { 0.0 }).collect();
    let (snrw_slope, snrw_intercept, snrw_r_sq) = linear_regression(&our_snrw, &pi_psfw);

    eprintln!();
    eprintln!("=== SNR WEIGHT vs PI PSF SIGNAL WEIGHT ===");
    eprintln!("  Our: median(flux)²/(noise²×bg);  PI: PSF Signal Weight");
    eprintln!();
    eprintln!("  Weight ratio (ours/PI):  mean={:.3}x  median={:.3}x  std={:.3}x",
        mean(&snrw_rats), median_f64(&snrw_rats), std_dev(&snrw_rats));
    eprintln!("  Regression: PI_PSFw = {:.6} × ours + {:.6}  (R²={:.4})",
        snrw_slope, snrw_intercept, snrw_r_sq);

    // Median per-star SNR summary
    let our_msnrs: Vec<f64> = results.iter().map(|r| r.ours.median_snr).collect();
    eprintln!();
    eprintln!("=== MEDIAN PER-STAR SNR ===");
    eprintln!("  mean={:.1}  median={:.1}  min={:.1}  max={:.1}",
        mean(&our_msnrs), median_f64(&our_msnrs),
        our_msnrs.iter().cloned().fold(f64::INFINITY, f64::min),
        our_msnrs.iter().cloned().fold(f64::NEG_INFINITY, f64::max));

    // ── Subgroup analysis (by exposure time pattern) ────────────────────────
    let groups: Vec<(&str, Box<dyn Fn(&FileResult) -> bool>)> = vec![
        ("300s exposures", Box::new(|r: &FileResult| r.fname.contains("300.00s"))),
        ("30s exposures", Box::new(|r: &FileResult| r.fname.contains("30.00s"))),
    ];

    for (label, filter) in &groups {
        let sub: Vec<&FileResult> = results.iter().filter(|r| filter(r)).collect();
        if sub.is_empty() { continue; }
        let fp: Vec<f64> = sub.iter().map(|r| r.fwhm_pct).collect();
        let ed: Vec<f64> = sub.iter().map(|r| r.ecc_diff).collect();
        let nr: Vec<f64> = sub.iter().map(|r| r.noise_ratio).collect();
        let sr: Vec<f64> = sub.iter().map(|r| r.star_ratio).collect();
        let snr_r: Vec<f64> = sub.iter().map(|r| r.snr_ratio).collect();
        let ms: Vec<f64> = sub.iter().map(|r| r.ours.median_snr).collect();
        eprintln!();
        eprintln!("--- {} (n={}) ---", label, sub.len());
        eprintln!("  FWHM:  mean={:+.3}%  median={:+.3}%  std={:.3}%",
            mean(&fp), median_f64(&fp), std_dev(&fp));
        eprintln!("  Ecc:   mean={:+.4}  median={:+.4}  std={:.4}",
            mean(&ed), median_f64(&ed), std_dev(&ed));
        eprintln!("  Noise: mean={:.3}x  median={:.3}x",
            mean(&nr), median_f64(&nr));
        eprintln!("  Stars: mean={:.2}x  median={:.2}x",
            mean(&sr), median_f64(&sr));
        eprintln!("  SNR:   mean={:.3}x  median={:.3}x  (frame_snr ratio)",
            mean(&snr_r), median_f64(&snr_r));
        eprintln!("  Star SNR: mean={:.1}  median={:.1}  (our median per-star SNR)",
            mean(&ms), median_f64(&ms));
    }

    // ── Worst outliers ──────────────────────────────────────────────────────
    let mut by_fwhm = results.iter().collect::<Vec<_>>();
    by_fwhm.sort_by(|a, b| b.fwhm_pct.abs().partial_cmp(&a.fwhm_pct.abs()).unwrap());
    eprintln!();
    eprintln!("=== WORST FWHM OUTLIERS ===");
    eprintln!("{:<48} {:>7} {:>7} {:>7} {:>7} {:>7}",
        "File", "FWHM%", "PI", "Ours", "PI_Ecc", "Ecc");
    for r in by_fwhm.iter().take(5) {
        eprintln!("{:<48} {:>+6.2}% {:>7.3} {:>7.3} {:>7.4} {:>7.4}",
            &r.fname[..r.fname.len().min(48)],
            r.fwhm_pct, r.pi.fwhm, r.ours.median_fwhm,
            r.pi.eccentricity, r.ours.median_ecc);
    }

    let mut by_ecc = results.iter().collect::<Vec<_>>();
    by_ecc.sort_by(|a, b| b.ecc_diff.abs().partial_cmp(&a.ecc_diff.abs()).unwrap());
    eprintln!();
    eprintln!("=== WORST ECCENTRICITY OUTLIERS ===");
    eprintln!("{:<48} {:>7} {:>7} {:>7} {:>7} {:>7}",
        "File", "Ecc_d", "PI", "Ours", "FWHM%", "Stars×");
    for r in by_ecc.iter().take(5) {
        eprintln!("{:<48} {:>+.4} {:>7.4} {:>7.4} {:>+6.2}% {:>6.2}x",
            &r.fname[..r.fname.len().min(48)],
            r.ecc_diff, r.pi.eccentricity, r.ours.median_ecc,
            r.fwhm_pct, r.star_ratio);
    }

    let mut by_snr = results.iter().collect::<Vec<_>>();
    by_snr.sort_by(|a, b| (b.snr_ratio - 1.0).abs().partial_cmp(&(a.snr_ratio - 1.0).abs()).unwrap());
    eprintln!();
    eprintln!("=== WORST FRAME SNR OUTLIERS ===");
    eprintln!("{:<48} {:>8} {:>8} {:>7} {:>7} {:>7}",
        "File", "PI_SNR", "OurSNR", "Ratio", "FWHM%", "Noise×");
    for r in by_snr.iter().take(5) {
        eprintln!("{:<48} {:>8.2} {:>8.2} {:>6.3}x {:>+6.2}% {:>6.3}x",
            &r.fname[..r.fname.len().min(48)],
            r.pi.snr, r.ours.frame_snr, r.snr_ratio,
            r.fwhm_pct, r.noise_ratio);
    }

    // Errors
    if !errors.is_empty() {
        eprintln!();
        eprintln!("=== ERRORS ({}) ===", errors.len());
        for (fname, err) in &errors {
            eprintln!("  {}: {}", fname, err);
        }
    }

    // TSV to stdout for downstream analysis
    println!("file\tpi_fwhm\tours_fwhm\tfwhm_pct\tpi_ecc\tours_ecc\tours_mean_ecc\tecc_diff\tpi_noise\tours_noise\tnoise_ratio\tpi_stars\tours_stars\tstar_ratio\tpi_snr\tours_frame_snr\tsnr_ratio\tours_snr_weight\tpi_psf_signal_weight\tours_median_snr");
    for r in &results {
        println!("{}\t{:.4}\t{:.4}\t{:+.3}\t{:.4}\t{:.4}\t{:.4}\t{:+.4}\t{:.2}\t{:.2}\t{:.4}\t{}\t{}\t{:.3}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.1}",
            r.fname, r.pi.fwhm, r.ours.median_fwhm, r.fwhm_pct,
            r.pi.eccentricity, r.ours.median_ecc, r.ours.mean_ecc, r.ecc_diff,
            r.pi.noise_adu, r.ours.noise, r.noise_ratio,
            r.pi.stars, r.ours.stars, r.star_ratio,
            r.pi.snr, r.ours.frame_snr, r.snr_ratio,
            r.ours.snr_weight, r.pi.psf_signal_weight, r.ours.median_snr);
    }

    Ok(())
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn percentile(sorted: &[f32], p: f32) -> f32 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = p * (sorted.len() - 1) as f32;
    let lo = idx.floor() as usize;
    let hi = (lo + 1).min(sorted.len() - 1);
    let frac = idx - lo as f32;
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

fn f32_min_max(data: &[f32]) -> (f32, f32) {
    let mut min = f32::INFINITY;
    let mut max = f32::NEG_INFINITY;
    for &v in data {
        if v.is_finite() {
            if v < min { min = v; }
            if v > max { max = v; }
        }
    }
    (min, max)
}
