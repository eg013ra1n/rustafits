use anyhow::{Context, Result};
use astroimage::{annotate_image, AnnotationConfig, ImageAnalyzer, ImageConverter};
use std::env;
use std::process;

fn print_usage(program: &str) {
    eprintln!("FITS/XISF to JPEG Converter");
    eprintln!();
    eprintln!("Usage: {} <input> <output.jpg> [OPTIONS]", program);
    eprintln!();
    eprintln!("Supported input formats:");
    eprintln!("  .fits, .fit   - FITS (Flexible Image Transport System)");
    eprintln!("  .xisf         - XISF (PixInsight native format)");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  --downscale <N>      Downscale factor (1 = no downscaling, default: 1)");
    eprintln!("  --quality <Q>        JPEG quality 1-100 (default: 95)");
    eprintln!("  --no-debayer         Disable Bayer pattern debayering");
    eprintln!("  --preview            Enable preview mode (2x2 binning for mono, faster)");
    eprintln!("  --annotate           Overlay star detection ellipses on the output image");
    eprintln!("  --max-stars <N>      Max stars for annotation analysis (default: 200)");
    eprintln!("  --log                Show detailed conversion information");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  {} image.fits image.jpg", program);
    eprintln!("  {} image.xisf image.jpg", program);
    eprintln!("  {} image.fits image.jpg --downscale 2 --quality 90", program);
    eprintln!("  {} raw.fits output.jpg --no-debayer --log", program);
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}

fn run() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        print_usage(&args[0]);
        return Err(anyhow::anyhow!("Missing required arguments"));
    }

    let input_path = &args[1];
    let output_path = &args[2];

    // Parse optional arguments
    let mut downscale = 1;
    let mut quality = 95;
    let mut apply_debayer = true;
    let mut preview_mode = false;
    let mut annotate = false;
    let mut max_stars: usize = 200;
    let mut log_enabled = false;

    let mut i = 3;
    while i < args.len() {
        match args[i].as_str() {
            "--downscale" => {
                if i + 1 >= args.len() {
                    return Err(anyhow::anyhow!("--downscale requires a value"));
                }
                downscale = args[i + 1]
                    .parse::<usize>()
                    .context("Invalid downscale factor")?;
                if downscale == 0 {
                    return Err(anyhow::anyhow!("Downscale factor must be >= 1"));
                }
                i += 2;
            }
            "--quality" => {
                if i + 1 >= args.len() {
                    return Err(anyhow::anyhow!("--quality requires a value"));
                }
                quality = args[i + 1]
                    .parse::<u8>()
                    .context("Invalid quality value")?;
                if quality < 1 || quality > 100 {
                    return Err(anyhow::anyhow!("Quality must be between 1 and 100"));
                }
                i += 2;
            }
            "--no-debayer" => {
                apply_debayer = false;
                i += 1;
            }
            "--preview" => {
                preview_mode = true;
                i += 1;
            }
            "--annotate" => {
                annotate = true;
                i += 1;
            }
            "--max-stars" => {
                if i + 1 >= args.len() {
                    return Err(anyhow::anyhow!("--max-stars requires a value"));
                }
                max_stars = args[i + 1]
                    .parse::<usize>()
                    .context("Invalid max-stars value")?;
                if max_stars == 0 {
                    return Err(anyhow::anyhow!("--max-stars must be >= 1"));
                }
                i += 2;
            }
            "--log" => {
                log_enabled = true;
                i += 1;
            }
            "--help" | "-h" => {
                print_usage(&args[0]);
                return Ok(());
            }
            _ => {
                return Err(anyhow::anyhow!("Unknown option: {}", args[i]));
            }
        }
    }

    if log_enabled {
        println!("Converting to JPEG...");
        println!("  Input:  {}", input_path);
        println!("  Output: {}", output_path);
        println!("  Downscale: {}x", downscale);
        println!("  Quality: {}", quality);
        println!("  Debayer: {}", if apply_debayer { "enabled" } else { "disabled" });
        println!("  Preview mode: {}", if preview_mode { "enabled" } else { "disabled" });
    }

    // Build converter with options
    let mut converter = ImageConverter::new()
        .with_downscale(downscale)
        .with_quality(quality);

    if !apply_debayer {
        converter = converter.without_debayer();
    }

    if preview_mode {
        converter = converter.with_preview_mode();
    }

    if annotate {
        // Process image (get mutable ProcessedImage for annotation)
        let mut image = converter.process(input_path)
            .context("Image processing failed")?;

        // Run analysis on the same input
        let result = ImageAnalyzer::new()
            .with_max_stars(max_stars)
            .analyze(input_path)
            .context("Analysis failed")?;

        if log_enabled {
            println!("Analysis: {} stars detected, median FWHM={:.2}, median ecc={:.3}",
                result.stars.len(), result.median_fwhm, result.median_eccentricity);
        }

        // Annotate
        annotate_image(&mut image, &result, &AnnotationConfig::default());

        // Save
        astroimage::ImageConverter::save_processed(&image, output_path, quality)
            .context("Image save failed")?;
    } else {
        // Standard conversion
        converter.convert(input_path, output_path)
            .context("Conversion failed")?;
    }

    if log_enabled {
        println!("Conversion successful!");
    }

    Ok(())
}
