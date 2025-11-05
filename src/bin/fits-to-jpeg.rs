use anyhow::{Context, Result};
use fits_converter::FitsConverter;
use std::env;
use std::process;

fn print_usage(program: &str) {
    eprintln!("FITS to JPEG Converter");
    eprintln!();
    eprintln!("Usage: {} <input.fits> <output.jpg> [OPTIONS]", program);
    eprintln!();
    eprintln!("Options:");
    eprintln!("  --downscale <N>      Downscale factor (1 = no downscaling, default: 1)");
    eprintln!("  --quality <Q>        JPEG quality 1-100 (default: 95)");
    eprintln!("  --no-debayer         Disable Bayer pattern debayering");
    eprintln!("  --log                Show detailed conversion information");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  {} image.fits image.jpg", program);
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
        println!("Converting FITS to JPEG...");
        println!("  Input:  {}", input_path);
        println!("  Output: {}", output_path);
        println!("  Downscale: {}x", downscale);
        println!("  Quality: {}", quality);
        println!("  Debayer: {}", if apply_debayer { "enabled" } else { "disabled" });
    }

    // Build converter with options
    let mut converter = FitsConverter::new()
        .with_downscale(downscale)
        .with_quality(quality);
      //  .with_logging(log_enabled);

    if !apply_debayer {
        converter = converter.without_debayer();
    }

    // Perform conversion
    converter.convert(input_path, output_path)
        .context("Failed to convert FITS to JPEG")?;

    if log_enabled {
        println!("Conversion successful!");
    }

    Ok(())
}
