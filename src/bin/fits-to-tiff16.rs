use anyhow::{Context, Result};
use fits_converter::{fits::FitsImage, processor::{process_fits, ProcessorConfig, ProcessedImage}, output::save_tiff_16bit};
use std::env;
use std::process;

fn print_usage(program: &str) {
    eprintln!("FITS to 16-bit TIFF Converter (debayered, unstretched)");
    eprintln!();
    eprintln!("Usage: {} <input.fits> <output.tiff>", program);
    eprintln!();
    eprintln!("This outputs the debayered RGB data in 16-bit TIFF format");
    eprintln!("without any stretching applied. Use this to verify debayering");
    eprintln!("and apply your own stretch in PixInsight or other tools.");
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

    println!("Converting FITS to 16-bit TIFF...");
    println!("  Input:  {}", input_path);
    println!("  Output: {}", output_path);

    // Load FITS
    let fits = FitsImage::open(input_path)
        .context("Failed to read FITS file")?;

    println!("  Dimensions: {}x{}", fits.metadata.width, fits.metadata.height);
    if let Some(pattern) = &fits.metadata.bayer_pattern {
        println!("  Bayer pattern: {:?}", pattern);
    }

    // Process with skip_stretch enabled
    let config = ProcessorConfig {
        downscale_factor: 1,
        manual_stretch: None,
        apply_debayer: true,
        skip_stretch: true,  // Don't stretch, return raw 16-bit data
    };

    let result = process_fits(fits, config)
        .context("Failed to process FITS file")?;

    // Save as 16-bit TIFF
    match result {
        ProcessedImage::RGB16(ref rgb_data) => {
            save_tiff_16bit(rgb_data, output_path)
                .context("Failed to save TIFF file")?;
        }
        _ => {
            return Err(anyhow::anyhow!("Expected RGB16 output for TIFF export"));
        }
    }

    println!("Conversion successful!");
    println!("Output is debayered 16-bit RGB TIFF, ready for stretching in PixInsight.");

    Ok(())
}
