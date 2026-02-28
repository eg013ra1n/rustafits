use astroimage::ImageAnalyzer;

fn main() {
    let dir = "tests/no_trails";
    let mut files: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map_or(false, |ext| ext == "fits"))
        .map(|e| e.path())
        .collect();
    files.sort();

    println!(
        "{:<12} {:>5} {:>5} {:>6} {:>6} {:>6} {:>6} {:>7} {:>7} {:>8} {:>6} {:>6} {:>6} {:>7}",
        "FILE", "DET", "KEPT", "FWHM", "ECC", "SNR", "HFR", "SNR_dB", "SNR_W", "PSF_SIG", "W", "H", "R²", "TRAIL?"
    );
    println!("{}", "-".repeat(125));

    for path in &files {
        let name = path.file_name().unwrap().to_str().unwrap();
        // Extract short name (frame number)
        let short = if let Some(pos) = name.rfind('_') {
            &name[pos + 1..name.len() - 5] // strip .fits
        } else {
            &name[..name.len().min(12)]
        };

        let analyzer = ImageAnalyzer::new();
        match analyzer.analyze(path) {
            Ok(r) => {
                println!(
                    "{:<12} {:>5} {:>5} {:>6.2} {:>6.3} {:>6.1} {:>6.2} {:>7.2} {:>7.3} {:>8.1} {:>6} {:>6} {:>6.3} {:>7}",
                    short,
                    r.stars_detected,
                    r.stars.len(),
                    r.median_fwhm,
                    r.median_eccentricity,
                    r.median_snr,
                    r.median_hfr,
                    r.snr_db,
                    r.snr_weight,
                    r.psf_signal,
                    r.width,
                    r.height,
                    r.trail_r_squared,
                    if r.possibly_trailed { "YES" } else { "no" },
                );
            }
            Err(e) => {
                println!("{:<12} ERROR: {}", short, e);
            }
        }
    }
}
