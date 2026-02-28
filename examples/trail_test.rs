use astroimage::ImageAnalyzer;

fn analyze_file(label: &str, path: &str) {
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("{label}");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let analyzer = ImageAnalyzer::new();
    let r = analyzer.analyze(path).unwrap();

    println!("  Stars detected:       {}", r.stars_detected);
    println!("  Stars after filter:   {}", r.stars.len());
    println!("  Median FWHM:          {:.3}", r.median_fwhm);
    println!("  Median eccentricity:  {:.3}", r.median_eccentricity);
    println!("  Median SNR:           {:.1}", r.median_snr);
    println!("  PSF signal:           {:.1}", r.psf_signal);
    println!("  Trail R²:             {:.4}", r.trail_r_squared);
    println!("  Possibly trailed:     {}", r.possibly_trailed);

    // Show theta distribution
    let n = r.stars.len();
    if n >= 5 {
        let mut thetas: Vec<f32> = r.stars.iter().map(|s| s.theta.to_degrees()).collect();
        thetas.sort_by(|a, b| a.total_cmp(b));

        let (sum_cos, sum_sin) = r.stars.iter().fold((0.0f64, 0.0f64), |(sc, ss), s| {
            let a = 2.0 * s.theta as f64;
            (sc + a.cos(), ss + a.sin())
        });
        let mean_theta = (sum_sin.atan2(sum_cos) * 0.5).to_degrees();

        println!("\n  Theta ({} measured stars):", n);
        println!("    Mean theta:   {mean_theta:.1}deg");
        println!(
            "    min={:.1}  p25={:.1}  median={:.1}  p75={:.1}  max={:.1}",
            thetas[0],
            thetas[n / 4],
            thetas[n / 2],
            thetas[3 * n / 4],
            thetas[n - 1]
        );
    }

    if !r.stars.is_empty() {
        println!("\n  Top 5 stars:");
        for (i, s) in r.stars.iter().take(5).enumerate() {
            println!(
                "    #{}: ecc={:.3} theta={:.1}deg fwhm={:.2} peak={:.0}",
                i + 1,
                s.eccentricity,
                s.theta.to_degrees(),
                s.fwhm,
                s.peak
            );
        }
    }
    println!();
}

fn main() {
    let good = "/Volumes/BigMac/Users/astrobureau/Pictures/Unsorted/edph/NINA/COPIED/CAMERA_DUO/Barnard 150/2024-09-28/LIGHT/2024-09-28_23-16-30_-10.00_180.00s_0281.fits";
    let s_trail = "/Volumes/BigMac/Users/astrobureau/Pictures/Unsorted/edph/NINA/COPIED/CAMERA_DUO/Barnard 150/2024-09-28/LIGHT/2024-09-29_05-41-37_-9.90_180.00s_0381.fits";
    let linear_trail = "/Volumes/BigMac/Users/astrobureau/Pictures/Unsorted/edph/NINA/COPIED/CAMERA_DUO/Barnard 150/2024-09-27/LIGHT/2024-09-28_05-41-28_-10.00_180.00s_0236.fits";

    analyze_file("GOOD STARS (0281)", good);
    analyze_file("S-SHAPED TRAILS (0381)", s_trail);
    analyze_file("LINEAR TRAILS (0236)", linear_trail);
}
