use astroimage::ImageAnalyzer;
use std::path::Path;

fn has_test_file(name: &str) -> bool {
    Path::new(&format!("tests/{name}")).exists()
}

#[test]
fn analyze_mono_fits() {
    if !has_test_file("mono.fits") {
        eprintln!("Skipping: tests/mono.fits not found");
        return;
    }

    let result = ImageAnalyzer::new()
        .with_detection_sigma(5.0)
        .with_max_stars(200)
        .analyze("tests/mono.fits")
        .unwrap();

    assert!(result.width > 0);
    assert!(result.height > 0);
    assert_eq!(result.source_channels, 1);
    assert!(result.background > 0.0, "bg: {}", result.background);
    assert!(result.noise > 0.0, "noise: {}", result.noise);
    assert!(
        result.stars_detected > 0,
        "Should detect at least 1 star in mono.fits"
    );
    assert!(!result.stars.is_empty());

    // FWHM should be near PixInsight's ~2.0 px (Moffat fit default)
    assert!(
        result.median_fwhm > 1.2 && result.median_fwhm < 3.5,
        "median FWHM {} out of range [1.2, 3.5]",
        result.median_fwhm
    );

    // Eccentricity: 0 to 1
    assert!(
        result.median_eccentricity >= 0.0 && result.median_eccentricity <= 1.0,
        "median eccentricity {} out of range",
        result.median_eccentricity
    );

    // SNR should be positive
    assert!(result.median_snr > 0.0, "median SNR: {}", result.median_snr);

    // HFR should be positive and correlate with FWHM
    assert!(result.median_hfr > 0.0, "median HFR: {}", result.median_hfr);

    // SNR weight and PSF signal should be positive
    assert!(result.snr_weight > 0.0, "SNR weight: {}", result.snr_weight);
    assert!(result.psf_signal > 0.0, "PSF signal: {}", result.psf_signal);

    eprintln!(
        "mono.fits: {}x{}, {} stars, FWHM={:.2}, ecc={:.3}, SNR={:.1}, HFR={:.2}, SNRw={:.1}, PSF={:.1}",
        result.width,
        result.height,
        result.stars_detected,
        result.median_fwhm,
        result.median_eccentricity,
        result.median_snr,
        result.median_hfr,
        result.snr_weight,
        result.psf_signal,
    );
}

#[test]
fn analyze_osc_fits() {
    if !has_test_file("osc.fits") {
        eprintln!("Skipping: tests/osc.fits not found");
        return;
    }

    let result = ImageAnalyzer::new()
        .with_detection_sigma(5.0)
        .with_max_stars(200)
        .analyze("tests/osc.fits")
        .unwrap();

    // Green interpolation keeps native resolution, mono green channel
    assert_eq!(
        result.source_channels, 1,
        "OSC green interpolation should produce 1 channel"
    );

    assert!(result.width > 0);
    assert!(result.height > 0);
    assert!(result.background > 0.0);
    assert!(result.noise > 0.0);
    assert!(
        result.stars_detected > 0,
        "Should detect at least 1 star in osc.fits"
    );

    // FWHM at native resolution, green-pixel-only fitting, Moffat default (~2.5 PI ref)
    assert!(
        result.median_fwhm > 1.5 && result.median_fwhm < 3.5,
        "median FWHM {} out of range [1.5, 3.5]",
        result.median_fwhm
    );

    eprintln!(
        "osc.fits: {}x{}, {} stars, FWHM={:.2}, ecc={:.3}, SNR={:.1}, HFR={:.2}",
        result.width,
        result.height,
        result.stars_detected,
        result.median_fwhm,
        result.median_eccentricity,
        result.median_snr,
        result.median_hfr,
    );
}

#[test]
fn analyze_xisf() {
    if !has_test_file("test.xisf") {
        eprintln!("Skipping: tests/test.xisf not found");
        return;
    }

    let result = ImageAnalyzer::new()
        .with_max_stars(100)
        .analyze("tests/test.xisf")
        .unwrap();

    assert!(result.width > 0);
    assert!(result.height > 0);
    assert!(result.background > 0.0 || result.background == 0.0); // XISF might have near-zero bg
    assert!(result.noise > 0.0);

    eprintln!(
        "test.xisf: {}x{}, {} stars, FWHM={:.2}, bg={:.1}, noise={:.3}",
        result.width,
        result.height,
        result.stars_detected,
        result.median_fwhm,
        result.background,
        result.noise,
    );
}

#[test]
fn analyze_mono_fit_method_populated() {
    if !has_test_file("mono.fits") {
        eprintln!("Skipping: tests/mono.fits not found");
        return;
    }

    let result = ImageAnalyzer::new()
        .with_max_stars(50)
        .analyze("tests/mono.fits")
        .unwrap();

    assert!(!result.stars.is_empty());
    assert!(result.median_fwhm > 0.5 && result.median_fwhm < 50.0);

    // Every star should have a valid fit method
    for star in &result.stars {
        assert!(
            star.fit_method == astroimage::FitMethod::FreeMoffat
            || star.fit_method == astroimage::FitMethod::FixedMoffat
            || star.fit_method == astroimage::FitMethod::Gaussian
            || star.fit_method == astroimage::FitMethod::Moments,
        );
    }

    // Most stars should use Moffat (FixedMoffat from pass 2)
    let moffat_count = result.stars.iter()
        .filter(|s| matches!(s.fit_method, astroimage::FitMethod::FreeMoffat | astroimage::FitMethod::FixedMoffat))
        .count();
    eprintln!(
        "Fit methods: {}/{} Moffat, FWHM={:.2}",
        moffat_count, result.stars.len(), result.median_fwhm
    );
}

#[test]
fn analyze_osc_without_debayer() {
    if !has_test_file("osc.fits") {
        eprintln!("Skipping: tests/osc.fits not found");
        return;
    }

    let result = ImageAnalyzer::new()
        .without_debayer()
        .with_max_stars(100)
        .analyze("tests/osc.fits")
        .unwrap();

    // Without debayer, raw CFA treated as mono at native resolution
    assert_eq!(
        result.source_channels, 1,
        "without_debayer should produce 1 channel"
    );

    // Native resolution preserved (no halving)
    assert!(result.width > 0);
    assert!(result.height > 0);

    eprintln!(
        "osc.fits (no debayer): {}x{}, {} stars, FWHM={:.2}",
        result.width, result.height, result.stars_detected, result.median_fwhm,
    );
}

#[test]
fn analyze_with_thread_pool() {
    if !has_test_file("mono.fits") {
        eprintln!("Skipping: tests/mono.fits not found");
        return;
    }

    let pool = std::sync::Arc::new(
        astroimage::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap(),
    );

    let result = ImageAnalyzer::new()
        .with_thread_pool(pool)
        .with_max_stars(50)
        .analyze("tests/mono.fits")
        .unwrap();

    assert!(result.stars_detected > 0);
    assert!(result.median_fwhm > 0.0);

    eprintln!(
        "Thread pool: {} stars, FWHM={:.2}",
        result.stars_detected, result.median_fwhm
    );
}

#[test]
fn analyze_auto_mesh_background() {
    if !has_test_file("mono.fits") {
        eprintln!("Skipping: tests/mono.fits not found");
        return;
    }

    // Mesh background is now always-on (auto-tuned cell size)
    let result = ImageAnalyzer::new()
        .with_max_stars(100)
        .analyze("tests/mono.fits")
        .unwrap();

    assert!(result.background > 0.0);
    assert!(result.noise > 0.0);
    assert!(result.stars_detected > 0);

    eprintln!(
        "Auto mesh background: {} stars, bg={:.1}, noise={:.3}, FWHM={:.2}",
        result.stars_detected, result.background, result.noise, result.median_fwhm
    );
}

#[test]
fn analyze_data_synthetic() {
    // Generate synthetic star field and analyze via analyze_data()
    let width = 200;
    let height = 200;
    let background = 1000.0_f32;
    let mut data = vec![background; width * height];

    // Add 5 Gaussian stars
    let stars = [
        (50.0_f32, 50.0_f32, 5000.0_f32, 3.0_f32),
        (100.0, 80.0, 3000.0, 3.0),
        (150.0, 120.0, 7000.0, 3.0),
        (30.0, 160.0, 4000.0, 3.0),
        (170.0, 40.0, 6000.0, 3.0),
    ];

    for &(sx, sy, amp, sigma) in &stars {
        let r = (4.0 * sigma) as i32;
        let inv_2s2 = 1.0 / (2.0 * sigma * sigma);
        for dy in -r..=r {
            for dx in -r..=r {
                let px = (sx as i32 + dx) as usize;
                let py = (sy as i32 + dy) as usize;
                if px < width && py < height {
                    data[py * width + px] +=
                        amp * (-inv_2s2 * (dx as f32 * dx as f32 + dy as f32 * dy as f32)).exp();
                }
            }
        }
    }

    // Add noise
    let mut rng = 42u64;
    for val in data.iter_mut() {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
        let noise = ((rng >> 33) as f32 / (1u64 << 31) as f32 - 0.5) * 100.0;
        *val += noise;
    }

    let result = ImageAnalyzer::new()
        .with_detection_sigma(5.0)
        .analyze_data(&data, width, height, 1)
        .unwrap();

    assert!(
        result.stars_detected >= 3,
        "Expected at least 3 synthetic stars, got {}",
        result.stars_detected
    );

    let expected_fwhm = 2.3548 * 3.0; // ~7.06
    // Moffat/Gaussian fit may under- or over-estimate depending on convergence
    assert!(
        result.median_fwhm > 3.0 && result.median_fwhm < 15.0,
        "FWHM {:.2} out of reasonable range for synthetic stars with sigma=3",
        result.median_fwhm,
    );

    eprintln!(
        "Synthetic: {} stars, FWHM={:.2} (expected {:.2}), ecc={:.3}, SNR={:.1}",
        result.stars_detected, result.median_fwhm, expected_fwhm, result.median_eccentricity, result.median_snr
    );
}

#[test]
fn analyze_cocoon_fits() {
    if !has_test_file("cocoon.fits") {
        eprintln!("Skipping: tests/cocoon.fits not found");
        return;
    }

    let result = ImageAnalyzer::new()
        .with_detection_sigma(5.0)
        .with_max_stars(500)
        .analyze("tests/cocoon.fits")
        .unwrap();

    assert!(result.width > 0);
    assert!(result.height > 0);
    assert!(result.background > 0.0, "bg: {}", result.background);
    assert!(result.noise > 0.0, "noise: {}", result.noise);

    // Dense star field should detect many stars (no cap on measurement)
    assert!(
        result.stars_detected >= 500,
        "Expected at least 500 stars in cocoon.fits (crowded field), got {}",
        result.stars_detected
    );

    // FWHM should be reasonable (1.5-5.0 range)
    assert!(
        result.median_fwhm > 1.5 && result.median_fwhm < 5.0,
        "median FWHM {} out of range [1.5, 5.0]",
        result.median_fwhm
    );

    eprintln!(
        "cocoon.fits: {}x{}, {} stars, FWHM={:.2}, ecc={:.3}, SNR={:.1}, kernel_fwhm={:.2}",
        result.width,
        result.height,
        result.stars_detected,
        result.median_fwhm,
        result.median_eccentricity,
        result.median_snr,
        result.measured_fwhm_kernel,
    );
}

#[test]
fn analyze_mono_moffat() {
    if !has_test_file("mono.fits") {
        eprintln!("Skipping: tests/mono.fits not found");
        return;
    }

    // Moffat is now always-on (unified pipeline)
    let result = ImageAnalyzer::new()
        .with_max_stars(100)
        .analyze("tests/mono.fits")
        .unwrap();

    assert!(!result.stars.is_empty());

    // Moffat should report beta values
    assert!(
        result.median_beta.is_some(),
        "Moffat fit should produce median_beta"
    );
    let beta = result.median_beta.unwrap();
    assert!(
        beta > 1.0 && beta < 10.0,
        "median beta {} out of reasonable range [1.0, 10.0]",
        beta
    );

    // Count how many stars got Moffat fits (have beta)
    let moffat_count = result.stars.iter().filter(|s| s.beta.is_some()).count();
    let total = result.stars.len();

    eprintln!(
        "Moffat: {} stars, {}/{} with Moffat fit, beta={:.2}, FWHM={:.2}",
        result.stars_detected,
        moffat_count,
        total,
        beta,
        result.median_fwhm,
    );
}

#[test]
fn analyze_with_mrs_layers() {
    if !has_test_file("mono.fits") {
        eprintln!("Skipping: tests/mono.fits not found");
        return;
    }

    let result = ImageAnalyzer::new()
        .with_mrs_layers(2)
        .with_max_stars(100)
        .analyze("tests/mono.fits")
        .unwrap();

    // Should detect stars with MRS noise
    assert!(result.stars_detected > 0);
    assert!(result.noise > 0.0);

    eprintln!(
        "MRS layers 2: {} stars, bg={:.1}, noise={:.3}, FWHM={:.2}",
        result.stars_detected, result.background, result.noise, result.median_fwhm,
    );
}
