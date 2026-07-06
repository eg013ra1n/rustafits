//! Real-data integration gate for `ImageAnalyzer::detect_fast`.
//!
//! Compares the fast detector against the precise `ImageAnalyzer::analyze`
//! pipeline on the real FITS/XISF files in `tests/`. The coverage question it
//! answers: does the fast pass-1 path (blind plate solving's star source)
//! still find the real astronomical sources the precise pipeline trusts?
//!
//! Method: fast runs with max_stars=2000 (the cap ranking itself is pinned by
//! `adaptive_detection::tests::cap_keeps_brightest_not_first_in_scan_order`),
//! slow with max_stars=500; for each of the brightest 100 slow stars we ask
//! whether a fast star sits within 2 px.
//!
//! Measured baselines (2026-07-06, after the full-scan-then-cap fix):
//!   - mono.fits  100/100 (nearest p90 0.21 px)
//!   - osc.fits   100/100 (p90 0.21 px)
//!   - test.xisf   92/100 (p90 1.66 px)
//!   - cocoon.fits 22/100 — KNOWN DEFICIENCY: on strong-nebula frames the
//!     hfd_at flux-weighted centroid drags 3–8 px on the nebula pedestal
//!     (its shrink-to-symmetry box sees the whole window as signal), and the
//!     flux ranking mixes saturated stars the slow path gates out. The 15%
//!     floor is a regression tripwire: the scan-order stripe truncation this
//!     gate caught scored 1–3%.
//!
//! Also asserts detection-count parity and that fast is not slower than slow.
//!
//! **Must be run in release mode** — debug builds make the matched-filter
//! convolution slow enough to dominate everything:
//!     cargo test -p rustafits --test fast_detect_real --release -- --nocapture

use std::path::Path;

use astroimage::{ImageAnalyzer, StarMetrics};

/// (file, min slow-top-100 coverage at 2 px)
const TEST_FILES: &[(&str, f32)] = &[
    ("tests/cocoon.fits", 0.15),
    ("tests/mono.fits", 0.85),
    ("tests/osc.fits", 0.85),
    ("tests/test.xisf", 0.85),
];

fn has(name: &str) -> bool {
    Path::new(name).exists()
}

/// Fraction of the brightest `top_n` slow stars that have at least one fast
/// star within `tol_px`. Measures real-star coverage — independent of
/// spurious extras the fast path may pick up from extended sources.
fn slow_coverage_by_fast(
    slow: &[StarMetrics],
    fast: &[astroimage::FastStar],
    top_n: usize,
    tol_px: f32,
) -> f32 {
    if slow.is_empty() {
        return 0.0;
    }
    let n = slow.len().min(top_n);
    let tol2 = tol_px * tol_px;
    let mut hits = 0usize;
    for s in slow.iter().take(n) {
        let mut best = f32::INFINITY;
        for f in fast {
            let dx = f.x - s.x;
            let dy = f.y - s.y;
            let d2 = dx * dx + dy * dy;
            if d2 < best {
                best = d2;
            }
        }
        if best < tol2 {
            hits += 1;
        }
    }
    hits as f32 / n as f32
}

struct Measurement {
    file: &'static str,
    fast_count: usize,
    slow_count: usize,
    fast_ms: f64,
    slow_ms: f64,
    matched_frac: f32,
}

fn run_one(file: &'static str, min_coverage: f32) -> Option<Measurement> {
    if !has(file) {
        eprintln!("SKIP: {} not found", file);
        return None;
    }

    // Production-shaped fast run (solver uses max_stars=600, registration 500):
    // count and speed gates measure THIS configuration.
    let fast_analyzer = ImageAnalyzer::new()
        .with_detection_sigma(5.0)
        .with_max_stars(600);
    let fast_res = fast_analyzer
        .detect_fast(file)
        .unwrap_or_else(|e| panic!("detect_fast failed on {}: {}", file, e));

    // Deep fast run — generous cap so the coverage gate measures detection
    // (can the fast path find the stars at all), not cap ranking.
    let deep_res = ImageAnalyzer::new()
        .with_detection_sigma(5.0)
        .with_max_stars(2000)
        .detect_fast(file)
        .unwrap_or_else(|e| panic!("detect_fast (deep) failed on {}: {}", file, e));

    // Slow path (reference)
    let slow_analyzer = ImageAnalyzer::new()
        .with_detection_sigma(5.0)
        .with_max_stars(500);
    let slow_res = slow_analyzer
        .analyze(file)
        .unwrap_or_else(|e| panic!("analyze failed on {}: {}", file, e));

    let fast_count = fast_res.stars.len();
    let slow_count = slow_res.stars.len();
    let fast_ms = fast_res.timing.total_ms;
    let slow_ms = slow_res.stage_timing.total_ms;
    let matched_frac = slow_coverage_by_fast(&slow_res.stars, &deep_res.stars, 100, 2.0);

    eprintln!(
        "{file}: fast={fast_count}★ {fast_ms:.0}ms  slow={slow_count}★ {slow_ms:.0}ms  ratio={:.2}×  slow_top100_covered_2px={:.1}% (gate {:.0}%)",
        slow_ms / fast_ms.max(0.001),
        matched_frac * 100.0,
        min_coverage * 100.0,
    );
    eprintln!(
        "  fast stages: prep={:.0} bg={:.0} detect={:.0} read={:.0}",
        fast_res.timing.prep_ms,
        fast_res.timing.background_ms,
        fast_res.timing.detection_ms,
        fast_res.timing.read_ms,
    );

    assert!(
        fast_count >= 10,
        "{}: expected ≥10 fast stars, got {}",
        file,
        fast_count
    );
    assert!(
        slow_count >= 10,
        "{}: expected ≥10 slow stars, got {}",
        file,
        slow_count
    );
    assert!(
        fast_count as f32 >= 0.8 * slow_count as f32,
        "{}: fast coverage too low: fast={} slow={} ({:.0}% < 80%)",
        file,
        fast_count,
        slow_count,
        100.0 * fast_count as f32 / slow_count as f32,
    );
    assert!(
        matched_frac >= min_coverage,
        "{}: slow-top-100 coverage regressed: {:.1}% < {:.0}%",
        file,
        matched_frac * 100.0,
        min_coverage * 100.0,
    );
    assert!(
        fast_ms <= slow_ms + 50.0,
        "{}: fast ({:.0} ms) is materially slower than slow ({:.0} ms) — regression",
        file,
        fast_ms,
        slow_ms,
    );

    Some(Measurement {
        file,
        fast_count,
        slow_count,
        fast_ms,
        slow_ms,
        matched_frac,
    })
}

#[test]
fn fast_detect_matches_analyze_on_real_data() {
    let mut measurements = Vec::new();
    for (file, min_coverage) in TEST_FILES {
        if let Some(m) = run_one(file, *min_coverage) {
            measurements.push(m);
        }
    }

    if measurements.is_empty() {
        eprintln!("SKIP: no real test files available");
        return;
    }

    let total_fast: f64 = measurements.iter().map(|m| m.fast_ms).sum();
    let total_slow: f64 = measurements.iter().map(|m| m.slow_ms).sum();
    let ratio = total_slow / total_fast.max(0.001);
    eprintln!(
        "\nAggregate: {} files, fast={:.0}ms  slow={:.0}ms  speedup={:.2}×",
        measurements.len(),
        total_fast,
        total_slow,
        ratio,
    );

    // Sanity: at least one file matches ≥ 90% — the best case should be
    // excellent on a frame where pass-1 and pass-2 detections align.
    let best = measurements
        .iter()
        .map(|m| m.matched_frac)
        .fold(0.0_f32, f32::max);
    assert!(
        best >= 0.90,
        "best-case correctness too low: {:.1}%",
        best * 100.0,
    );

    for m in &measurements {
        let _ = (m.file, m.fast_count, m.slow_count);
    }
}
