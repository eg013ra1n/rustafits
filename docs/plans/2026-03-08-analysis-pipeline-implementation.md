# Analysis Pipeline Redesign — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Unify the analysis pipeline into one codepath with Moffat-primary PSF fitting, auto-tuned mesh background, MRS wavelet noise, and iterative PSF subtraction for deblending.

**Architecture:** Two-pass detection — pass 1 discovers the PSF model from bright calibration stars (free-beta Moffat), pass 2 applies fixed-beta Moffat to all stars with iterative subtraction for blends. Background is always mesh-grid with auto-tuned cell size and MRS noise. Config knobs control thresholds only, not algorithms.

**Tech Stack:** Rust, rayon (parallelism), f64 fitting arithmetic, SIMD (AVX2/NEON) for convolution

**Design doc:** `docs/plans/2026-03-08-analysis-pipeline-redesign.md`

---

### Task 1: Add FitMethod Enum

**Files:**
- Modify: `src/analysis/mod.rs:50-77` (StarMetrics struct)

**Step 1: Add FitMethod enum above StarMetrics**

Add before line 50 in `src/analysis/mod.rs`:

```rust
/// Method used to measure this star's PSF.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FitMethod {
    /// Free-beta Moffat (8 params) — highest accuracy.
    FreeMoffat,
    /// Fixed-beta Moffat (7 params) — field median beta.
    FixedMoffat,
    /// Gaussian fallback (7 params).
    Gaussian,
    /// Windowed moments — lowest accuracy, flagged unreliable.
    Moments,
}
```

**Step 2: Add `fit_method` field to StarMetrics**

Add after the `beta` field (line 76):

```rust
    /// Method used to fit this star's PSF.
    pub fit_method: FitMethod,
```

**Step 3: Add FitMethod to MeasuredStar in metrics.rs**

In `src/analysis/metrics.rs`, add to the `MeasuredStar` struct (after `beta` field, line 22):

```rust
    /// Method used for this measurement.
    pub fit_method: super::FitMethod,
```

**Step 4: Update all MeasuredStar construction sites**

Every place that creates a `MeasuredStar` in `metrics.rs` must now include `fit_method`. Search for `MeasuredStar {` — there are construction sites in `measure_with_moments`, `measure_with_gaussian_fit`, and `measure_with_moffat_fit`. Set accordingly:
- `measure_with_moments` → `fit_method: super::FitMethod::Moments`
- `measure_with_gaussian_fit` → `fit_method: super::FitMethod::Gaussian`
- `measure_with_moffat_fit` → `fit_method: super::FitMethod::FixedMoffat` or `FreeMoffat` depending on whether `fixed_beta` is `Some`

**Step 5: Update StarMetrics construction in mod.rs**

In `run_analysis()` (around line 600-650 in mod.rs), where `StarMetrics` is constructed from `MeasuredStar`, pass through `fit_method: m.fit_method`.

**Step 6: Update lib.rs re-export**

In `src/lib.rs` line 16, add `FitMethod` to the re-export:

```rust
pub use analysis::{AnalysisConfig, AnalysisResult, FitMethod, ImageAnalyzer, StarMetrics};
```

**Step 7: Build and test**

Run: `cargo test`
Expected: All existing tests pass (FitMethod is additive, no breaking changes yet).

**Step 8: Commit**

```
feat: add FitMethod enum to StarMetrics for PSF fit provenance
```

---

### Task 2: Fitting Engine Convergence Improvements

**Files:**
- Modify: `src/analysis/fitting.rs:1-531`

**Step 1: Write failing test for Cholesky retry**

Add to the test module in `fitting.rs` (after line 778):

```rust
#[test]
fn test_cholesky_retry_on_near_singular() {
    // Very faint star (low SNR) that may cause near-singular Jacobian.
    // Should still converge with lambda retry.
    let mut samples = Vec::new();
    let sigma = 2.0_f64;
    let amplitude = 5.0_f64; // Very faint
    let noise_bg = 100.0_f64;
    for y in -8..=8 {
        for x in -8..=8 {
            let r2 = (x * x + y * y) as f64;
            let val = noise_bg + amplitude * (-0.5 * r2 / (sigma * sigma)).exp();
            samples.push(PixelSample { x: x as f64, y: y as f64, value: val });
        }
    }
    let result = fit_gaussian_2d(&samples, noise_bg, amplitude, 0.0, 0.0, sigma, sigma, 0.0);
    assert!(result.converged, "Faint star should converge with retry");
    assert!((result.sigma_x - sigma).abs() < 0.5, "sigma_x should be close");
}
```

**Step 2: Run test to verify it passes or establishes baseline**

Run: `cargo test test_cholesky_retry -- --nocapture`

**Step 3: Update constants**

In `fitting.rs`, change lines 4-5:

```rust
const MAX_ITER: usize = 50;
const CONV_TOL: f64 = 1e-6;
```

**Step 4: Add Cholesky retry to `lm_solve_2d`**

Replace lines 150-153 (the Cholesky call + break) with:

```rust
        let delta = {
            let mut result = cholesky_solve(&mat, &jtr, np);
            // Retry with increased lambda if Cholesky fails (near-singular)
            if result.is_none() {
                for _ in 0..3 {
                    lambda *= 10.0;
                    for p in 0..np {
                        mat[p * np + p] = jtj[p * np + p] + lambda * jtj[p * np + p].max(1e-12);
                    }
                    result = cholesky_solve(&mat, &jtr, np);
                    if result.is_some() { break; }
                }
            }
            match result {
                Some(d) => d,
                None => break,
            }
        };
```

**Step 5: Add dual convergence criterion to `lm_solve_2d`**

Replace lines 192-198 (convergence check) with:

```rust
        // Dual convergence: parameter norm OR residual improvement
        let param_norm = params[..np].iter().map(|p| p * p).sum::<f64>().sqrt();
        let delta_norm = delta.iter().map(|d| d * d).sum::<f64>().sqrt();
        let param_converged = delta_norm / param_norm.max(1e-12) < CONV_TOL;
        let residual_converged = if best_cost > 0.0 {
            (prev_cost - best_cost).abs() / best_cost < 1e-4
        } else {
            false
        };
        if param_converged || residual_converged {
            converged = true;
            break;
        }
```

Note: Need to track `prev_cost`. Add `let mut prev_cost = best_cost;` before the loop, and update `prev_cost = best_cost;` before the cost comparison in the gain-ratio block.

**Step 6: Add centroid drift and amplitude bounds to `lm_solve_2d`**

After the sigma positivity clamps (lines 160-165), add:

```rust
        // Centroid drift guard
        let dx = new_params[2] - params[2];
        let dy = new_params[3] - params[3];
        if dx * dx + dy * dy > 4.0 { // > 2 px drift
            break;
        }
```

**Step 7: Tighten sigma bounds**

After the sigma positivity check, add validation:

```rust
        // Sigma bounds: reject unphysical values
        if new_params[4] < 0.5 || new_params[5] < 0.5 {
            break;
        }
```

The upper bound (5× initial) is checked post-convergence, not per-iteration, to allow the optimizer room.

**Step 8: Apply same changes to `lm_solve_moffat_impl`**

Apply identical changes to `lm_solve_moffat_impl` (lines 399-531):
- Cholesky retry (replace lines 480-483)
- Dual convergence criterion (replace lines 522-527)
- Centroid drift guard (after alpha positivity clamps)
- Alpha minimum bound: 0.5 px
- Beta bounds: change 0.5/25.0 (lines 494-495) to 1.5/10.0

**Step 9: Run all tests**

Run: `cargo test`
Expected: All tests pass. Existing tests use well-conditioned data, so relaxed tolerance still passes.

**Step 10: Commit**

```
feat: improve LM convergence — Cholesky retry, dual criterion, tighter bounds
```

---

### Task 3: Dilated A Trous Wavelet (MRS Layers 2+)

**Files:**
- Modify: `src/analysis/convolution.rs:302-357`
- Modify: `src/analysis/background.rs:216-269`

**Step 1: Write failing test for dilated B3 spline**

Add to the test module in `convolution.rs`:

```rust
#[test]
fn test_b3_spline_dilated_layer2() {
    // Layer 2 should smooth more than layer 1 (wider kernel).
    let w = 64;
    let h = 64;
    let mut data = vec![0.0_f32; w * h];
    // Point source at center
    data[32 * w + 32] = 1000.0;

    let mut out1 = vec![0.0_f32; w * h];
    b3_spline_smooth(&data, w, h, &mut out1);
    let peak1 = out1[32 * w + 32];

    let mut out2 = vec![0.0_f32; w * h];
    b3_spline_smooth_dilated(&data, w, h, &mut out2, 2);
    let peak2 = out2[32 * w + 32];

    // Dilated (wider) kernel should produce a lower peak
    assert!(peak2 < peak1, "Layer 2 peak {} should be < layer 1 peak {}", peak2, peak1);
}
```

**Step 2: Run test to verify it fails**

Run: `cargo test test_b3_spline_dilated_layer2`
Expected: FAIL — `b3_spline_smooth_dilated` does not exist.

**Step 3: Implement `b3_spline_smooth_dilated`**

Add after `b3_spline_smooth` in `convolution.rs`:

```rust
/// Dilated B3-spline smoothing for à trous wavelet layers 2+.
///
/// Same kernel [1/16, 1/4, 3/8, 1/4, 1/16] but with spacing `2^(layer-1)`
/// between taps. Layer 1 has spacing 1 (use `b3_spline_smooth`).
/// Layer 2 has spacing 2, layer 3 has spacing 4, etc.
pub fn b3_spline_smooth_dilated(
    data: &[f32],
    width: usize,
    height: usize,
    output: &mut [f32],
    layer: usize,
) {
    use rayon::prelude::*;

    const K: [f32; 5] = [1.0 / 16.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 4.0, 1.0 / 16.0];
    let spacing = 1_i32 << (layer - 1); // 2^(layer-1)

    let mut hpass = vec![0.0_f32; width * height];

    // Horizontal pass with reflected boundary
    hpass
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            let src = &data[y * width..(y + 1) * width];
            for x in 0..width {
                let mut sum = 0.0_f32;
                for (ki, &kv) in K.iter().enumerate() {
                    let sx = x as i32 + (ki as i32 - 2) * spacing;
                    let sx = if sx < 0 {
                        -sx
                    } else if sx >= width as i32 {
                        2 * width as i32 - 2 - sx
                    } else {
                        sx
                    } as usize;
                    sum += src[sx.min(width - 1)] * kv;
                }
                row[x] = sum;
            }
        });

    // Vertical pass with reflected boundary
    output
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            for x in 0..width {
                let mut sum = 0.0_f32;
                for (ki, &kv) in K.iter().enumerate() {
                    let sy = y as i32 + (ki as i32 - 2) * spacing;
                    let sy = if sy < 0 {
                        -sy
                    } else if sy >= height as i32 {
                        2 * height as i32 - 2 - sy
                    } else {
                        sy
                    } as usize;
                    sum += hpass[sy.min(height - 1) * width + x] * kv;
                }
                row[x] = sum;
            }
        });
}
```

**Step 4: Run test to verify it passes**

Run: `cargo test test_b3_spline_dilated_layer2`
Expected: PASS

**Step 5: Update `estimate_noise_mrs` in background.rs to support multiple layers**

Current implementation (lines 216-269) only does layer 1. Modify to iterate:

```rust
pub fn estimate_noise_mrs(
    data: &[f32],
    width: usize,
    height: usize,
    noise_layers: usize,
) -> f32 {
    let n = width * height;
    let layers = noise_layers.max(1);

    // Working buffers
    let mut current = data.to_vec();
    let mut smoothed = vec![0.0_f32; n];
    let mut wavelet_coeffs = vec![0.0_f32; n];

    for layer in 1..=layers {
        if layer == 1 {
            super::convolution::b3_spline_smooth(&current, width, height, &mut smoothed);
        } else {
            super::convolution::b3_spline_smooth_dilated(&current, width, height, &mut smoothed, layer);
        }

        // Wavelet coefficients = current - smoothed
        for i in 0..n {
            wavelet_coeffs[i] = current[i] - smoothed[i];
        }

        // For next layer, the smoothed version becomes input
        current.copy_from_slice(&smoothed);
    }

    // Noise from the finest layer's wavelet coefficients
    // (always use layer 1 coefficients for noise — higher layers capture structure)
    // Re-compute layer 1 if we went deeper
    if layers > 1 {
        super::convolution::b3_spline_smooth(data, width, height, &mut smoothed);
        for i in 0..n {
            wavelet_coeffs[i] = data[i] - smoothed[i];
        }
    }

    // σ = 1.4826 × MAD(wavelet_coefficients)
    let mut abs_coeffs: Vec<f32> = wavelet_coeffs.iter().map(|&v| v.abs()).collect();
    abs_coeffs.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_abs = abs_coeffs[abs_coeffs.len() / 2];
    1.4826 * median_abs
}
```

Note: The multi-layer decomposition is computed but currently only layer 1's coefficients are used for noise. Higher layers are available for future use (e.g., structure detection). The key improvement is that with layers > 1, the smoothed output separates structure at multiple scales — useful for background estimation in images with nebulosity.

**Step 6: Run all tests**

Run: `cargo test`
Expected: PASS

**Step 7: Commit**

```
feat: add dilated à trous wavelet for MRS layers 2+
```

---

### Task 4: Unified Background Pipeline

**Files:**
- Modify: `src/analysis/background.rs`
- Modify: `src/analysis/mod.rs`

**Step 1: Write test for auto-tuned cell size**

Add to tests in `background.rs`:

```rust
#[test]
fn test_auto_cell_size() {
    assert_eq!(auto_cell_size(1024, 768), 32);   // 1024/32 = 32
    assert_eq!(auto_cell_size(4096, 2160), 128);  // 4096/32 = 128
    assert_eq!(auto_cell_size(256, 256), 16);     // max(16, 256/32=8) = 16
}
```

**Step 2: Run test to verify it fails**

Run: `cargo test test_auto_cell_size`
Expected: FAIL — function does not exist.

**Step 3: Implement `auto_cell_size`**

Add as a public function in `background.rs`:

```rust
/// Compute auto-tuned mesh cell size from image dimensions.
/// Targets ~32 cells across the longest axis, minimum 16 pixels per cell.
pub fn auto_cell_size(width: usize, height: usize) -> usize {
    let longest = width.max(height);
    (longest / 32).max(16)
}
```

**Step 4: Run test to verify it passes**

Run: `cargo test test_auto_cell_size`
Expected: PASS

**Step 5: Make MRS noise always-on in `run_analysis`**

In `mod.rs`, `run_analysis()` (around line 430-450), find where background estimation branches on `self.config.background_mesh_size` and `self.config.noise_layers`. Replace with:

- Always call `estimate_background_mesh` with `auto_cell_size(width, height)`
- Always call `estimate_noise_mrs` (use `self.config.noise_layers.max(1)`)
- Remove the `estimate_background` (global-only) call path

The global `background` and `noise` values are derived from the mesh result:
```rust
let cell_size = background::auto_cell_size(width, height);
let bg_result = background::estimate_background_mesh(&lum, width, height, cell_size);
let noise = background::estimate_noise_mrs(&lum, width, height, self.config.noise_layers.max(1));
let background = bg_result.background; // Already the median of the grid
```

**Step 6: Run tests**

Run: `cargo test`
Expected: PASS (mesh background produces similar results to global for test images).

**Step 7: Commit**

```
feat: unify background to always mesh-grid with auto-tuned cell size + MRS noise
```

---

### Task 5: Moffat Evaluation Function

**Files:**
- Modify: `src/analysis/fitting.rs`

**Step 1: Write failing test**

Add to tests in `fitting.rs`:

```rust
#[test]
fn test_evaluate_moffat_2d() {
    // Evaluate a known Moffat at center — should equal B + A
    let params = [100.0, 5000.0, 0.0, 0.0, 3.0, 3.0, 0.0, 3.0]; // B, A, x0, y0, ax, ay, theta, beta
    let val = evaluate_moffat_2d(&params, 0.0, 0.0);
    assert!((val - 5100.0).abs() < 0.01, "At center: B + A = 5100, got {}", val);

    // At distance: should be less than B + A
    let val_off = evaluate_moffat_2d(&params, 3.0, 0.0);
    assert!(val_off < val, "Off-center should be less than center");
    assert!(val_off > 100.0, "Off-center should be above background");
}
```

**Step 2: Run test to verify it fails**

Run: `cargo test test_evaluate_moffat_2d`
Expected: FAIL — function does not exist.

**Step 3: Implement `evaluate_moffat_2d`**

Add as a public function in `fitting.rs`:

```rust
/// Evaluate the 2D Moffat model at a single point.
/// params: [B, A, x0, y0, alpha_x, alpha_y, theta, beta]
pub fn evaluate_moffat_2d(params: &[f64], x: f64, y: f64) -> f64 {
    let (cos_t, sin_t) = (params[6].cos(), params[6].sin());
    let dx = x - params[2];
    let dy = y - params[3];
    let u = dx * cos_t + dy * sin_t;
    let v = -dx * sin_t + dy * cos_t;
    let ax = params[4];
    let ay = params[5];
    let q = (u / ax).powi(2) + (v / ay).powi(2);
    params[0] + params[1] * (1.0 + q).powf(-params[7])
}
```

**Step 4: Run test to verify it passes**

Run: `cargo test test_evaluate_moffat_2d`
Expected: PASS

**Step 5: Commit**

```
feat: add evaluate_moffat_2d for PSF subtraction
```

---

### Task 6: Iterative PSF Subtraction in Detection

**Files:**
- Modify: `src/analysis/detection.rs`

**Step 1: Write failing test**

Add to tests in `detection.rs`:

```rust
#[test]
fn test_iterative_psf_subtraction() {
    // Two stars close together, one bright, one faint
    let w = 64;
    let h = 64;
    let mut data = vec![100.0_f32; w * h]; // background = 100

    // Bright star at (20, 32)
    let sigma_bright = 3.0_f32;
    for y in 0..h {
        for x in 0..w {
            let dx = x as f32 - 20.0;
            let dy = y as f32 - 32.0;
            data[y * w + x] += 8000.0 * (-0.5 * (dx * dx + dy * dy) / (sigma_bright * sigma_bright)).exp();
        }
    }
    // Faint star at (26, 32) — only 6 pixels away, overlapping
    let sigma_faint = 3.0_f32;
    for y in 0..h {
        for x in 0..w {
            let dx = x as f32 - 26.0;
            let dy = y as f32 - 32.0;
            data[y * w + x] += 2000.0 * (-0.5 * (dx * dx + dy * dy) / (sigma_faint * sigma_faint)).exp();
        }
    }

    let params = super::detection::DetectionParams {
        detection_sigma: 3.0,
        min_star_area: 3,
        max_star_area: 500,
        saturation_limit: 60000.0,
    };
    let stars = super::detection::detect_stars(&data, w, h, 100.0, None, None, &params);
    assert!(stars.len() >= 2, "Should detect both stars, got {}", stars.len());

    // Check both centroids are reasonable
    let mut found_bright = false;
    let mut found_faint = false;
    for s in &stars {
        if (s.x - 20.0).abs() < 2.0 && (s.y - 32.0).abs() < 2.0 { found_bright = true; }
        if (s.x - 26.0).abs() < 2.0 && (s.y - 32.0).abs() < 2.0 { found_faint = true; }
    }
    assert!(found_bright, "Should find bright star near (20, 32)");
    assert!(found_faint, "Should find faint star near (26, 32)");
}
```

**Step 2: Run test to establish baseline**

Run: `cargo test test_iterative_psf_subtraction`

This test may already pass with the existing deblending. If so, the test validates the current behavior and the refactoring to iterative subtraction must maintain it.

**Step 3: Add `subtract_moffat_model` helper to detection.rs**

```rust
/// Subtract a fitted Moffat model from a stamp region of the image.
/// Returns a modified copy of the data within the stamp bounds.
fn subtract_moffat_model(
    data: &mut [f32],
    width: usize,
    height: usize,
    params: &[f64], // [B, A, x0, y0, ax, ay, theta, beta]
) {
    let cx = params[2];
    let cy = params[3];
    let ax = params[4].max(1.0);
    let ay = params[5].max(1.0);
    let radius = (5.0 * ax.max(ay)) as usize;

    let x_min = (cx as isize - radius as isize).max(0) as usize;
    let x_max = (cx as usize + radius + 1).min(width);
    let y_min = (cy as isize - radius as isize).max(0) as usize;
    let y_max = (cy as usize + radius + 1).min(height);

    for y in y_min..y_max {
        for x in x_min..x_max {
            let model_val = super::fitting::evaluate_moffat_2d(params, x as f64, y as f64);
            // Subtract only the star component (model - background)
            let star_component = (model_val - params[0]) as f32;
            data[y * width + x] -= star_component.max(0.0);
        }
    }
}
```

**Step 4: Integrate iterative subtraction into deblending**

In `detect_stars` (or `process_component`), when a multi-peak blob is encountered (the current Voronoi path), replace with:

1. Sort peaks by convolution value (already done)
2. For brightest peak: extract stamp, fit Moffat via `fitting::fit_moffat_2d`, subtract model
3. For next peak: re-extract stamp from subtracted data, fit again
4. Fall back to existing Voronoi if Moffat fit fails

This is a refactor of the existing deblending path in `process_component` (around lines 330-420 in detection.rs). The exact integration depends on how `process_component` handles multi-peak blobs — the key change is replacing Voronoi pixel assignment with subtract-and-remeasure.

**Step 5: Run tests**

Run: `cargo test`
Expected: All detection tests pass. The close-pair test should produce equal or better centroid accuracy.

**Step 6: Commit**

```
feat: replace Voronoi deblending with iterative PSF subtraction
```

---

### Task 7: Unified `measure_star` with Fallback Chain

**Files:**
- Modify: `src/analysis/metrics.rs`

**Step 1: Write failing test for fallback chain**

Add to tests in `metrics.rs`:

```rust
#[test]
fn test_measure_star_fallback_chain() {
    // Create a well-formed Gaussian star — should succeed with Moffat
    let w = 64;
    let h = 64;
    let sigma = 3.0_f32;
    let (data, _) = make_gaussian_stamp(w, h, 32.0, 32.0, sigma, sigma, 0.0, 5000.0, 100.0);

    let star = super::detection::DetectedStar {
        x: 32.0, y: 32.0, peak: 5000.0, flux: 50000.0,
        area: 50, theta: 0.0, eccentricity: 0.0,
    };

    let measured = measure_stars(&data, w, h, &[star], 100.0, None, true, true, None, None);
    assert_eq!(measured.len(), 1);
    // Should use Moffat (free or fixed depending on config)
    assert!(
        measured[0].fit_method == super::FitMethod::FreeMoffat
        || measured[0].fit_method == super::FitMethod::FixedMoffat,
        "Expected Moffat fit, got {:?}", measured[0].fit_method
    );
}
```

**Step 2: Refactor `measure_single_star` to use fallback chain**

Replace the current branching logic in `measure_single_star` (lines 58-143) which checks `use_gaussian_fit` and `use_moffat` flags. Instead:

```rust
fn measure_single_star(
    data: &[f32],
    width: usize,
    height: usize,
    star: &DetectedStar,
    background: f32,
    bg_map: Option<&[f32]>,
    field_beta: Option<f64>,   // field-wide median beta from pass 1
    green_mask: Option<&[bool]>,
    local_noise: f32,          // for amplitude validation
) -> Option<MeasuredStar> {
    // Try free-beta Moffat first (if no field_beta provided)
    // Then fixed-beta Moffat
    // Then Gaussian
    // Then moments (last resort)

    if field_beta.is_none() {
        if let Some(result) = measure_with_moffat_fit(data, width, height, star, background, bg_map, green_mask, None) {
            if result.fwhm > 0.0 {
                return Some(result); // fit_method = FreeMoffat
            }
        }
    }

    if let Some(beta) = field_beta {
        if let Some(result) = measure_with_moffat_fit(data, width, height, star, background, bg_map, green_mask, Some(beta)) {
            if result.fwhm > 0.0 {
                return Some(result); // fit_method = FixedMoffat
            }
        }
    }

    if let Some(result) = measure_with_gaussian_fit(data, width, height, star, background, bg_map, green_mask) {
        if result.fwhm > 0.0 {
            return Some(result); // fit_method = Gaussian
        }
    }

    // Last resort: moments
    measure_with_moments(data, width, height, star, background, bg_map, green_mask)
    // fit_method = Moments
}
```

**Step 3: Update `measure_stars` signature**

Remove the `use_gaussian_fit`, `use_moffat`, and `fixed_beta` parameters. Replace with `field_beta: Option<f64>`:

```rust
pub fn measure_stars(
    data: &[f32],
    width: usize,
    height: usize,
    stars: &[DetectedStar],
    background: f32,
    bg_map: Option<&[f32]>,
    green_mask: Option<&[bool]>,
    field_beta: Option<f64>,
    noise: f32,
) -> Vec<MeasuredStar> {
```

Pass 1 calls with `field_beta: None` (free Moffat on calibration stars).
Pass 2 calls with `field_beta: Some(median_beta)`.

**Step 4: Update all callers of `measure_stars`**

In `mod.rs` `run_analysis()` and `debug.rs`, update call sites to match the new signature.

**Step 5: Run tests**

Run: `cargo test`
Expected: All tests pass.

**Step 6: Commit**

```
feat: unify measure_star with Moffat→Gaussian→moments fallback chain
```

---

### Task 8: Calibration Star Selection & Two-Pass Pipeline

**Files:**
- Modify: `src/analysis/mod.rs` (run_analysis)

This is the largest task — restructuring the orchestration in `run_analysis()`.

**Step 1: Add calibration star selection function**

Add to `mod.rs`:

```rust
/// Select bright, clean, isolated stars for PSF model calibration.
fn select_calibration_stars(
    stars: &[detection::DetectedStar],
    max_count: usize,
) -> Vec<&detection::DetectedStar> {
    stars
        .iter()
        .filter(|s| {
            s.eccentricity < 0.5  // Not a trail/cosmic
            && s.area >= 5        // Not a hot pixel
        })
        // Stars are already sorted by flux (descending) from detection
        .take(max_count)
        .collect()
}
```

**Step 2: Restructure `run_analysis` into two-pass flow**

Replace the body of `run_analysis()` with:

```
// 1. Background (mesh + MRS)
let cell_size = background::auto_cell_size(width, height);
let bg_result = background::estimate_background_mesh(&lum, width, height, cell_size);
let noise = background::estimate_noise_mrs(&lum, width, height, self.config.noise_layers.max(1));
let background = bg_result.background;
let bg_map = bg_result.background_map.as_deref();
let noise_map = bg_result.noise_map.as_deref();

// 2. Pass 1: Discovery
let kernel_sigma = 3.0_f32 / FWHM_FACTOR;  // Initial conservative kernel
let detection_params = DetectionParams { ... };
let pass1_stars = detection::detect_stars(&lum, width, height, background, bg_map, noise_map, &detection_params);

// Select calibration stars
let calibration = select_calibration_stars(&pass1_stars, 100);

// Free-beta Moffat on calibration stars
let calibration_measured = metrics::measure_stars(
    &lum, width, height, &calibration_as_owned, background, bg_map, green_mask, None, noise,
);

// Derive field-wide PSF model
let field_fwhm = sigma_clipped_median(&calibration_measured.iter().map(|s| s.fwhm).collect::<Vec<_>>());
let field_beta = sigma_clipped_median(&calibration_measured.iter().filter_map(|s| s.beta).collect::<Vec<_>>());

// Source-mask background re-estimation
let source_mask = build_source_mask(&calibration_measured, width, height, field_fwhm);
let bg_result = background::estimate_background_mesh_masked(&lum, width, height, cell_size, &source_mask);
let background = bg_result.background;
let bg_map = bg_result.background_map.as_deref();

// 3. Pass 2: Full detection with refined kernel
let pass2_sigma = field_fwhm / FWHM_FACTOR;
let pass2_stars = detection::detect_stars(&lum, width, height, background, bg_map, noise_map, &detection_params_pass2);

// Fixed-beta Moffat on all stars
let measured = metrics::measure_stars(
    &lum, width, height, &pass2_stars, background, bg_map, green_mask, Some(field_beta as f64), noise,
);

// 4. Metrics
// ... SNR computation, statistics, trail detection (existing code)
```

**Step 3: Add `build_source_mask` helper**

```rust
fn build_source_mask(
    stars: &[metrics::MeasuredStar],
    width: usize,
    height: usize,
    median_fwhm: f32,
) -> Vec<bool> {
    let mut mask = vec![false; width * height];
    let mask_radius = (2.5 * median_fwhm) as isize;
    for star in stars {
        let cx = star.x as isize;
        let cy = star.y as isize;
        for dy in -mask_radius..=mask_radius {
            for dx in -mask_radius..=mask_radius {
                if dx * dx + dy * dy <= mask_radius * mask_radius {
                    let px = cx + dx;
                    let py = cy + dy;
                    if px >= 0 && px < width as isize && py >= 0 && py < height as isize {
                        mask[py as usize * width + px as usize] = true;
                    }
                }
            }
        }
    }
    mask
}
```

**Step 4: Run tests**

Run: `cargo test`

**Step 5: Commit**

```
feat: restructure analysis into two-pass pipeline with calibration stars
```

---

### Task 9: Remove SNR dB

**Files:**
- Modify: `src/analysis/snr.rs` — delete `compute_snr_db` function (lines 176-203)
- Modify: `src/analysis/mod.rs` — remove `snr_db` field from `AnalysisResult` (line 106), remove the call to `compute_snr_db` in `run_analysis`
- Modify: `src/bin/debug.rs` — remove any references to `snr_db`

**Step 1: Remove `snr_db` field from `AnalysisResult`**

Delete line 105-106 from `mod.rs`:
```rust
    /// Image-wide SNR in decibels: 20 × log10(mean_signal / noise).
    pub snr_db: f32,
```

**Step 2: Remove `compute_snr_db` from snr.rs**

Delete the function (lines 176-203) and its test.

**Step 3: Remove all references**

Search for `snr_db` across the codebase and remove:
- Assignment in `run_analysis` in `mod.rs`
- Display in `debug.rs` (`cmd_pipeline`, `cmd_measure`, `cmd_compare`)

**Step 4: Run tests**

Run: `cargo test`
Expected: PASS

**Step 5: Commit**

```
refactor: remove redundant SNR dB metric
```

---

### Task 10: Simplify Config API

**Files:**
- Modify: `src/analysis/mod.rs`

**Step 1: Update `AnalysisConfig` struct**

Replace the struct (lines 128-147) with:

```rust
pub struct AnalysisConfig {
    detection_sigma: f32,
    min_star_area: usize,
    max_star_area: usize,
    saturation_fraction: f32,
    max_stars: usize,
    apply_debayer: bool,
    trail_r_squared_threshold: f32,
    /// MRS wavelet noise layers (default 1).
    noise_layers: usize,
}
```

Removed fields: `max_measure`, `use_gaussian_fit`, `background_mesh_size`, `use_moffat_fit`, `iterative_background`, `moffat_beta`, `max_distortion`.

**Step 2: Remove deprecated builder methods**

Delete these methods from the `impl ImageAnalyzer` block:
- `without_gaussian_fit` (lines 224-227)
- `with_background_mesh` (lines 230-233)
- `with_iterative_background` (lines 253-256)
- `with_moffat_beta` (lines 271-274)
- `with_max_distortion` (lines 279-282)
- `with_moffat_fit` (lines 287-290)
- `without_moffat_fit` (lines 293-296)
- `with_max_measure` (lines 207-210)

**Step 3: Rename `with_mrs_noise` to `with_mrs_layers`**

```rust
pub fn with_mrs_layers(mut self, layers: usize) -> Self {
    self.config.noise_layers = layers;
    self
}
```

**Step 4: Update `new()` defaults**

```rust
pub fn new() -> Self {
    Self {
        config: AnalysisConfig {
            detection_sigma: 5.0,
            min_star_area: 5,
            max_star_area: 2000,
            saturation_fraction: 0.95,
            max_stars: 200,
            apply_debayer: true,
            trail_r_squared_threshold: 0.5,
            noise_layers: 1,
        },
        thread_pool: None,
    }
}
```

**Step 5: Update debug.rs to match**

Update all config construction in `debug.rs` to use the new API. Remove CLI flags that map to deleted methods (`--no-moffat`, `--fixed-beta`, `--max-ecc`, `--mesh`). The mesh is now auto-tuned. Moffat is always on.

Keep `--mrs` flag, mapping to `.with_mrs_layers(n)`.

**Step 6: Run tests**

Run: `cargo test`
Expected: PASS

**Step 7: Commit**

```
refactor: simplify config API — remove algorithm-selection knobs
```

---

### Task 11: Update Debug CLI

**Files:**
- Modify: `src/bin/debug.rs`

**Step 1: Remove deprecated CLI flags**

Remove from `Options` and `parse_args()`:
- `--no-moffat` → Moffat is always on
- `--fixed-beta` → Beta is auto-derived
- `--max-ecc` → Internal to calibration selection
- `--mesh` → Auto-tuned

**Step 2: Update subcommand handlers**

- `cmd_background`: Use `auto_cell_size` instead of `opts.mesh`
- `cmd_measure`: Remove `use_moffat` / `fixed_beta` branching, use unified `measure_stars`
- `cmd_fit`: Keep as-is (this is for single-star deep inspection, both Gaussian and Moffat are shown)
- `cmd_pipeline`: Reflect the new two-pass structure
- `cmd_dump`: Remove `snr_db` column
- `cmd_compare`: Remove `snr_db` from comparison

**Step 3: Run build**

Run: `cargo build --features debug-pipeline`
Expected: Compiles clean.

**Step 4: Commit**

```
refactor: update debug CLI for unified pipeline
```

---

### Task 12: Integration Tests on All Test Files

**Files:**
- Create: `tests/analysis_integration.rs`

**Step 1: Write integration tests**

```rust
use astroimage::ImageAnalyzer;

#[test]
fn test_analysis_cocoon_fits() {
    let result = ImageAnalyzer::new()
        .with_detection_sigma(5.0)
        .analyze("tests/cocoon.fits")
        .expect("Failed to analyze cocoon.fits");

    assert!(result.stars.len() > 10, "Should detect stars in cocoon");
    assert!(result.median_fwhm > 0.0, "FWHM should be positive");
    assert!(result.noise > 0.0, "Noise should be positive");
    assert!(result.snr_weight > 0.0, "SNR weight should be positive");
    assert!(result.median_beta.is_some(), "Should have Moffat beta");
}

#[test]
fn test_analysis_mono_fits() {
    let result = ImageAnalyzer::new()
        .analyze("tests/mono.fits")
        .expect("Failed to analyze mono.fits");

    assert!(result.stars.len() > 0, "Should detect stars in mono");
    assert!(result.median_fwhm > 0.0);
    // Mono comet image — may have extended object, fewer clean stars
}

#[test]
fn test_analysis_osc_fits() {
    let result = ImageAnalyzer::new()
        .analyze("tests/osc.fits")
        .expect("Failed to analyze osc.fits");

    assert!(result.stars.len() > 0, "Should detect stars in OSC");
    assert!(result.source_channels == 3, "OSC should produce 3-channel");
}

#[test]
fn test_analysis_xisf() {
    let result = ImageAnalyzer::new()
        .analyze("tests/test.xisf")
        .expect("Failed to analyze test.xisf");

    assert!(result.stars.len() > 0, "Should detect stars in XISF");
    assert!(result.median_fwhm > 0.0);
}

#[test]
fn test_fit_method_populated() {
    let result = ImageAnalyzer::new()
        .analyze("tests/cocoon.fits")
        .expect("Failed");

    for star in &result.stars {
        // Every star should have a fit method
        assert!(
            star.fit_method == astroimage::FitMethod::FreeMoffat
            || star.fit_method == astroimage::FitMethod::FixedMoffat
            || star.fit_method == astroimage::FitMethod::Gaussian
            || star.fit_method == astroimage::FitMethod::Moments,
        );
    }

    // Most stars should be FixedMoffat (pass 2)
    let moffat_count = result.stars.iter()
        .filter(|s| s.fit_method == astroimage::FitMethod::FixedMoffat)
        .count();
    assert!(
        moffat_count > result.stars.len() / 2,
        "Most stars should use FixedMoffat, got {}/{}", moffat_count, result.stars.len()
    );
}
```

**Step 2: Run integration tests**

Run: `cargo test --test analysis_integration`
Expected: All PASS

**Step 3: Commit**

```
test: add integration tests for unified analysis pipeline
```

---

### Task 13: Final Cleanup & Verification

**Step 1: Run full test suite**

Run: `cargo test`
Expected: All PASS

**Step 2: Run clippy**

Run: `cargo clippy`
Expected: No warnings related to analysis modules

**Step 3: Build release**

Run: `cargo build --release`
Expected: Clean build

**Step 4: Build debug binary**

Run: `cargo build --features debug-pipeline`
Expected: Clean build

**Step 5: Manual verification on all test files**

```bash
cargo run --features debug-pipeline -- pipeline tests/cocoon.fits --save /tmp/analysis-test
cargo run --features debug-pipeline -- pipeline tests/mono.fits --save /tmp/analysis-test
cargo run --features debug-pipeline -- pipeline tests/osc.fits --save /tmp/analysis-test
cargo run --features debug-pipeline -- pipeline tests/test.xisf --save /tmp/analysis-test
```

Visually inspect the outputs for each file.

**Step 6: Commit**

```
chore: final cleanup and verification of analysis pipeline redesign
```
