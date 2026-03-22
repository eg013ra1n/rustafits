# Analysis Speed & Robustness Optimization — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce per-frame analysis time from ~1.25s to ~0.5-0.6s on M4 and improve star selection robustness by rejecting nebula knots, edge stars, and blended stars before expensive PSF fitting.

**Architecture:** Add three detection-stage morphological filters (sharpness, concentration index, edge margin) that run in Pass 2 only, a moments-based pre-screening gate before LM fitting in the measurement stage, and spatially-balanced star selection with a lower default measure_cap of 500.

**Tech Stack:** Rust, rayon (parallelism), existing DAOFIND/CCL/LM infrastructure

**Spec:** `docs/superpowers/specs/2026-03-22-analysis-speed-robustness-design.md`

**Test files (ALWAYS test on ALL):**
- `tests/cocoon.fits` — dense star field
- `tests/mono.fits` — monochrome (comet, bright extended object)
- `tests/osc.fits` — OSC Bayer pattern
- `tests/test.xisf` — XISF format

---

## Chunk 1: Detection-Stage Filters

### Task 1: Add `field_fwhm` parameter to `detect_stars()`

**Files:**
- Modify: `src/analysis/detection.rs:48-58` (function signature)
- Modify: `src/analysis/mod.rs:421-425` (Pass 1 call site)
- Modify: `src/analysis/mod.rs:526-530` (Pass 2 call site)

- [ ] **Step 1: Add `Option<f32>` parameter to `detect_stars` signature**

In `src/analysis/detection.rs`, add `field_fwhm: Option<f32>` as the last parameter:

```rust
pub fn detect_stars(
    data: &[f32],
    width: usize,
    height: usize,
    background: f32,
    noise: f32,
    bg_map: Option<&[f32]>,
    noise_map: Option<&[f32]>,
    params: &DetectionParams,
    fwhm: f32,
    field_fwhm: Option<f32>,  // NEW: enables Pass 2 filters when Some
) -> Vec<DetectedStar> {
```

- [ ] **Step 2: Update all call sites to pass `None` / `Some`**

In `src/analysis/mod.rs`, Pass 1 call (~line 421):
```rust
detection::detect_stars(
    &lum, width, height,
    bg_result.background, bg_result.noise,
    bg_map, noise_map, &det_params, initial_fwhm,
    None,  // Pass 1: no field_fwhm yet
)
```

Pass 2 call (~line 526):
```rust
detection::detect_stars(
    &lum, width, height,
    bg_result.background, bg_result.noise,
    bg_map, noise_map, &det_params, final_fwhm,
    Some(field_fwhm),  // Pass 2: enable sharpness/CI/edge filters
)
```

Update test call sites in `detection.rs` tests to pass `None` (all 5 tests: `test_detect_synthetic_stars`, `test_reject_hot_pixels`, `test_skip_blended_close_stars`, `test_no_duplicate_detections`, `test_no_deblend_extended_object`).

Also, derive `Clone` on `DetectedStar` (add `#[derive(Clone)]` before the struct definition at line 4) — this avoids fragile manual field-by-field copying when the struct is used in spatial grid selection later.

- [ ] **Step 3: Verify compilation**

Run: `cargo check`
Expected: Compiles cleanly.

- [ ] **Step 4: Run existing tests**

Run: `cargo test`
Expected: All existing tests pass (no behavioral change yet).

- [ ] **Step 5: Commit**

```bash
git add src/analysis/detection.rs src/analysis/mod.rs
git commit -m "refactor: add field_fwhm parameter to detect_stars for Pass 2 filters"
```

---

### Task 2: Add DAOFIND sharpness filter

**Files:**
- Modify: `src/analysis/detection.rs:244-255` (main validation loop)
- Test: `src/analysis/detection.rs` (tests module)

- [ ] **Step 1: Write failing test for sharpness rejection**

Add to the `tests` module in `detection.rs`:

```rust
#[test]
fn test_sharpness_rejects_diffuse_blob() {
    // A broad, diffuse blob (sigma=15) should be rejected by sharpness filter
    // when field_fwhm is provided (Pass 2 mode).
    let width = 200;
    let height = 200;
    let background = 1000.0;
    let noise = 30.0;

    // Real star + diffuse nebula knot
    let star_defs = vec![
        (50.0, 50.0, 5000.0, 2.0),   // real star, sigma=2
        (130.0, 130.0, 3000.0, 15.0), // diffuse blob, sigma=15
    ];
    let data = make_star_field(width, height, &star_defs, background, noise);

    let params = DetectionParams {
        detection_sigma: 3.0,
        min_star_area: 5,
        max_star_area: 2000,
        saturation_limit: 0.95 * 65535.0,
    };

    let field_fwhm = 2.0 * 2.3548; // sigma=2 → fwhm≈4.7

    // With field_fwhm (Pass 2): blob should be rejected by sharpness
    let stars_pass2 = detect_stars(
        &data, width, height, background, noise,
        None, None, &params, 3.0, Some(field_fwhm),
    );

    // Real star should be found
    let has_real = stars_pass2.iter().any(|s| {
        let dx = s.x - 50.0;
        let dy = s.y - 50.0;
        (dx * dx + dy * dy).sqrt() < 3.0
    });
    assert!(has_real, "Real star at (50,50) should be detected");

    // Diffuse blob should NOT be found
    let has_blob = stars_pass2.iter().any(|s| {
        let dx = s.x - 130.0;
        let dy = s.y - 130.0;
        (dx * dx + dy * dy).sqrt() < 20.0
    });
    assert!(!has_blob, "Diffuse blob at (130,130) should be rejected by sharpness filter");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test test_sharpness_rejects_diffuse_blob -- --nocapture`
Expected: FAIL — blob is not yet rejected.

- [ ] **Step 3: Implement sharpness filter**

In `detect_stars()`, after the component validation loop (after line 252, where `process_component` returns), add the sharpness check. The filter goes between `process_component` returning and pushing to `stars`:

Replace the loop body at lines 244-255:

```rust
let mut stars = Vec::new();
for (label, pixels) in &components {
    let peaks = &peak_map[label];
    let max_conv = peaks.iter().map(|p| p.2).fold(0.0_f32, f32::max);
    let significant_peaks = peaks.iter().filter(|p| p.2 > 0.3 * max_conv).count();
    if significant_peaks > 1 {
        continue;
    }
    if let Some(star) = process_component(pixels, data, width, height, background, bg_map, params) {
        // ── Pass 2 filters (when field_fwhm is provided) ──
        if let Some(_ff) = field_fwhm {
            // Sharpness: (conv[peak] - mean(conv[8_neighbors])) / conv[peak]
            // Use the strongest peak in this component for the lookup.
            let &(px, py, _) = peaks.iter().max_by(|a, b| a.2.total_cmp(&b.2)).unwrap();
            if px >= 1 && py >= 1 && px < width - 1 && py < height - 1 {
                let cp = conv[py * width + px];
                if cp > 0.0 {
                    let neighbors_sum =
                        conv[(py - 1) * width + px - 1]
                        + conv[(py - 1) * width + px]
                        + conv[(py - 1) * width + px + 1]
                        + conv[py * width + px - 1]
                        + conv[py * width + px + 1]
                        + conv[(py + 1) * width + px - 1]
                        + conv[(py + 1) * width + px]
                        + conv[(py + 1) * width + px + 1];
                    let sharpness = (cp - neighbors_sum / 8.0) / cp;
                    if sharpness < 0.2 || sharpness > 0.9 {
                        continue;
                    }
                }
            }
        }

        stars.push(star);
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test test_sharpness_rejects_diffuse_blob -- --nocapture`
Expected: PASS

- [ ] **Step 5: Run all existing tests**

Run: `cargo test`
Expected: All tests pass. Existing tests use `field_fwhm: None` so sharpness filter is inactive.

- [ ] **Step 6: Commit**

```bash
git add src/analysis/detection.rs
git commit -m "feat: add DAOFIND sharpness filter for Pass 2 detection"
```

---

### Task 3: Add concentration index filter

**Files:**
- Modify: `src/analysis/detection.rs:244-255` (main validation loop, after sharpness)

- [ ] **Step 1: Write failing test for concentration index**

Add to `detection.rs` tests:

```rust
#[test]
fn test_concentration_index_rejects_extended_source() {
    // A moderately extended source (sigma=8) that passes sharpness
    // but fails concentration index.
    let width = 200;
    let height = 200;
    let background = 1000.0;
    let noise = 30.0;

    let star_defs = vec![
        (50.0, 50.0, 5000.0, 2.0),  // real star
        (130.0, 130.0, 4000.0, 8.0), // extended source (sigma=8)
    ];
    let data = make_star_field(width, height, &star_defs, background, noise);

    let params = DetectionParams {
        detection_sigma: 3.0,
        min_star_area: 5,
        max_star_area: 2000,
        saturation_limit: 0.95 * 65535.0,
    };

    let field_fwhm = 2.0 * 2.3548;

    let stars = detect_stars(
        &data, width, height, background, noise,
        None, None, &params, 3.0, Some(field_fwhm),
    );

    // Real star found
    let has_real = stars.iter().any(|s| {
        let dx = s.x - 50.0;
        let dy = s.y - 50.0;
        (dx * dx + dy * dy).sqrt() < 3.0
    });
    assert!(has_real, "Real star at (50,50) should be detected");

    // Extended source should be rejected
    let has_extended = stars.iter().any(|s| {
        let dx = s.x - 130.0;
        let dy = s.y - 130.0;
        (dx * dx + dy * dy).sqrt() < 10.0
    });
    assert!(!has_extended, "Extended source at (130,130) should be rejected by CI filter");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test test_concentration_index_rejects_extended_source -- --nocapture`
Expected: FAIL — extended source passes.

- [ ] **Step 3: Implement concentration index filter**

After the sharpness check in the main loop, add CI computation. Insert inside the `if let Some(_ff) = field_fwhm` block, after the sharpness check:

```rust
// Concentration index: flux(1σ) / flux(3σ) from raw bg-subtracted pixels
let ff = _ff;  // rename for use
let field_sigma = ff / 2.3548;
let ci_inner_r = field_sigma;
let ci_outer_r = 3.0 * field_sigma;
let ci_inner_r_sq = ci_inner_r * ci_inner_r;
let ci_outer_r_sq = ci_outer_r * ci_outer_r;
let ci_radius = ci_outer_r.ceil() as i32;
let scx = star.x.round() as i32;
let scy = star.y.round() as i32;

let mut flux_inner = 0.0_f64;
let mut flux_outer = 0.0_f64;
for dy in -ci_radius..=ci_radius {
    let py = scy + dy;
    if py < 0 || py >= height as i32 { continue; }
    for dx in -ci_radius..=ci_radius {
        let px = scx + dx;
        if px < 0 || px >= width as i32 { continue; }
        let r_sq = (dx * dx + dy * dy) as f32;
        if r_sq > ci_outer_r_sq { continue; }
        let bg = bg_map.map_or(background, |m| m[py as usize * width + px as usize]);
        let val = (data[py as usize * width + px as usize] - bg).max(0.0) as f64;
        flux_outer += val;
        if r_sq <= ci_inner_r_sq {
            flux_inner += val;
        }
    }
}
if flux_outer > 0.0 {
    let ci = flux_inner / flux_outer;
    if ci < 0.15 || ci > 0.75 {
        continue;
    }
}
```

Note: the `_ff` rename needs to replace the `_ff` in the sharpness section. Change `if let Some(_ff) = field_fwhm` to `if let Some(ff) = field_fwhm` and use `ff` for both sharpness and CI.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test test_concentration_index_rejects_extended_source -- --nocapture`
Expected: PASS

- [ ] **Step 5: Run all tests**

Run: `cargo test`
Expected: All pass.

- [ ] **Step 6: Commit**

```bash
git add src/analysis/detection.rs
git commit -m "feat: add concentration index filter for Pass 2 detection"
```

---

### Task 4: Add FWHM-proportional edge margin

**Files:**
- Modify: `src/analysis/detection.rs` (main validation loop, inside the `field_fwhm` block)

- [ ] **Step 1: Write failing test for edge margin**

Add to `detection.rs` tests:

```rust
#[test]
fn test_edge_margin_rejects_near_border_stars() {
    // Stars within 2*FWHM of the edge should be rejected in Pass 2.
    let width = 100;
    let height = 100;
    let background = 1000.0;
    let noise = 30.0;

    let field_fwhm = 3.0 * 2.3548; // sigma=3 → fwhm≈7.06
    let margin = (2.0 * field_fwhm).max(8.0); // ≈14.1

    // Star near edge (at x=10, well within margin of ~14.1) and one safely inside
    let star_defs = vec![
        (10.0, 50.0, 5000.0, 3.0),  // near left edge
        (50.0, 50.0, 5000.0, 3.0),  // safely inside
    ];
    let data = make_star_field(width, height, &star_defs, background, noise);

    let params = DetectionParams::default();

    // Pass 2 (with field_fwhm): edge star should be rejected
    let stars = detect_stars(
        &data, width, height, background, noise,
        None, None, &params, field_fwhm, Some(field_fwhm),
    );

    let has_center = stars.iter().any(|s| {
        let dx = s.x - 50.0;
        let dy = s.y - 50.0;
        (dx * dx + dy * dy).sqrt() < 3.0
    });
    assert!(has_center, "Center star at (50,50) should be detected");

    let has_edge = stars.iter().any(|s| {
        let dx = s.x - 10.0;
        let dy = s.y - 50.0;
        (dx * dx + dy * dy).sqrt() < 3.0
    });
    assert!(!has_edge, "Edge star at (10,50) should be rejected by edge margin (margin={:.1})", margin);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test test_edge_margin_rejects_near_border_stars -- --nocapture`
Expected: FAIL — edge star currently passes (only pixel-level border check in `process_component`).

- [ ] **Step 3: Implement edge margin filter**

Add inside the `if let Some(ff) = field_fwhm` block, after CI and before `stars.push(star)`:

```rust
// Edge margin: reject stars whose centroid is within 2*FWHM of any edge
let margin = (2.0 * ff).max(8.0);
if star.x < margin || star.y < margin
    || star.x > (width as f32 - margin)
    || star.y > (height as f32 - margin)
{
    continue;
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test test_edge_margin_rejects_near_border_stars -- --nocapture`
Expected: PASS

- [ ] **Step 5: Run all tests**

Run: `cargo test`
Expected: All pass.

- [ ] **Step 6: Commit**

```bash
git add src/analysis/detection.rs
git commit -m "feat: add FWHM-proportional edge margin for Pass 2 detection"
```

---

## Chunk 2: Moments Pre-Screening Gate

### Task 5: Add moments pre-screening to `measure_single_star`

**Files:**
- Modify: `src/analysis/metrics.rs:38-59` (measure_stars signature)
- Modify: `src/analysis/metrics.rs:61-146` (measure_single_star)
- Modify: `src/analysis/mod.rs:446-453` (calibration call site)
- Modify: `src/analysis/mod.rs:596-603` (main measurement call site)

- [ ] **Step 1: Write failing test for moments screening**

Note: `measure_single_star` takes a full image (`data`, `width`, `height`) and extracts its own stamp using the star's absolute coordinates. Tests must embed stars in a full-sized synthetic image, not pass raw stamps.

Add a new test module section at the bottom of `metrics.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    /// Generate a synthetic image with a single Gaussian star on a flat background.
    fn make_star_image(
        width: usize, height: usize,
        star_x: f32, star_y: f32,
        amp: f32, sigma: f32,
        background: f32,
    ) -> Vec<f32> {
        let mut data = vec![background; width * height];
        let r = (4.0 * sigma).ceil() as i32;
        let inv_2s2 = 1.0 / (2.0 * sigma * sigma);
        for dy in -r..=r {
            for dx in -r..=r {
                let px = star_x as i32 + dx;
                let py = star_y as i32 + dy;
                if px >= 0 && px < width as i32 && py >= 0 && py < height as i32 {
                    let ddx = px as f32 - star_x;
                    let ddy = py as f32 - star_y;
                    data[py as usize * width + px as usize] +=
                        amp * (-inv_2s2 * (ddx * ddx + ddy * ddy)).exp();
                }
            }
        }
        data
    }

    #[test]
    fn test_moments_screening_rejects_extended() {
        // An extended source (sigma=10) embedded in a 200x200 image.
        // With field_fwhm ~4.7 (sigma=2), this source has moment-FWHM >> 2.5 * field_fwhm.
        let width = 200;
        let height = 200;
        let background = 1000.0;
        let field_fwhm = 4.7_f32;

        let data = make_star_image(width, height, 100.0, 100.0, 5000.0, 10.0, background);

        let star = DetectedStar {
            x: 100.0, y: 100.0, peak: 5000.0, flux: 50000.0,
            area: 200, theta: 0.0, eccentricity: 0.0,
        };

        let result = measure_single_star(
            &data, width, height, &star, background, None, None,
            Some(1.0), 25, 1e-4, 5,
            Some(field_fwhm), false,
        );

        assert!(result.is_none(), "Extended source should be rejected by moments screening");
    }

    #[test]
    fn test_moments_screening_accepts_real_star() {
        // A normal star (sigma=2) matching the field_fwhm.
        let width = 200;
        let height = 200;
        let background = 1000.0;
        let field_fwhm = 4.7_f32; // sigma=2 → fwhm≈4.7

        let data = make_star_image(width, height, 100.0, 100.0, 5000.0, 2.0, background);

        let star = DetectedStar {
            x: 100.0, y: 100.0, peak: 5000.0, flux: 30000.0,
            area: 50, theta: 0.0, eccentricity: 0.1,
        };

        let result = measure_single_star(
            &data, width, height, &star, background, None, None,
            Some(1.0), 25, 1e-4, 5,
            Some(field_fwhm), false,
        );

        assert!(result.is_some(), "Real star should pass moments screening");
    }
}
```

Note: This test won't compile yet because `measure_single_star` doesn't have the new parameters. That's expected — it tests the interface change.

- [ ] **Step 2: Update `measure_stars` and `measure_single_star` signatures**

In `src/analysis/metrics.rs`, update `measure_stars`:

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
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
    field_fwhm: Option<f32>,       // NEW
    possibly_trailed: bool,         // NEW
) -> Vec<MeasuredStar> {
    use rayon::prelude::*;

    let result: Vec<MeasuredStar> = stars
        .par_iter()
        .filter_map(|star| {
            measure_single_star(
                data, width, height, star, background, bg_map, green_mask,
                field_beta, max_iter, conv_tol, max_rejects,
                field_fwhm, possibly_trailed,
            )
        })
        .collect();

    // Safety fallback: if moments screening rejected >50% of candidates,
    // the eccentricity gate may be too aggressive (extreme wind shake).
    // Re-run with screening disabled (field_fwhm = None).
    if field_fwhm.is_some() && !possibly_trailed
        && result.len() < stars.len() / 2
        && stars.len() >= 10
    {
        return stars
            .par_iter()
            .filter_map(|star| {
                measure_single_star(
                    data, width, height, star, background, bg_map, green_mask,
                    field_beta, max_iter, conv_tol, max_rejects,
                    None, false,  // disable screening
                )
            })
            .collect();
    }

    result
}
```

Update `measure_single_star`:

```rust
fn measure_single_star(
    data: &[f32],
    width: usize,
    height: usize,
    star: &DetectedStar,
    background: f32,
    bg_map: Option<&[f32]>,
    green_mask: Option<&[bool]>,
    field_beta: Option<f64>,
    max_iter: usize,
    conv_tol: f64,
    max_rejects: usize,
    field_fwhm: Option<f32>,       // NEW
    possibly_trailed: bool,         // NEW
) -> Option<MeasuredStar> {
```

- [ ] **Step 3: Update call sites in `mod.rs`**

Calibration call (~line 446):
```rust
let cal_measured = metrics::measure_stars(
    &lum, width, height, &cal_owned,
    bg_result.background,
    bg_result.background_map.as_deref(),
    green_mask,
    None, // free-beta Moffat
    50, 1e-6, 5,
    None,  // no screening for calibration
    false, // not trailed (calibration doesn't care)
);
```

Main measurement call (~line 596):
```rust
let mut measured = metrics::measure_stars(
    &lum, width, height, to_measure,
    bg_result.background, bg_map_ref,
    green_mask, field_beta,
    self.config.fit_max_iter,
    self.config.fit_tolerance,
    self.config.fit_max_rejects,
    Some(field_fwhm),    // enable moments screening
    possibly_trailed,     // bypass ecc gate on trailed frames
);
```

- [ ] **Step 4: Implement the moments screening gate**

In `measure_single_star`, after stamp extraction (after line ~108, the `stamp` vec is built) and before the fallback chain (before line ~133), insert:

```rust
// ── Moments pre-screening gate ──────────────────────────────────
// Cheap moments check to reject non-stellar sources before expensive LM fitting.
// Only active when field_fwhm is provided (Pass 2 measurement).
if let Some(ff) = field_fwhm {
    // Quick moments: 2-iteration centroid + second-order moments
    let mut mcx = rel_cx as f64;
    let mut mcy = rel_cy as f64;
    let fit_r = 5.0_f64.max(4.0 * estimated_sigma as f64);
    let fit_r_sq = fit_r * fit_r;

    let mut mxx = 0.0_f64;
    let mut myy = 0.0_f64;
    let mut mxy = 0.0_f64;

    for _ in 0..2 {
        mxx = 0.0; myy = 0.0; mxy = 0.0;
        let mut sw = 0.0_f64;
        let mut swx = 0.0_f64;
        let mut swy = 0.0_f64;

        for sy in 0..stamp_h {
            for sx in 0..stamp_w {
                let dx = sx as f64 - mcx;
                let dy = sy as f64 - mcy;
                if dx * dx + dy * dy > fit_r_sq { continue; }
                let val = stamp[sy * stamp_w + sx].max(0.0) as f64;
                if val <= 0.0 { continue; }
                sw += val;
                swx += val * sx as f64;
                swy += val * sy as f64;
                mxx += val * dx * dx;
                myy += val * dy * dy;
                mxy += val * dx * dy;
            }
        }

        if sw < 1e-10 { return None; }
        mcx = swx / sw;
        mcy = swy / sw;
        mxx /= sw;
        myy /= sw;
        mxy /= sw;
    }

    // Derive moment-FWHM and eccentricity
    let trace = mxx + myy;
    let det = mxx * myy - mxy * mxy;
    let disc = (trace * trace - 4.0 * det).max(0.0).sqrt();
    let l1 = (trace + disc) * 0.5;
    let l2 = (trace - disc) * 0.5;

    if l1 > 0.0 && l2 > 0.0 {
        let moment_fwhm = (2.3548 * (l1.sqrt() * l2.sqrt())) as f32;
        let moment_ecc = (1.0 - l2 / l1).max(0.0).sqrt() as f32;

        // Reject if moment-FWHM outside [0.7, 2.5] * field_fwhm
        if moment_fwhm < 0.7 * ff || moment_fwhm > 2.5 * ff {
            return None;
        }

        // Reject high eccentricity (unless trailed).
        // Note: the 50% fallback (spec requirement) is handled at the
        // measure_stars level, not here. If >50% of candidates are rejected
        // by this gate, measure_stars re-runs with the ecc gate disabled.
        if !possibly_trailed && moment_ecc > 0.8 {
            return None;
        }

        // Peak/mean sharpness
        let mean_flux = star.flux / star.area as f32;
        if mean_flux > 0.0 {
            let peak_sharpness = star.peak / mean_flux;
            if peak_sharpness < 0.3 || peak_sharpness > 5.0 {
                return None;
            }
        }
    }
}
```

- [ ] **Step 5: Run tests**

Run: `cargo test`
Expected: All pass including new moments screening tests.

- [ ] **Step 6: Commit**

```bash
git add src/analysis/metrics.rs src/analysis/mod.rs
git commit -m "feat: add moments-based pre-screening gate before LM fitting"
```

---

## Chunk 3: Reduced Measure Cap with Spatial Sampling

### Task 6: Lower default measure_cap to 500

**Files:**
- Modify: `src/analysis/mod.rs:188` (default value)

- [ ] **Step 1: Change default**

At line 188, change:
```rust
measure_cap: 2000,
```
to:
```rust
measure_cap: 500,
```

- [ ] **Step 2: Run tests**

Run: `cargo test`
Expected: All pass.

- [ ] **Step 3: Commit**

```bash
git add src/analysis/mod.rs
git commit -m "perf: lower default measure_cap from 2000 to 500"
```

---

### Task 7: Add spatial grid selection before measurement

**Files:**
- Modify: `src/analysis/mod.rs:585-603` (before measure_stars call)

- [ ] **Step 1: Write test for spatial distribution**

Add a test in `mod.rs` (or a new test file) that verifies spatial distribution. Since this is an integration-level behavior, we'll verify it works correctly via the existing test files in the final integration test. For now, just implement and verify no regression.

- [ ] **Step 2: Implement spatial grid selection**

In `src/analysis/mod.rs`, replace the measure cap slicing logic (~lines 588-594):

```rust
// Apply measure cap with spatial grid balancing.
// Divide image into 4×4 grid, round-robin select from each cell
// to ensure spatial coverage across the field.
let effective_cap = if self.config.measure_cap == 0 {
    detected.len()
} else {
    self.config.measure_cap
};

let to_measure: Vec<&detection::DetectedStar> = if detected.len() <= effective_cap {
    // Fewer stars than cap — use all
    detected.iter().collect()
} else {
    // Spatial grid: 4×4 cells, round-robin by flux rank
    debug_assert!(
        detected.windows(2).all(|w| w[0].flux >= w[1].flux),
        "detected stars must be sorted by flux descending"
    );
    const GRID_N: usize = 4;
    let cell_w = width as f32 / GRID_N as f32;
    let cell_h = height as f32 / GRID_N as f32;
    let mut buckets: Vec<Vec<&detection::DetectedStar>> =
        vec![Vec::new(); GRID_N * GRID_N];

    for star in &detected {
        let gx = ((star.x / cell_w) as usize).min(GRID_N - 1);
        let gy = ((star.y / cell_h) as usize).min(GRID_N - 1);
        buckets[gy * GRID_N + gx].push(star);
    }

    let mut selected: Vec<&detection::DetectedStar> = Vec::with_capacity(effective_cap);
    let mut idx = vec![0usize; GRID_N * GRID_N];
    loop {
        let mut added_any = false;
        for cell in 0..(GRID_N * GRID_N) {
            if selected.len() >= effective_cap { break; }
            if idx[cell] < buckets[cell].len() {
                selected.push(buckets[cell][idx[cell]]);
                idx[cell] += 1;
                added_any = true;
            }
        }
        if !added_any || selected.len() >= effective_cap { break; }
    }
    selected
};

// Convert references to owned for measure_stars (which expects &[DetectedStar]).
// Uses Clone derived on DetectedStar (added in Task 1).
let to_measure_owned: Vec<detection::DetectedStar> = to_measure
    .iter()
    .map(|s| (*s).clone())
    .collect();

let mut measured = metrics::measure_stars(
    &lum, width, height, &to_measure_owned,
    bg_result.background, bg_map_ref,
    green_mask, field_beta,
    self.config.fit_max_iter,
    self.config.fit_tolerance,
    self.config.fit_max_rejects,
    Some(field_fwhm),
    possibly_trailed,
);
```

- [ ] **Step 3: Run tests**

Run: `cargo test`
Expected: All pass.

- [ ] **Step 4: Commit**

```bash
git add src/analysis/mod.rs
git commit -m "feat: add spatial grid selection for balanced star measurement"
```

---

## Chunk 4: Integration Testing & Validation

### Task 8: Run full integration tests on all test files

**Files:**
- Test files: `tests/cocoon.fits`, `tests/mono.fits`, `tests/osc.fits`, `tests/test.xisf`

- [ ] **Step 1: Build release binary**

Run: `cargo build --release --features debug-pipeline`
Expected: Builds cleanly.

- [ ] **Step 2: Run on all 4 test files and capture metrics**

Run each and compare output metrics against pre-optimization baseline:

```bash
cargo run --release --features debug-pipeline -- measure tests/cocoon.fits 2>/dev/null
cargo run --release --features debug-pipeline -- measure tests/mono.fits 2>/dev/null
cargo run --release --features debug-pipeline -- measure tests/osc.fits 2>/dev/null
cargo run --release --features debug-pipeline -- measure tests/test.xisf 2>/dev/null
```

Expected: median_fwhm, median_eccentricity, median_hfr, median_snr within 5% of previous values (allowing for improved robustness to reject junk stars that were skewing stats).

- [ ] **Step 3: Generate annotated outputs for visual inspection**

```bash
cargo run --release -- tests/cocoon.fits /tmp/cocoon_opt.jpg --annotate
cargo run --release -- tests/mono.fits /tmp/mono_opt.jpg --annotate
cargo run --release -- tests/osc.fits /tmp/osc_opt.jpg --annotate
cargo run --release -- tests/test.xisf /tmp/xisf_opt.jpg --annotate
```

Visually verify: no good stars missing from annotations, no obvious nebula knots or edge artifacts annotated.

- [ ] **Step 4: Run all unit tests**

Run: `cargo test`
Expected: All pass.

- [ ] **Step 5: Commit any test adjustments needed**

Only if threshold adjustments were needed during validation:
```bash
git add src/analysis/ && git commit -m "test: adjust thresholds after optimization validation"
```

---

### Task 9: Performance measurement

- [ ] **Step 1: Benchmark before vs after**

If you have access to a batch of FITS frames (M42 dataset):
```bash
time for f in /path/to/frames/*.fits; do cargo run --release --features debug-pipeline -- measure "$f" 2>/dev/null; done
```

Otherwise, benchmark on the 4 test files with timing:
```bash
time cargo run --release --features debug-pipeline -- measure tests/cocoon.fits 2>/dev/null
```

Expected: ~2x speedup or better compared to pre-optimization timing.

- [ ] **Step 2: Document results**

Add timing comparison as a comment in the commit or note for the user.

- [ ] **Step 3: Final commit**

```bash
git add src/analysis/ && git commit -m "perf: analysis pipeline speed & robustness optimization complete"
```
