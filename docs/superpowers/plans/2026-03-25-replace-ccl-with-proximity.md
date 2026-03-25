# Replace CCL with Peak Proximity Blend Detection

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the expensive Connected Component Labeling (CCL) in star detection with a lightweight peak-proximity blend test, eliminating ~200ms per frame on dense fields while maintaining identical blend rejection behavior.

**Architecture:** CCL currently does two full-image scans (26M pixels each) to label connected regions, then rejects components with multiple peaks. The replacement uses the NMS spatial grid (already built) to check if any other peak is within 2×FWHM — if so, both are blended. Per-star metrics (area, flux, centroid, theta, ecc) move from CCL pixel lists to stamp-based computation around each peak position.

**Tech Stack:** Rust, `cargo test`, `rustafits-debug pipeline`

---

## File Map

| File | Changes |
|------|---------|
| `src/analysis/detection.rs` | Replace CCL (lines 151-240) + component processing (248-325) with proximity blend test + stamp-based star processing |
| `src/analysis/metrics.rs` | Change `estimated_sigma` from area-based to FWHM-based (line 104) |

No API changes — `DetectedStar` struct keeps the same fields. No test file changes needed (no tests reference `DetectedStar` fields directly).

---

## Task 1: Replace CCL with proximity-based blend detection

**Files:**
- Modify: `src/analysis/detection.rs:150-345`

The current flow after NMS (line 149) is:
1. CCL first pass — label 26M pixels (lines 151-196)
2. CCL second pass — resolve labels (lines 198-203)
3. Map peaks to labels (lines 205-229)
4. Collect pixels per component — scan 26M pixels again (lines 232-240)
5. Reject multi-peak components (lines 248-254)
6. Process each component: centroid, area, flux, theta, ecc from CCL pixel list (lines 255, 352-489)
7. Pass 2 filters (lines 256-321)
8. Dedup (lines 327-343)

The new flow after NMS:
1. **Proximity blend test** — for each peak, check spatial grid for any other peak within `blend_radius = 2.0 * fwhm`. Mark as blended if found.
2. **Stamp-based star processing** — for each non-blended peak, extract a small stamp and compute centroid, peak, flux, area, theta, ecc directly from the stamp (same logic as current `process_component` lines 374-478, but using stamp pixels instead of CCL pixel list).
3. Pass 2 filters (unchanged)
4. Dedup (unchanged)

- [ ] **Step 1: Write the proximity blend check**

Replace everything from line 151 (CCL start) through line 325 (end of Pass 2 filters, before dedup) with:

```rust
    // ── Stage 2: Proximity-based blend rejection + stamp metrics ────────
    // Two peaks within blend_radius share overlapping PSFs — skip both.
    // This replaces full-image CCL (O(W×H)) with O(N) spatial grid queries.
    let blend_radius = 2.0 * fwhm;
    let blend_radius_sq = blend_radius * blend_radius;
    let blend_cell = blend_radius.ceil() as usize;
    let blend_grid_w = (width + blend_cell - 1) / blend_cell;
    let blend_grid_h = (height + blend_cell - 1) / blend_cell;

    // Build spatial grid of all surviving peaks for blend queries
    let mut blend_grid: Vec<Vec<usize>> = vec![Vec::new(); blend_grid_w * blend_grid_h];
    for (i, &(px, py, _)) in peak_positions.iter().enumerate() {
        let gx = (px / blend_cell).min(blend_grid_w - 1);
        let gy = (py / blend_cell).min(blend_grid_h - 1);
        blend_grid[gy * blend_grid_w + gx].push(i);
    }

    // Mark blended peaks
    let mut blended = vec![false; peak_positions.len()];
    for (i, &(px, py, _)) in peak_positions.iter().enumerate() {
        if blended[i] { continue; }
        let gx = px / blend_cell;
        let gy = py / blend_cell;
        let search_r = (blend_radius / blend_cell as f32).ceil() as usize + 1;
        let gx_lo = gx.saturating_sub(search_r);
        let gy_lo = gy.saturating_sub(search_r);
        let gx_hi = (gx + search_r + 1).min(blend_grid_w);
        let gy_hi = (gy + search_r + 1).min(blend_grid_h);

        for ngy in gy_lo..gy_hi {
            for ngx in gx_lo..gx_hi {
                for &j in &blend_grid[ngy * blend_grid_w + ngx] {
                    if j <= i { continue; }
                    let (jx, jy, _) = peak_positions[j];
                    let dx = px as f32 - jx as f32;
                    let dy = py as f32 - jy as f32;
                    if dx * dx + dy * dy <= blend_radius_sq {
                        blended[i] = true;
                        blended[j] = true;
                    }
                }
            }
        }
    }

    // Process non-blended peaks: stamp-based metrics
    let stamp_r = ((2.0 * fwhm) as i32).max(5);
    let mut stars = Vec::new();
    for (i, &(px, py, _conv_val)) in peak_positions.iter().enumerate() {
        if blended[i] { continue; }
        if let Some(star) = process_peak_stamp(
            px, py, stamp_r, data, width, height, background, bg_map, params,
        ) {
            // ── Pass 2 filters (when field_fwhm is provided) ──
            if let Some(ff) = field_fwhm {
                // 1. Sharpness
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
                        if sharpness < 0.1 || sharpness > 0.9 {
                            continue;
                        }
                    }
                }

                // 2. Concentration index
                let field_sigma = ff / 2.3548;
                let ci_inner_r_sq = field_sigma * field_sigma;
                let ci_outer_r = 3.0 * field_sigma;
                let ci_outer_r_sq = ci_outer_r * ci_outer_r;
                let ci_radius = ci_outer_r.ceil() as i32;
                let scx = star.x.round() as i32;
                let scy = star.y.round() as i32;

                let mut flux_inner = 0.0_f64;
                let mut flux_outer = 0.0_f64;
                for dy in -ci_radius..=ci_radius {
                    let py2 = scy + dy;
                    if py2 < 0 || py2 >= height as i32 { continue; }
                    for dx in -ci_radius..=ci_radius {
                        let px2 = scx + dx;
                        if px2 < 0 || px2 >= width as i32 { continue; }
                        let r_sq = (dx * dx + dy * dy) as f32;
                        if r_sq > ci_outer_r_sq { continue; }
                        let bg = bg_map.map_or(background, |m| m[py2 as usize * width + px2 as usize]);
                        let val = (data[py2 as usize * width + px2 as usize] - bg).max(0.0) as f64;
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

                // 3. Edge margin
                let margin = (2.0 * ff).max(8.0);
                if star.x < margin || star.y < margin
                    || star.x > (width as f32 - margin)
                    || star.y > (height as f32 - margin)
                {
                    continue;
                }
            }

            stars.push(star);
        }
    }
```

- [ ] **Step 2: Write the `process_peak_stamp` function**

Replace the old `process_component` (lines 352-489) with a new function that takes a peak position + stamp radius instead of a CCL pixel list:

```rust
/// Compute star metrics from a stamp around a peak position.
/// Replaces CCL-based process_component — uses the continuous image directly.
fn process_peak_stamp(
    peak_x: usize,
    peak_y: usize,
    stamp_r: i32,
    data: &[f32],
    width: usize,
    height: usize,
    background: f32,
    bg_map: Option<&[f32]>,
    params: &DetectionParams,
) -> Option<DetectedStar> {
    let cx_i = peak_x as i32;
    let cy_i = peak_y as i32;

    // Bounds check: stamp must fit in image
    if cx_i - stamp_r <= 0 || cy_i - stamp_r <= 0
        || cx_i + stamp_r >= width as i32 - 1
        || cy_i + stamp_r >= height as i32 - 1
    {
        return None;
    }

    // Compute peak, centroid, flux, area from stamp
    let bg_at = |x: usize, y: usize| bg_map.map_or(background, |m| m[y * width + x]);
    let low_thresh = background * 0.02; // ~2% of background as floor

    let mut peak = 0.0_f32;
    let mut raw_peak = 0.0_f32;
    let mut sum_w = 0.0_f64;
    let mut sum_wx = 0.0_f64;
    let mut sum_wy = 0.0_f64;
    let mut flux = 0.0_f64;
    let mut area = 0_usize;

    for dy in -stamp_r..=stamp_r {
        let py = (cy_i + dy) as usize;
        for dx in -stamp_r..=stamp_r {
            let px = (cx_i + dx) as usize;
            let raw = data[py * width + px];
            let bg = bg_at(px, py);
            let val = raw - bg;
            if val <= low_thresh { continue; }

            area += 1;
            if val > peak { peak = val; }
            if raw > raw_peak { raw_peak = raw; }

            let w = (val.max(0.0) as f64).powi(2);
            sum_w += w;
            sum_wx += w * px as f64;
            sum_wy += w * py as f64;
            flux += val.max(0.0) as f64;
        }
    }

    // Area filter
    if area < params.min_star_area || area > params.max_star_area {
        return None;
    }

    // Saturation check
    if raw_peak > params.saturation_limit {
        return None;
    }

    if sum_w < 1e-10 {
        return None;
    }

    let cx = (sum_wx / sum_w) as f32;
    let cy = (sum_wy / sum_w) as f32;

    // Stamp-based I-weighted second moments for theta and eccentricity
    let mut sf = 0.0_f64;
    let mut six = 0.0_f64;
    let mut siy = 0.0_f64;
    let mut sixx = 0.0_f64;
    let mut siyy = 0.0_f64;
    let mut sixy = 0.0_f64;
    for dy in -stamp_r..=stamp_r {
        let py = (cy_i + dy) as usize;
        for dx in -stamp_r..=stamp_r {
            let px = (cx_i + dx) as usize;
            let bg = bg_at(px, py);
            let v = (data[py * width + px] - bg).max(0.0) as f64;
            sf += v;
            six += v * px as f64;
            siy += v * py as f64;
            sixx += v * (px as f64) * (px as f64);
            siyy += v * (py as f64) * (py as f64);
            sixy += v * (px as f64) * (py as f64);
        }
    }
    let (theta, ecc) = if sf > 1e-10 {
        let icx = six / sf;
        let icy = siy / sf;
        let mxx = sixx / sf - icx * icx;
        let myy = siyy / sf - icy * icy;
        let mxy = sixy / sf - icx * icy;
        let t = (0.5 * (2.0 * mxy).atan2(mxx - myy)) as f32;
        let trace = mxx + myy;
        let det = mxx * myy - mxy * mxy;
        let disc = (trace * trace - 4.0 * det).max(0.0);
        let l1 = (trace + disc.sqrt()) * 0.5;
        let l2 = (trace - disc.sqrt()) * 0.5;
        let e = if l1 > 0.0 { (1.0 - l2 / l1).max(0.0).sqrt() as f32 } else { 0.0 };
        (t, e)
    } else {
        (0.0, 0.0)
    };

    Some(DetectedStar {
        x: cx,
        y: cy,
        peak,
        flux: flux as f32,
        area,
        theta,
        eccentricity: ecc,
    })
}
```

- [ ] **Step 3: Remove old CCL code**

Delete: `process_component` function (old lines 352-489), `find` function, `union` function, and the union-find helpers. These are no longer used.

- [ ] **Step 4: Update module doc comment**

Line 1: change "Star detection: DAOFIND-inspired matched filter + connected component labeling." to "Star detection: DAOFIND-inspired matched filter + proximity blend rejection."

- [ ] **Step 5: Run `cargo test`**

```bash
cargo test
```

All 86 tests must pass.

- [ ] **Step 6: Commit**

```text
perf: replace CCL with proximity-based blend detection

Eliminates two full-image scans (26M pixels each) and 104MB labels
array. Blend detection now uses O(N) spatial grid query instead of
O(W×H) connected component labeling. ~200ms faster on dense fields.
```

---

## Task 2: Switch area-based sigma estimate to FWHM-based in metrics

**Files:**
- Modify: `src/analysis/metrics.rs:103-105`

The stamp radius in `measure_single_star` currently uses CCL area to estimate sigma:
```rust
let estimated_sigma = (star.area as f32 / std::f32::consts::PI).sqrt() * 0.5;
```

With stamp-based area (instead of CCL area), the value might differ slightly. More robust: use `field_fwhm` when available, fall back to area-based estimate.

- [ ] **Step 1: Update sigma estimation**

Replace lines 103-105:
```rust
    // Estimate sigma from the star's area: area ≈ π*(2σ)² for a Gaussian at low_threshold
    let estimated_sigma = (star.area as f32 / std::f32::consts::PI).sqrt() * 0.5;
    let estimated_sigma = estimated_sigma.max(1.0).min(20.0);
```

with:
```rust
    // Estimate sigma: prefer field_fwhm if available, fall back to area-based estimate
    let estimated_sigma = if let Some(ff) = field_fwhm {
        ff / 2.3548
    } else {
        (star.area as f32 / std::f32::consts::PI).sqrt() * 0.5
    };
    let estimated_sigma = estimated_sigma.max(1.0).min(20.0);
```

Note: `field_fwhm` is already passed to `measure_single_star` (as parameter at line 100). This change makes the stamp radius consistent regardless of how `area` was computed.

- [ ] **Step 2: Run `cargo test`**

```bash
cargo test
```

- [ ] **Step 3: Commit**

```text
refactor: prefer field_fwhm for stamp radius estimation
```

---

## Task 3: Validate accuracy and performance

- [ ] **Step 1: Run PI comparison**

```bash
cargo build --features debug-pipeline --release
target/release/rustafits-debug compare tests/pix_tests/m42_compare.csv
```

Check:
- Ecc R² ≥ 0.94 (was 0.9465)
- FWHM offset within ±1% (was +0.494%)
- Noise ratio ~1.0 (was 1.001)
- Star count should be similar (may change slightly due to different blend rejection)

- [ ] **Step 2: Benchmark all reference files**

```bash
for f in tests/cocoon.fits tests/mono.fits tests/osc.fits tests/test.xisf; do
    echo "=== $(basename $f) ===" && target/release/rustafits-debug pipeline "$f" 2>&1 | grep -E "Time:|Total pipeline|detections|Measured"
done
```

Expected: detection time for cocoon should drop significantly (was 607ms). Other files proportional.

- [ ] **Step 3: Run all four reference files and compare metrics**

```bash
for f in tests/cocoon.fits tests/mono.fits tests/osc.fits tests/test.xisf; do
    echo "=== $(basename $f) ===" && target/release/rustafits-debug pipeline "$f" 2>&1 | grep -E "Median FWHM|Median ecc|Median beta|Moffat"
done
```

Values should be close to previous (small changes acceptable due to different blend rejection at boundaries).

- [ ] **Step 4: Run full test suite**

```bash
cargo test
```

- [ ] **Step 5: Generate comparison chart**

```bash
target/release/rustafits-debug compare tests/pix_tests/m42_compare.csv 2>/dev/null > /tmp/data.tsv
python3 tests/pix_tests/plot_compare.py /tmp/data.tsv /tmp/proximity_compare.png
```
