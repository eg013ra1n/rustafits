# Multi-Layer MRS Noise Estimation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement proper iterative multiresolution support (MRS) noise estimation matching PixInsight's algorithm, where `noise_layers` controls the number of wavelet scales used to build a significance mask that excludes structure from the noise estimate.

**Architecture:** The à trous wavelet transform decomposes the image into wavelet coefficient layers at increasing scales. At each layer, pixels with significant structure (|w_j| > 3σ) are masked out. After processing all layers, the noise is re-estimated from layer-1 coefficients using only unmasked (structure-free) pixels. This iterates until convergence.

**Tech Stack:** Pure Rust, existing `b3_spline_smooth` and `b3_spline_smooth_dilated` functions, rayon for parallelism.

---

### Task 1: Un-gate `b3_spline_smooth_dilated` from `#[cfg(test)]`

**Files:**
- Modify: `src/analysis/convolution.rs:368` — remove `#[cfg(test)]` attribute

**Step 1: Remove the test-only gate**

Change line 368 from:
```rust
#[cfg(test)]
pub fn b3_spline_smooth_dilated(
```
to:
```rust
pub fn b3_spline_smooth_dilated(
```

The function needs `use rayon::prelude::*;` at the top — it currently has it inside the function body. That's fine, leave it.

**Step 2: Build to verify it compiles**

Run: `cargo build`
Expected: SUCCESS (no warnings about the function being unused — it will be used in Task 2)

**Step 3: Commit**

```bash
git add src/analysis/convolution.rs
git commit -m "feat: expose b3_spline_smooth_dilated for multi-layer MRS"
```
i
---

### Task 2: Write failing unit test for multi-layer MRS

**Files:**
- Modify: `src/analysis/background.rs` — add test in `mod tests`

**Step 1: Write the failing test**

Add this test to the `mod tests` block at the bottom of `background.rs`:

```rust
#[test]
fn test_mrs_significance_masking() {
    // Image with bright star on flat background.
    // Multi-layer MRS should give a LOWER noise estimate than single-layer
    // because it masks out the star's structure from the noise sample.
    let w = 200;
    let h = 200;
    let bg_level = 1000.0_f32;
    let mut data = vec![bg_level; w * h];

    // Add Gaussian noise (deterministic pseudo-random)
    let true_sigma = 30.0_f32;
    for i in 0..data.len() {
        // Simple LCG for reproducibility
        let r = ((i as u64 * 6364136223846793005 + 1442695040888963407) >> 33) as f32;
        data[i] += true_sigma * (r / (1u64 << 31) as f32 - 0.5) * 2.0;
    }

    // Add a bright extended structure (simulated nebula patch)
    for y in 60..140 {
        for x in 60..140 {
            let dx = (x as f32 - 100.0) / 30.0;
            let dy = (y as f32 - 100.0) / 30.0;
            data[y * w + x] += 2000.0 * (-0.5 * (dx * dx + dy * dy)).exp();
        }
    }

    let noise_1layer = estimate_noise_mrs(&data, w, h, 1);
    let noise_4layer = estimate_noise_mrs(&data, w, h, 4);

    // With bright structure, 1-layer estimate will be inflated.
    // 4-layer should be closer to the true noise (30.0).
    eprintln!("1-layer noise: {:.2}, 4-layer noise: {:.2}, true: {:.2}",
        noise_1layer, noise_4layer, true_sigma);

    assert!(
        noise_4layer <= noise_1layer,
        "4-layer ({:.2}) should be <= 1-layer ({:.2}) with bright structure",
        noise_4layer, noise_1layer,
    );
    // 4-layer should be reasonably close to truth
    assert!(
        (noise_4layer - true_sigma).abs() < true_sigma * 0.5,
        "4-layer noise {:.2} should be within 50% of true sigma {:.2}",
        noise_4layer, true_sigma,
    );
}
```

**Step 2: Run to verify it fails**

Run: `cargo test test_mrs_significance_masking -- --nocapture`
Expected: FAIL — because `estimate_noise_mrs` ignores the layers parameter, both estimates will be identical.

---

### Task 3: Implement multi-layer MRS with significance masking

**Files:**
- Modify: `src/analysis/background.rs:172-229` — rewrite `estimate_noise_mrs`

**Step 1: Implement the algorithm**

Replace the entire `estimate_noise_mrs` function with:

```rust
/// Noise estimation via iterative Multiresolution Support (MRS).
///
/// Algorithm (matches PixInsight's approach):
/// 1. Compute à trous wavelet coefficients at layers 1..noise_layers
/// 2. At each layer, identify "significant" pixels (|w_j| > 3σ_j)
/// 3. Build a cumulative significance mask across all layers
/// 4. Re-estimate noise from layer-1 coefficients using only unmasked pixels
/// 5. Iterate until convergence (σ changes < 0.1%)
///
/// `noise_layers`: number of wavelet scales for significance detection (1-6).
/// Layer 1 = finest scale (pixel noise), higher layers = coarser structures.
/// More layers → better structure rejection → purer noise estimate.
/// Default: 4 (matches PixInsight). Minimum: 1.
pub fn estimate_noise_mrs(
    data: &[f32],
    width: usize,
    height: usize,
    noise_layers: usize,
) -> f32 {
    use super::convolution::{b3_spline_smooth, b3_spline_smooth_dilated};

    if data.is_empty() {
        return 0.001;
    }

    let total = width * height;
    let noise_layers = noise_layers.clamp(1, 6);

    // ── Layer 1: compute wavelet coefficients w1 = data - smooth(data) ──
    let mut smoothed = vec![0.0_f32; total];
    b3_spline_smooth(data, width, height, &mut smoothed);

    let mut w1 = vec![0.0_f32; total];
    for i in 0..total {
        w1[i] = data[i] - smoothed[i];
    }

    // If only 1 layer requested, use simple MAD on w1 (original behavior)
    if noise_layers == 1 {
        return simple_mad_noise(&w1, width, height);
    }

    // ── Build significance mask across all wavelet layers ──
    // true = significant (contains structure), excluded from noise estimate
    let mut sig_mask = vec![false; total];

    // Layer 1: mark significant pixels in w1
    let sigma1 = simple_mad_noise(&w1, width, height);
    let thresh1 = 3.0 * sigma1;
    for i in 0..total {
        if w1[i].abs() > thresh1 {
            sig_mask[i] = true;
        }
    }

    // Layers 2..N: compute dilated wavelet coefficients, mark significant pixels
    // Each layer operates on the smoothed output of the previous layer
    let mut prev_smooth = smoothed;

    for layer in 2..=noise_layers {
        let mut next_smooth = vec![0.0_f32; total];
        b3_spline_smooth_dilated(&prev_smooth, width, height, &mut next_smooth, layer);

        // Wavelet coefficients at this layer: w_j = c_{j-1} - c_j
        // Subsample for sigma estimation
        let target_samples = 500_000usize;
        let stride = ((total as f64 / target_samples as f64).sqrt() as usize).max(1);
        let border = (2_usize << (layer - 1)).min(width / 4);

        let mut wj_samples = Vec::with_capacity(target_samples);
        let mut y = border;
        while y < height.saturating_sub(border) {
            let mut x = border;
            while x < width.saturating_sub(border) {
                let idx = y * width + x;
                if !sig_mask[idx] {
                    let coeff = prev_smooth[idx] - next_smooth[idx];
                    if coeff.is_finite() {
                        wj_samples.push(coeff);
                    }
                }
                x += stride;
            }
            y += stride;
        }

        if wj_samples.len() >= 100 {
            let med = find_median(&mut wj_samples);
            let mut abs_devs: Vec<f32> = wj_samples.iter().map(|&v| (v - med).abs()).collect();
            let mad = find_median(&mut abs_devs);
            let sigma_j = 1.4826 * mad;
            let thresh_j = 3.0 * sigma_j;

            // Mark significant pixels at this scale
            for idx in 0..total {
                if !sig_mask[idx] {
                    let coeff = prev_smooth[idx] - next_smooth[idx];
                    if coeff.abs() > thresh_j {
                        sig_mask[idx] = true;
                    }
                }
            }
        }

        prev_smooth = next_smooth;
    }

    // ── Re-estimate noise from layer-1 coefficients, excluding masked pixels ──
    // Iterate until convergence
    let target_samples = 500_000usize;
    let stride = ((total as f64 / target_samples as f64).sqrt() as usize).max(1);
    let border = 2usize;

    let mut sigma = sigma1;
    for _iteration in 0..4 {
        let thresh = 3.0 * sigma;

        let mut samples = Vec::with_capacity(target_samples);
        let mut y = border;
        while y < height.saturating_sub(border) {
            let mut x = border;
            while x < width.saturating_sub(border) {
                let idx = y * width + x;
                if !sig_mask[idx] && w1[idx].abs() <= thresh {
                    samples.push(w1[idx]);
                }
                x += stride;
            }
            y += stride;
        }

        if samples.len() < 100 {
            break;
        }

        let med = find_median(&mut samples);
        let mut abs_devs: Vec<f32> = samples.iter().map(|&v| (v - med).abs()).collect();
        let mad = find_median(&mut abs_devs);
        let new_sigma = 1.4826 * mad;

        if new_sigma <= 0.0 || (new_sigma - sigma).abs() / sigma < 0.001 {
            sigma = new_sigma;
            break;
        }
        sigma = new_sigma;
    }

    sigma.max(0.001)
}

/// Simple MAD-based noise estimate from wavelet coefficients (subsampled).
fn simple_mad_noise(w1: &[f32], width: usize, height: usize) -> f32 {
    let total = width * height;
    let target_samples = 500_000usize;
    let stride = ((total as f64 / target_samples as f64).sqrt() as usize).max(1);
    let border = 2usize;

    let mut samples = Vec::with_capacity(target_samples);
    let mut y = border;
    while y < height.saturating_sub(border) {
        let mut x = border;
        while x < width.saturating_sub(border) {
            let idx = y * width + x;
            let coeff = w1[idx];
            if coeff.is_finite() {
                samples.push(coeff);
            }
            x += stride;
        }
        y += stride;
    }

    if samples.len() < 100 {
        return 0.001;
    }

    let median = find_median(&mut samples);
    let mut abs_devs: Vec<f32> = samples.iter().map(|&v| (v - median).abs()).collect();
    let mad = find_median(&mut abs_devs);
    (1.4826 * mad).max(0.001)
}
```

**Step 2: Run unit test to verify it passes**

Run: `cargo test test_mrs_significance_masking -- --nocapture`
Expected: PASS — 4-layer noise should be lower (closer to true sigma) than 1-layer.

**Step 3: Run all unit tests**

Run: `cargo test --lib`
Expected: All pass

**Step 4: Commit**

```bash
git add src/analysis/background.rs
git commit -m "feat: implement multi-layer MRS noise estimation with significance masking"
```

---

### Task 4: Update default `noise_layers` from 1 to 4

**Files:**
- Modify: `src/analysis/mod.rs:173` — change default from 1 to 4
- Modify: `src/bin/debug.rs:91` — change default from 1 to 4

**Step 1: Update library default**

In `src/analysis/mod.rs`, line 173, change:
```rust
noise_layers: 1,
```
to:
```rust
noise_layers: 4,
```

**Step 2: Update debug CLI default**

In `src/bin/debug.rs`, line 91, change:
```rust
mrs_layers: 1,
```
to:
```rust
mrs_layers: 4,
```

**Step 3: Build to verify**

Run: `cargo build`
Expected: SUCCESS

**Step 4: Commit**

```bash
git add src/analysis/mod.rs src/bin/debug.rs
git commit -m "feat: default noise_layers to 4 (matches PixInsight)"
```

---

### Task 5: Run full integration tests

**Files:** None (test-only)

**Step 1: Run all tests**

Run: `cargo test`
Expected: All 72+ tests pass. Noise values may change (they should be slightly lower or equal for images with structure).

**Step 2: Run debug CLI on test files to verify reasonable output**

Run: `cargo build --features debug-pipeline && target/debug/rustafits-debug background tests/cocoon.fits`
Verify: MRS noise is reported with `(4L)` and value is reasonable.

Run: `target/debug/rustafits-debug background tests/mono.fits`
Verify: Same — noise should be ≤ the old 1-layer value for this image (bright comet structure).

**Step 3: Commit test results if any assertion values need updating**

Only needed if integration test assertions are too tight.

---

### Task 6: Update documentation

**Files:**
- Modify: `docs/usage.md` — update `with_mrs_layers` default description
- Modify: `docs/analysis-pipeline.md` — document MRS algorithm
- Modify: `CLAUDE.md` — update if needed

**Step 1: Update docs**

In `docs/usage.md`, update the `with_mrs_layers` entry to mention default is now 4 and describe the algorithm.

In `docs/analysis-pipeline.md`, add a section describing the iterative MRS significance masking algorithm.

**Step 2: Commit**

```bash
git add docs/
git commit -m "docs: update MRS noise estimation description"
```
