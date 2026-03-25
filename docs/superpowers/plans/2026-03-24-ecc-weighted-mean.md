# Eccentricity Weighted Mean Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Switch eccentricity aggregation from sigma-clipped weighted median to sigma-clipped weighted mean, and remove the ecc≤0.8 pre-filter for eccentricity — matching PixInsight SubframeSelector's approach where fit-residual weighting alone handles outliers.

**Architecture:** Two targeted changes: (1) replace `sigma_clipped_weighted_median` with `sigma_clipped_weighted_mean` for eccentricity only (FWHM/HFR stay as weighted median), (2) remove the `ecc <= 0.8` pre-filter for eccentricity (keep it for FWHM where it prevents elongated profiles from inflating geometric-mean FWHM). Both changes apply in the library (`mod.rs`) and debug binary (`debug.rs`). The `FWHM_ECC_MAX` constant and its filter stay for FWHM/HFR — only the eccentricity usage is removed. The `moment_ecc > 0.8` detection-stage gate in `metrics.rs:198` is unrelated (pre-screening) and stays.

**Tech Stack:** Rust, `cargo test`, `rustafits-debug compare`

---

## Task 1: Add unit tests for `sigma_clipped_weighted_mean`

The function already exists at `src/analysis/mod.rs:861` but has no dedicated test.

**Files:**

- Modify: `src/analysis/mod.rs:1096-1103` (add tests next to existing `test_sigma_clipped_weighted_median_basic`)

- [ ] **Step 1: Write the tests**

Add after the existing `test_sigma_clipped_weighted_median_basic` test:

```rust
#[test]
fn test_sigma_clipped_weighted_mean_basic() {
    // With an outlier, sigma clipping should reject it, then compute weighted mean
    let vals = [3.0_f32, 3.1, 3.0, 3.2, 3.0, 100.0]; // 100.0 is outlier
    let wts = [1.0_f32; 6];
    let scwm = sigma_clipped_weighted_mean(&vals, &wts);
    assert!(scwm < 4.0, "Outlier should be clipped, got {}", scwm);
    // Weighted mean of [3.0, 3.1, 3.0, 3.2, 3.0] with equal weights = 3.06
    assert!((scwm - 3.06).abs() < 0.05, "Expected ~3.06, got {}", scwm);
}

#[test]
fn test_sigma_clipped_weighted_mean_skewed_weights() {
    // High weight on one value should pull the mean toward it
    let vals = [0.3_f32, 0.4, 0.5, 0.6];
    let wts = [10.0_f32, 1.0, 1.0, 1.0];
    let scwm = sigma_clipped_weighted_mean(&vals, &wts);
    assert!(scwm < 0.4, "Heavy weight on 0.3 should pull mean below 0.4, got {}", scwm);
}

#[test]
fn test_sigma_clipped_weighted_mean_empty() {
    assert_eq!(sigma_clipped_weighted_mean(&[], &[]), 0.0);
}

#[test]
fn test_sigma_clipped_weighted_mean_single() {
    let scwm = sigma_clipped_weighted_mean(&[0.42], &[1.0]);
    assert!((scwm - 0.42).abs() < 0.001);
}
```

- [ ] **Step 2: Run tests to verify they pass** (function already exists)

Run: `cargo test test_sigma_clipped_weighted_mean -- --nocapture`

Expected: PASS (the function is implemented, just untested)

- [ ] **Step 3: Remove `#[allow(dead_code)]` from both functions**

In `src/analysis/mod.rs`, remove `#[allow(dead_code)]` from `sigma_clipped_weighted_mean` (line 860) and `weighted_mean` (line 899). They will no longer be dead code after Task 2 calls them from library code.

- [ ] **Step 4: Run `cargo test` to confirm all pass, no dead_code warnings**

Run: `cargo test`

Expected: all pass

- [ ] **Step 5: Commit**

```text
feat: add unit tests for sigma_clipped_weighted_mean
```

---

## Task 2: Switch eccentricity to weighted mean and remove ecc≤0.8 filter in library

**Files:**

- Modify: `src/analysis/mod.rs:685-706` (eccentricity aggregation block)

- [ ] **Step 1: Replace the eccentricity aggregation block**

In `src/analysis/mod.rs`, replace lines 685-700:

```rust
        // Eccentricity: on normal frames, ecc ≤ 0.8 cutoff removes noise from
        // faint detections. On trailed frames, elongation IS the signal — bypass
        // the cutoff so the reported ecc reflects actual frame quality.
        let ecc_use_all = possibly_trailed;
        let ecc_filtered: Vec<&metrics::MeasuredStar> = if ecc_use_all {
            measured.iter().collect()
        } else {
            let filtered: Vec<&metrics::MeasuredStar> = measured.iter()
                .filter(|s| s.eccentricity <= FWHM_ECC_MAX)
                .collect();
            if filtered.len() >= 3 { filtered } else { measured.iter().collect() }
        };
        let ecc_vals: Vec<f32> = ecc_filtered.iter().map(|s| s.eccentricity).collect();
        let ecc_weights: Vec<f32> = ecc_filtered.iter()
            .map(|s| 1.0 / (1.0 + s.fit_residual))
            .collect();
```

with:

```rust
        // Eccentricity: fit-residual weighting handles outlier rejection
        // (matching PI SubframeSelector's residual-weighted mean approach).
        // No ecc pre-filter — all measured stars contribute, weighted by fit quality.
        let ecc_vals: Vec<f32> = measured.iter().map(|s| s.eccentricity).collect();
        let ecc_weights: Vec<f32> = measured.iter()
            .map(|s| 1.0 / (1.0 + s.fit_residual))
            .collect();
```

- [ ] **Step 2: Change the aggregation call from median to mean**

On the line (currently ~706) that reads:

```rust
        let median_eccentricity = sigma_clipped_weighted_median(&ecc_vals, &ecc_weights);
```

Change to:

```rust
        let median_eccentricity = sigma_clipped_weighted_mean(&ecc_vals, &ecc_weights);
```

(Field name `median_eccentricity` kept to avoid breaking public API.)

- [ ] **Step 3: Run `cargo test`**

Run: `cargo test`

Expected: all pass

- [ ] **Step 4: Commit**

```text
feat: switch eccentricity to weighted mean, remove ecc≤0.8 filter
```

---

## Task 3: Match changes in debug binary

**Files:**

- Modify: `src/bin/debug.rs:1131-1143` (ecc block in `cmd_pipeline`)
- Modify: `src/bin/debug.rs:1575-1592` (ecc block in `analyze_one_file`)
- Modify: `src/bin/debug.rs:1325-1341` (remove `wmean_ecc` from `OurResult`)
- Modify: `src/bin/debug.rs:1720-1793` (remove wmean diagnostic lines in `cmd_compare`)

- [ ] **Step 1: Update `cmd_pipeline` eccentricity block**

In `src/bin/debug.rs`, replace the ecc block at ~lines 1131-1143:

```rust
    // Ecc: trail-aware cutoff
    let ecc_filtered: Vec<&metrics::MeasuredStar> = if possibly_trailed {
        measured.iter().collect()
    } else {
        let filtered: Vec<&metrics::MeasuredStar> = measured.iter()
            .filter(|s| s.eccentricity <= FWHM_ECC_MAX)
            .collect();
        if filtered.len() >= 3 { filtered } else { measured.iter().collect() }
    };
    let ecc_vals: Vec<f32> = ecc_filtered.iter().map(|s| s.eccentricity).collect();
    let ecc_weights: Vec<f32> = ecc_filtered.iter()
        .map(|s| 1.0 / (1.0 + s.fit_residual)).collect();
    let median_ecc = sigma_clipped_weighted_median(&ecc_vals, &ecc_weights);
```

with:

```rust
    // Eccentricity: fit-residual weighted mean (no ecc pre-filter).
    let ecc_vals: Vec<f32> = measured.iter().map(|s| s.eccentricity).collect();
    let ecc_weights: Vec<f32> = measured.iter()
        .map(|s| 1.0 / (1.0 + s.fit_residual)).collect();
    let median_ecc = sigma_clipped_weighted_mean(&ecc_vals, &ecc_weights);
```

- [ ] **Step 2: Update `analyze_one_file` eccentricity block**

In `src/bin/debug.rs`, replace the ecc block at ~lines 1575-1592:

```rust
    // Eccentricity: on normal frames, ecc ≤ 0.8 cutoff removes noise from
    // faint detections. On trailed frames, elongation IS the signal — bypass
    // the cutoff so the reported ecc reflects actual frame quality.
    let ecc_use_all = possibly_trailed;
    let ecc_filtered: Vec<&metrics::MeasuredStar> = if ecc_use_all {
        measured.iter().collect()
    } else {
        let filtered: Vec<&metrics::MeasuredStar> = measured.iter()
            .filter(|s| s.eccentricity <= FWHM_ECC_MAX)
            .collect();
        if filtered.len() >= 3 { filtered } else { measured.iter().collect() }
    };
    let (ecc_vals, ecc_weights) = (
        ecc_filtered.iter().map(|s| s.eccentricity).collect::<Vec<f32>>(),
        ecc_filtered.iter().map(|s| 1.0 / (1.0 + s.fit_residual)).collect::<Vec<f32>>(),
    );
    let median_ecc = sigma_clipped_weighted_median(&ecc_vals, &ecc_weights);
    let wmean_ecc = sigma_clipped_weighted_mean(&ecc_vals, &ecc_weights);
```

with:

```rust
    // Eccentricity: fit-residual weighted mean (no ecc pre-filter).
    let ecc_vals: Vec<f32> = measured.iter().map(|s| s.eccentricity).collect();
    let ecc_weights: Vec<f32> = measured.iter()
        .map(|s| 1.0 / (1.0 + s.fit_residual)).collect();
    let median_ecc = sigma_clipped_weighted_mean(&ecc_vals, &ecc_weights);
```

- [ ] **Step 3: Remove `wmean_ecc` from `OurResult` and all usages**

Now that `median_ecc` IS the weighted mean, the separate `wmean_ecc` field is redundant. Remove it from:

1. `OurResult` struct (~line 1331): remove `wmean_ecc: f64,`
2. Early-return constructor (~line 1539): remove `wmean_ecc: 0.0,`
3. Main constructor (~line 1603): remove `wmean_ecc: wmean_ecc as f64,`
4. `wmean_ecc_diffs` computation (~line 1722-1723): remove
5. "Using our wmean" eprintln (~line 1763-1764): remove
6. `our_wmean_eccs` regression (~lines 1786-1787): remove
7. "Ecc regression (wmean)" eprintln (~line 1792-1793): remove
8. TSV header (~line 1928): remove `\tours_wmean_ecc`
9. TSV data row (~line 1932): remove `r.ours.wmean_ecc,`

- [ ] **Step 4: Build debug binary**

Run: `cargo build --features debug-pipeline --release`

Expected: clean build, no warnings

- [ ] **Step 5: Commit**

```text
feat: match ecc weighted mean in debug binary, remove redundant wmean_ecc
```

---

## Task 4: Validate with PI comparison and reference files

- [ ] **Step 1: Run PI comparison**

```bash
target/release/rustafits-debug compare tests/pix_tests/m42_compare.csv
```

Check the aggregate output for:

- Ecc offset should be closer to 0 than before (was -0.073 mean)
- R² should stay ≥ 0.97
- FWHM should be unchanged

- [ ] **Step 2: Run cocoon test**

```bash
target/release/rustafits-debug pipeline tests/cocoon.fits
```

Check that median ecc is higher than before (was 0.426, should be closer to PI's 0.61).

- [ ] **Step 3: Run all four reference files**

```bash
for f in tests/cocoon.fits tests/mono.fits tests/osc.fits tests/test.xisf; do
    echo "=== $f ===" && target/release/rustafits-debug pipeline "$f" 2>&1 | grep -E "Median ecc|Median FWHM"
done
```

- [ ] **Step 4: Run full test suite**

```bash
cargo test
```

Expected: all pass

- [ ] **Step 5: Generate comparison chart**

```bash
target/release/rustafits-debug compare tests/pix_tests/m42_compare.csv 2>/dev/null > /tmp/data.tsv
python3 tests/pix_tests/plot_compare.py /tmp/data.tsv /tmp/ecc_wmean_compare.png
```

---

## Task 5: Update documentation and memory

**Files:**

- Modify: `docs/athenaeum-integration.md` (lines 166, 170-178, 374, 921)
- Modify: `docs/analysis-pipeline.md` (line 132)
- Modify: `docs/trail-rejection.md` (lines 162, 174-176, 179, 234, 249, 265)

- [ ] **Step 1: Update `docs/athenaeum-integration.md`**

Update the accuracy table (line 166) with new R² and offset values from Task 4.

Update the eccentricity methodology note (lines 170-178) to mention weighted mean instead of weighted median.

Update line 374: change "ecc ≤ 0.8 filter is **bypassed** for FWHM, ecc, HFR" to note ecc no longer uses the filter at all.

Update line 921: change "FWHM/ecc/HFR statistics bypass ecc ≤ 0.8 filter on trailed frames" to note eccentricity never uses the filter (only FWHM/HFR use it, and they bypass on trailed frames).

- [ ] **Step 2: Update `docs/analysis-pipeline.md`**

Line 132: change "Trail-aware: bypass ecc<=0.8 filter when trailed" to "Residual-weighted mean (no ecc filter)".

- [ ] **Step 3: Update `docs/trail-rejection.md`**

Lines 162, 174-176: Update the "Effect on Statistics" table. For eccentricity, both normal and trailed frames now use "residual-weighted mean, no ecc filter". FWHM and HFR keep the current behavior (filter on normal, bypass on trailed).

Lines 179, 234, 249, 265: Update any text that says eccentricity uses the ecc filter.

- [ ] **Step 4: Update project memory**

Update `MEMORY.md` line about "Stats ecc filter (ecc <= 0.8)" to reflect that eccentricity no longer uses the filter. Also update the PI comparison results line with new numbers.

- [ ] **Step 5: Commit**

```text
docs: update ecc methodology after weighted mean switch
```
