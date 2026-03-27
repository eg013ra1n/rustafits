# Stage Timing + Batch Analysis API Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add per-stage timing instrumentation to `AnalysisResult` and a batch analysis API for processing multiple images in parallel.

**Architecture:** `StageTiming` struct tracks milliseconds per pipeline stage via `Instant::now()` (always-on, ~20ns overhead). `analyze_batch()` processes paths in chunks of `concurrency` using `par_iter()`, calling `progress(completed, total, path)` after each frame.

**Tech Stack:** Rust, rayon, std::time::Instant

---

## Task 1: Add `StageTiming` struct and field to `AnalysisResult`

**Files:**
- Modify: `src/analysis/mod.rs:97-147` (add struct + field)
- Modify: `src/lib.rs:16` (export `StageTiming`)

- [ ] **Step 1: Add `StageTiming` struct**

Add before `AnalysisResult` (around line 97):

```rust
/// Per-stage timing in milliseconds for the analysis pipeline.
pub struct StageTiming {
    /// Background mesh + noise estimation.
    pub background_ms: f64,
    /// Pass 1 star detection.
    pub detection_pass1_ms: f64,
    /// Free-beta Moffat calibration on bright stars.
    pub calibration_ms: f64,
    /// Pass 2 star detection with refined kernel.
    pub detection_pass2_ms: f64,
    /// PSF measurement on measured stars.
    pub measurement_ms: f64,
    /// Per-star SNR computation.
    pub snr_ms: f64,
    /// Statistics aggregation (medians, trail detection, etc).
    pub statistics_ms: f64,
    /// Total pipeline wall time.
    pub total_ms: f64,
}
```

- [ ] **Step 2: Add field to `AnalysisResult`**

Add to `AnalysisResult` struct after `median_beta`:

```rust
    /// Per-stage timing breakdown.
    pub stage_timing: StageTiming,
```

- [ ] **Step 3: Export from lib.rs**

Change line 16 in `src/lib.rs`:
```rust
pub use analysis::{AnalysisConfig, AnalysisResult, FitMethod, ImageAnalyzer, StageTiming, StarMetrics};
```

- [ ] **Step 4: Fix compilation — update all `AnalysisResult` constructions**

Three places in `run_analysis()` construct `AnalysisResult`:

1. `make_zero_result` closure (~line 548): add `stage_timing: StageTiming { background_ms: 0.0, detection_pass1_ms: 0.0, calibration_ms: 0.0, detection_pass2_ms: 0.0, measurement_ms: 0.0, snr_ms: 0.0, statistics_ms: 0.0, total_ms: 0.0 }`

2. Main result (~line 702): add `stage_timing` field (will be populated in Task 2)

- [ ] **Step 5: Build**

Run: `cargo build`

- [ ] **Step 6: Commit**

```text
feat: add StageTiming struct to AnalysisResult
```

---

## Task 2: Instrument `run_analysis()` with timing

**Files:**
- Modify: `src/analysis/mod.rs` — `run_analysis()` function (lines ~388-712)

- [ ] **Step 1: Add timing instrumentation**

At the start of `run_analysis()`, add:
```rust
let pipeline_start = std::time::Instant::now();
```

Wrap each stage with `Instant::now()` / `.elapsed()`:

```rust
// Stage 1: Background
let t = std::time::Instant::now();
let cell_size = background::auto_cell_size(width, height);
// ... existing background code ...
let background_ms = t.elapsed().as_secs_f64() * 1000.0;

// Stage 2 Pass 1: Detection
let t = std::time::Instant::now();
// ... existing pass 1 code ...
let detection_pass1_ms = t.elapsed().as_secs_f64() * 1000.0;

// Calibration
let t = std::time::Instant::now();
// ... existing calibration code ...
let calibration_ms = t.elapsed().as_secs_f64() * 1000.0;

// Stage 2 Pass 2: Detection
let t = std::time::Instant::now();
// ... existing pass 2 code ...
let detection_pass2_ms = t.elapsed().as_secs_f64() * 1000.0;

// Stage 3: Measurement
let t = std::time::Instant::now();
// ... existing measurement code ...
let measurement_ms = t.elapsed().as_secs_f64() * 1000.0;

// SNR
let t = std::time::Instant::now();
snr::compute_star_snr(...);
let snr_ms = t.elapsed().as_secs_f64() * 1000.0;

// Stage 4: Statistics
let t = std::time::Instant::now();
// ... existing stats code ...
let statistics_ms = t.elapsed().as_secs_f64() * 1000.0;

let total_ms = pipeline_start.elapsed().as_secs_f64() * 1000.0;
```

- [ ] **Step 2: Add `stage_timing` to the result constructions**

In the main `Ok(AnalysisResult { ... })` at ~line 702, add:
```rust
stage_timing: StageTiming {
    background_ms, detection_pass1_ms, calibration_ms,
    detection_pass2_ms, measurement_ms, snr_ms,
    statistics_ms, total_ms,
},
```

Update `make_zero_result` to capture the timing accumulated so far (use `pipeline_start.elapsed()`).

- [ ] **Step 3: Run tests**

Run: `cargo test`

- [ ] **Step 4: Commit**

```text
feat: instrument analysis pipeline with per-stage timing
```

---

## Task 3: Add `analyze_batch()` method

**Files:**
- Modify: `src/analysis/mod.rs` — add method to `ImageAnalyzer` impl block

- [ ] **Step 1: Add the method**

Add to the `impl ImageAnalyzer` block (after `analyze_raw`, before the closing `}`):

```rust
    /// Analyze multiple images in parallel.
    ///
    /// `concurrency` controls how many frames are analyzed simultaneously.
    /// `progress` is called after each frame completes with (completed, total, path).
    /// Returns results in completion order (not necessarily input order).
    pub fn analyze_batch<P, F>(
        &self,
        paths: &[P],
        concurrency: usize,
        progress: F,
    ) -> Vec<(std::path::PathBuf, Result<AnalysisResult>)>
    where
        P: AsRef<std::path::Path> + Sync,
        F: Fn(usize, usize, &std::path::Path) + Send + Sync,
    {
        use std::sync::atomic::{AtomicUsize, Ordering};

        let total = paths.len();
        let completed = AtomicUsize::new(0);
        let concurrency = concurrency.max(1);

        let do_batch = || {
            let mut results = Vec::with_capacity(total);
            for chunk in paths.chunks(concurrency) {
                let chunk_results: Vec<_> = chunk
                    .into_par_iter()
                    .map(|p| {
                        let path = p.as_ref();
                        let result = self.analyze_impl(path);
                        let n = completed.fetch_add(1, Ordering::Relaxed) + 1;
                        progress(n, total, path);
                        (path.to_path_buf(), result)
                    })
                    .collect();
                results.extend(chunk_results);
            }
            results
        };

        match &self.thread_pool {
            Some(pool) => pool.install(do_batch),
            None => do_batch(),
        }
    }
```

- [ ] **Step 2: Run tests**

Run: `cargo test`

- [ ] **Step 3: Commit**

```text
feat: add analyze_batch() for parallel multi-image analysis
```

---

## Task 4: Validate

- [ ] **Step 1: Test stage timing output**

```bash
cargo build --features debug-pipeline --release
target/release/rustafits-debug pipeline tests/cocoon.fits 2>&1 | head -20
```

Add a temporary print in the debug binary to show `result.stage_timing` fields, or add a test that verifies timing fields are non-zero.

- [ ] **Step 2: Test batch API**

Write a quick test or use the debug binary to batch-analyze test files:

```rust
#[test]
fn test_analyze_batch() {
    let analyzer = ImageAnalyzer::new();
    let paths = vec!["tests/cocoon.fits", "tests/mono.fits"];
    let results = analyzer.analyze_batch(&paths, 2, |done, total, path| {
        eprintln!("[{}/{}] {}", done, total, path.display());
    });
    assert_eq!(results.len(), 2);
    for (path, result) in &results {
        assert!(result.is_ok(), "Failed: {}", path.display());
        let r = result.as_ref().unwrap();
        assert!(r.stage_timing.total_ms > 0.0);
    }
}
```

- [ ] **Step 3: Run full test suite**

Run: `cargo test`

- [ ] **Step 4: Commit**

```text
test: validate stage timing and batch analysis
```
