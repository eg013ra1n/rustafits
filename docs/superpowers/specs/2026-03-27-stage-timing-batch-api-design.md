# Per-Stage Timing + Batch Analysis API

## Per-Stage Timing

Add `StageTiming` to `AnalysisResult` — always computed (no config flag, `Instant::now()` is ~20ns per call).

```rust
pub struct StageTiming {
    pub background_ms: f64,
    pub detection_pass1_ms: f64,
    pub calibration_ms: f64,
    pub detection_pass2_ms: f64,
    pub measurement_ms: f64,
    pub snr_ms: f64,
    pub statistics_ms: f64,
    pub total_ms: f64,
}
```

Added to `AnalysisResult`:
```rust
pub struct AnalysisResult {
    // ... existing fields ...
    pub stage_timing: StageTiming,
}
```

Implementation: wrap each stage in `run_analysis()` with `Instant::now()` / `.elapsed()`.

## Batch Analysis API

```rust
impl ImageAnalyzer {
    /// Analyze multiple images in parallel.
    /// `concurrency` controls how many frames run simultaneously.
    /// `progress` is called after each frame completes: (completed, total, path).
    /// Returns results in completion order (not input order).
    pub fn analyze_batch<P, F>(
        &self,
        paths: &[P],
        concurrency: usize,
        progress: F,
    ) -> Vec<(PathBuf, Result<AnalysisResult>)>
    where
        P: AsRef<Path> + Sync,
        F: Fn(usize, usize, &Path) + Send + Sync,
}
```

Implementation: process paths in chunks of `concurrency`, each chunk runs `par_iter().map(analyze).collect()`. Progress callback fires per-frame. Results accumulate in completion order.

```rust
let mut results = Vec::with_capacity(total);
for chunk in paths.chunks(concurrency) {
    let chunk_results: Vec<_> = chunk.par_iter().map(|p| {
        let result = self.analyze(p);
        let n = completed.fetch_add(1, Ordering::Relaxed) + 1;
        progress(n, total, p.as_ref());
        (p.as_ref().to_path_buf(), result)
    }).collect();
    results.extend(chunk_results);
}
```

Uses the analyzer's thread pool if set, otherwise rayon's global pool.

## Files to modify

| File | Change |
|------|--------|
| `src/analysis/mod.rs` | Add `StageTiming` struct, add field to `AnalysisResult`, instrument `run_analysis()`, add `analyze_batch()` |
| `src/lib.rs` | Export `StageTiming` if needed |
