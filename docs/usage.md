# Analysis API Usage Guide

The crate has two analyzers.

- **`star_measure::analysis_mode`** is the measurement behind Athenaeum's Analysis
  tab, and the one `rustafits analyze` and `rustafits convert -a` run. Use it.
- **`ImageAnalyzer::analyze*`** is the legacy analyzer. It is deprecated since 1.8.0
  and removed in 2.0. Its guide is [further down](#legacy-analyzer-deprecated-since-180-removed-in-20).

## Star measurement: `analysis_mode`

Reading a file, measuring it and reading the figures are three calls. Nothing is
configured by a builder: `FrameOptions` has three settings.

```rust
use std::path::Path;
use astroimage::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
use astroimage::star_measure::raw::read_raw_frame;

let mut frame = read_raw_frame(Path::new("light_001.fits"))?;
frame.drop_unmeasured_planes();

let opts = FrameOptions::default();               // 6.0 σ, Moffat 4, 500 stars
let pm = measure_raw_frame(&frame, &opts, None);  // None: rayon's global pool
let s = FrameSummary::from_measurement(&pm, &opts);

println!("Stars detected: {}", s.stars_detected);
println!("Median FWHM:    {:.2} px", s.median_fwhm);
println!("Median HFR:     {:.2} px", s.median_hfr);
println!("Median SNR:     {:.1}", s.median_snr);
println!("Eccentricity:   {:.3}", s.median_eccentricity);
println!("SNR weight:     {:.1}", s.snr_weight);
println!("Frame SNR:      {:.1}", s.frame_snr);
```

`raw` and `analysis_mode` are the stable entry of `star_measure`. `FrameOptions` and
`FrameSummary` are `#[non_exhaustive]`, so 1.x can add fields to them. The rest of
`star_measure` is public but outside that promise.

### 1. Read: `raw::read_raw_frame`

`read_raw_frame(path)` reads a FITS or XISF file as a `RawFrame`:

- `planes`: one plane for a mono file, the mosaic plus its Bayer pattern for a
  one-shot-color file, or three planes for a file that is already debayered. Values
  are the file's own ADU, with no stretch and no debayer.
- `w`, `h`: the size of the stored image. A mosaic keeps its full size.
- `flip_vertical`: the file is bottom-up. The rows are still in file order.
- `domain`: the full scale, the clip from `DATAMAX` and the gain from `EGAIN`.
- `channels`: 1 for mono and mosaic, 3 for a debayered file.

The header facts are read by Athenaeum's own rules, so a malformed value gives the
same number here as in the Analysis tab. A file that cannot be read is an error that
names the file.

`frame.drop_unmeasured_planes()` frees the red and blue planes of a debayered file
(about two thirds of its memory). The measurement uses the green plane only. It does
nothing for mono and mosaic frames.

### 2. Measure: `measure_raw_frame`

`measure_raw_frame(&frame, &opts, pool)` returns a `PlaneMeasurement`: the star list
(`pm.stars`), the sky, the noise and the frame statistics. Mono frames are measured as
they are. A mosaic is measured on its green sites, at its native size, with no
interpolation. A debayered file is measured on its green plane.

| `FrameOptions` field | Default | Meaning |
|----------------------|---------|---------|
| `detection_threshold` | 6.0 | Matched-filter threshold, σ above the noise. Athenaeum accepts 2 to 20 |
| `psf_model` | `PsfModel::Moffat4` | `Moffat4` fixes beta at 4. `Auto` fits the frame's own beta |
| `max_stars` | 500 | The budget: the brightest good stars measured over the field |

The struct is `#[non_exhaustive]`, so set fields on the default:

```rust
use astroimage::star_measure::PsfModel;

let mut opts = FrameOptions::default();
opts.detection_threshold = 4.0;
opts.psf_model = PsfModel::Auto;
opts.max_stars = 1000;
```

The last argument is `Option<&Arc<rayon::ThreadPool>>`. `None` runs on the global pool.
To keep several frames from fighting over the cores, pass one shared pool:

```rust
use std::sync::Arc;
use astroimage::ThreadPoolBuilder;

let pool = Arc::new(ThreadPoolBuilder::new().num_threads(4).build()?);
let pm = measure_raw_frame(&frame, &opts, Some(&pool));
```

### 3. Summarize: `FrameSummary`

`FrameSummary::from_measurement(&pm, &opts)` gives the frame figures:

| Field | Type | Description |
|-------|------|-------------|
| `stars_detected` | `usize` | Stars the matched filter found, extended objects left out |
| `median_fwhm`, `median_hfr` | `f64` | Medians in pixels. Good seeing: FWHM 2 to 3 px |
| `median_eccentricity` | `f64` | 0 = round. Above 0.5 suggests tracking trouble |
| `median_snr` | `f64` | Median per-star SNR |
| `median_beta` | `Option<f64>` | The frame's fitted beta under `PsfModel::Auto`, else `None` |
| `snr_weight` | `f64` | `median(star flux)² / (noise² × background)`, for stacking weights |
| `psf_signal` | `f64` | Median star amplitude over the noise |
| `frame_snr` | `f64` | `background / noise` |
| `background`, `noise` | `f64` | Sky median and noise σ, in ADU |
| `trail_r_squared` | `f64` | Rayleigh statistic of the star angles. 0 = no preferred direction, 1 = all aligned |
| `stars_refused` | `usize` | Stars the fit refused |
| `verdict` | `Option<Verdict>` | `Trailed`, `Defocused` or `NoStars`, else `None` |

The medians, `snr_weight` and `psf_signal` read 0.0 when no fit stands in the frame
statistic. A frame with a `verdict` has no medians worth showing:
`rustafits analyze` prints `-` for them.

### The star list and the overlay

`pm.stars` holds one fit per measured star, in file row order (not flipped). A display
flips them when `frame.flip_vertical` is true. To draw them on a rendered image, build
a `StarOverlay`:

```rust
use astroimage::{annotate_image_with, AnnotationConfig, ImageConverter, StarOverlay};

let summary = FrameSummary::from_measurement(&pm, &opts);
let mut image = ImageConverter::new().process("light_001.fits")?;

let overlay = StarOverlay::from_measurement(frame.w, frame.h, &pm, &summary);
annotate_image_with(&mut image, &overlay, &AnnotationConfig::default());
ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```

`frame.w` and `frame.h` are the measured grid. For a mosaic the render is half that
size, and the annotation functions scale the stars onto it. The overlay carries the
frame's median FWHM too, so `ColorScheme::Fwhm` can color each star by its ratio to
it. See the [Annotation Documentation](annotation.md) for the three tiers
(`compute_overlay`, `create_overlay_layer`, `annotate_image_with`) and for building a
`StarOverlay` from your own star list.

### Many frames

There is no batch call. Run one read, one measurement and one summary per file on your
own workers, and pass the same shared pool to every `measure_raw_frame`. Keep the
`FrameSummary` and the star list of each frame, not the `RawFrame`: a 26 MP frame is
hundreds of MB. On 26 MP frames the CLI peaks near 0.9 GB with one frame at a time and
near 1.5 GB with two, which is why its `-j` default is 2.
`examples/analyze_batch.rs` shows the pattern.

### Fast detection

When you need only `(x, y, flux)` centroids, `ImageAnalyzer::detect_fast*` is faster
and is not deprecated. See the README's "Fast star detection".

---

## Legacy analyzer (deprecated since 1.8.0, removed in 2.0)

`ImageAnalyzer::analyze`, `analyze_data`, `analyze_raw` and `analyze_batch` are
`#[deprecated]`. They keep their old numbers and their old API in 1.x. The numbers are
not the Analysis tab's: a different detector, a different fit and a different
statistic. Everything below describes them.

### Quick Start

```rust
use astroimage::ImageAnalyzer;

let result = ImageAnalyzer::new()
    .analyze("light_001.fits")?;

println!("Stars detected: {}", result.stars_detected);
println!("Median FWHM:    {:.2} px", result.median_fwhm);
println!("Median HFR:     {:.2} px", result.median_hfr);
println!("Median SNR:     {:.1}", result.median_snr);
println!("Eccentricity:   {:.3}", result.median_eccentricity);
println!("Moffat beta:    {:.2}", result.median_beta.unwrap_or_default());
println!("SNR weight:     {:.1}", result.snr_weight);
println!("PSF signal:     {:.1}", result.psf_signal);
println!("Frame SNR:      {:.1}", result.frame_snr);
```

`analyze()` accepts any path to a FITS or XISF file. It handles byte-swap,
u16-to-f32 conversion, and debayering (for OSC/Bayer images) automatically.

### Configuration

`ImageAnalyzer` uses a builder pattern. All settings have sensible defaults;
override only what you need.

```rust
let analyzer = ImageAnalyzer::new()
    .with_detection_sigma(3.0)
    .with_min_star_area(8)
    .with_max_stars(500);

let result = analyzer.analyze("light_001.fits")?;
```

#### Builder Methods

| Method | Default | Description |
|--------|---------|-------------|
| `with_detection_sigma(f32)` | 5.0 | Detection threshold in sigma above background. Lower values (2-3) find fainter stars but increase false positives. Clamped to minimum 1.0. |
| `with_min_star_area(usize)` | 5 | Minimum star area in pixels. Filters hot pixels and noise. Clamped to minimum 1. |
| `with_max_star_area(usize)` | 2000 | Maximum star area in pixels. Filters extended objects like galaxies and nebulae. |
| `with_saturation_fraction(f32)` | 0.95 | Reject stars with peak above this fraction of 65535. Clamped to 0.5-1.0. |
| `with_max_stars(usize)` | 200 | Keep only the brightest N stars. Clamped to minimum 1. |
| `with_mrs_layers(usize)` | 0 | MRS wavelet noise layers. 0 = fast MAD noise (default). Set to 4 for MRS wavelet with B3-spline à trous significance masking (more robust on nebula-rich fields, slower). |
| `without_debayer()` | Debayer enabled | Skip green-channel interpolation for OSC images. By default, OSC images get native-resolution green interpolation. This flag disables it, running analysis directly on the raw Bayer mosaic (faster but less accurate). |
| `with_trail_threshold(f32)` | 0.5 | R² threshold for trail detection (Path A). Lower = more aggressive. The raw `trail_r_squared` is always reported regardless of this setting. Clamped to 0.0-1.0. |
| `with_measure_cap(usize)` | 500 | Max stars to PSF-fit for statistics. Stars are sorted by flux (brightest first) before capping. Set to 0 to measure all (catalog export mode). Dense fields (100K+ stars) benefit most from capping. |
| `with_fit_max_iter(usize)` | 25 | LM max iterations for pass-2 measurement fits. Calibration pass always uses 50. |
| `with_fit_tolerance(f64)` | 1e-4 | LM convergence tolerance for pass-2 measurement fits. Calibration pass always uses 1e-6. |
| `with_fit_max_rejects(usize)` | 5 | Consecutive LM step rejects before early bailout. Stars that bail out fall through to the next method in the Moffat→Gaussian→Moments chain. |
| `with_optics(f64, f64)` | None | `with_optics(focal_mm, pixel_um)` enables arcsecond output for FWHM and HFR. |
| `with_centroid_refine(bool)` | `false` | Refine each fast-detection centroid with the Moffat LM and fill `sx`/`sy`/`fwhm` on `FastStar`. Off by default so `detect_fast` output stays byte-identical; per-star gates (SNR < 10, non-physical fit, centre moved > 2 px) keep the pass-1 centroid. |
| `with_thread_pool(Arc<ThreadPool>)` | Global rayon pool | Route all parallel work through a custom rayon thread pool. |

#### When to Adjust

- **Crowded fields** — lower `detection_sigma` to 3.0, raise `max_stars` to 500+.
- **Sparse fields / short exposures** — lower `detection_sigma` to 3.0-4.0 to pick up faint stars.
- **Hot pixel problems** — raise `min_star_area` to 8-10.
- **Nebula regions** — lower `max_star_area` to avoid detecting bright nebula knots as stars.
- **Speed over accuracy** — use `without_debayer()`.
- **MRS wavelet layers** — default 0 uses fast MAD noise. Set to 4 for MRS wavelet on nebula-rich fields where MAD may overestimate noise.
- **OSC images** — green-channel interpolation and green-pixel-only fitting are applied automatically. No configuration needed.

### Analyzing Raw Data

If you already have pixel data in memory, use `analyze_data()` to skip file I/O:

```rust
use astroimage::ImageAnalyzer;

// Planar f32 data in ADU range [0, 65535]
// For 3-channel: layout is RRRGGGBBB (channel planes stored contiguously)
let data: Vec<f32> = load_my_data();
let width = 4096;
let height = 2048;
let channels = 1; // 1 = mono, 3 = RGB

let result = ImageAnalyzer::new()
    .analyze_data(&data, width, height, channels)?;
```

`analyze_data()` skips file reading and debayering — the data must already be
in planar f32 format with values in the 0-65535 ADU range. Channel count must
be 1 (mono) or 3 (RGB in planar RRRGGGBBB layout).

### Batch Analysis

Use `analyze_batch()` to process multiple files with internal worker orchestration:

```rust
use astroimage::ImageAnalyzer;

let analyzer = ImageAnalyzer::new()
    .with_thread_pool(pool.clone());

let paths = vec!["light_001.fits", "light_002.fits", "light_003.fits"];
let results = analyzer.analyze_batch(
    paths.iter(),
    2,                           // concurrency
    |completed, total| {         // progress callback
        println!("{}/{}", completed, total);
    },
);

for result in results {
    match result {
        Ok(r) => println!("FWHM: {:.2}", r.median_fwhm),
        Err(e) => eprintln!("Error: {}", e),
    }
}
```

### Arcsecond Output (with_optics)

When optical parameters are known, FWHM and HFR can be reported in arcseconds:

```rust
let result = ImageAnalyzer::new()
    .with_optics(1000.0, 3.76)   // focal_length_mm, pixel_size_um
    .analyze("light_001.fits")?;

// median_fwhm and median_hfr are now in arcseconds
println!("FWHM: {:.2}\"", result.median_fwhm);
```

### Reading Results

#### `AnalysisResult`

| Field | Type | Description |
|-------|------|-------------|
| `width` | `usize` | Image width after debayer (if applicable). |
| `height` | `usize` | Image height after debayer (if applicable). |
| `source_channels` | `usize` | 1 = mono, 3 = color (after debayer). |
| `background` | `f32` | Global background level in ADU. |
| `noise` | `f32` | Background noise sigma in ADU. |
| `detection_threshold` | `f32` | Actual detection threshold used (ADU above background). |
| `stars_detected` | `usize` | Total stars detected (raw detection count, before measure cap). |
| `stars` | `Vec<StarMetrics>` | Per-star metrics, sorted by flux descending, capped at `max_stars`. |
| `median_fwhm` | `f32` | Median FWHM across all measured stars (pixels). Good seeing: 2-3 px. |
| `median_eccentricity` | `f32` | Median eccentricity. 0 = perfectly round. Values > 0.5 suggest tracking issues. |
| `median_snr` | `f32` | Median per-star SNR. |
| `median_hfr` | `f32` | Median half-flux radius (pixels). Similar to FWHM but more robust to non-Gaussian profiles. |
| `snr_weight` | `f32` | SNR weight for stacking: `median(star_flux)² / (noise² × background)`. Star-based, gradient-immune. |
| `frame_snr` | `f32` | Per-frame SNR: `background / noise`. Use for stacking prediction: `stacked_snr = sqrt(sum(frame_snr_i²))`. |
| `median_beta` | `Option<f32>` | Median Moffat beta. Always present with default pipeline. Typical: 2-5. |
| `measured_fwhm_kernel` | `f32` | FWHM used for final matched filter kernel in pixels. |
| `psf_signal` | `f32` | PSF signal strength: `median(star_peaks) / noise`. Higher = better signal. |
| `trail_r_squared` | `f32` | Rayleigh R² statistic for directional coherence of star position angles. 0.0 = uniform (no trail), 1.0 = all aligned (strong trail). |
| `possibly_trailed` | `bool` | True if the image is likely trailed (R² above threshold or eccentricity-gated Rayleigh fires). |
| `stage_timing` | `StageTiming` | Per-stage timing in milliseconds (background, detection, calibration, measurement, SNR, etc.). |

#### `StarMetrics`

| Field | Type | Description |
|-------|------|-------------|
| `x`, `y` | `f32` | Subpixel centroid position. |
| `peak` | `f32` | Background-subtracted peak value (ADU). |
| `flux` | `f32` | Total background-subtracted flux (ADU). |
| `fwhm_x` | `f32` | FWHM along major axis (pixels). |
| `fwhm_y` | `f32` | FWHM along minor axis (pixels). |
| `fwhm` | `f32` | Geometric mean FWHM: `sqrt(fwhm_x * fwhm_y)`. |
| `eccentricity` | `f32` | 0 = round, approaching 1 = elongated. |
| `snr` | `f32` | Per-star aperture photometry SNR. |
| `hfr` | `f32` | Half-flux radius (pixels). |
| `theta` | `f32` | PSF position angle in radians, counter-clockwise from +X axis. Orientation of the major axis (`fwhm_x` direction). |
| `beta` | `Option<f32>` | Moffat beta shape parameter. `None` if Gaussian/moments fit was used. |
| `fit_method` | `FitMethod` | Which PSF fitting method produced this measurement: `FreeMoffat`, `FixedMoffat`, `Gaussian`, or `Moments`. |

### Thread Pool Usage

By default, parallel work uses the global rayon thread pool. To control
concurrency (e.g., to limit CPU usage or run multiple analyses concurrently),
pass a custom pool:

```rust
use std::sync::Arc;
use astroimage::{ImageAnalyzer, ThreadPoolBuilder};

let pool = Arc::new(
    ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("Failed to create thread pool")
);

let analyzer = ImageAnalyzer::new()
    .with_thread_pool(pool.clone());

let result = analyzer.analyze("light_001.fits")?;
```

`ThreadPool` and `ThreadPoolBuilder` are re-exported from `rayon` by
`rustafits`, so you don't need to add `rayon` as a direct dependency.

### Star Annotation Overlay

Analysis results can be visualized as color-coded ellipses on converted images.
See the [Annotation Documentation](annotation.md) for the full API reference,
including three tiers of integration (raw geometry, RGBA overlay layer, burn-in)
and configurable color thresholds.

Quick example:

```rust
use astroimage::{ImageConverter, ImageAnalyzer, annotate_image, AnnotationConfig};

let mut image = ImageConverter::new().process("light.fits")?;
let result = ImageAnalyzer::new().analyze("light.fits")?;

annotate_image(&mut image, &result, &AnnotationConfig::default());
ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```
