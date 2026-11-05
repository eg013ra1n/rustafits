# rustafits

High-performance FITS/XISF to JPEG/PNG converter for astronomical images with auto-stretch, Bayer debayering, and SIMD acceleration.

## Features

- **FITS & XISF Support**: Native readers for both formats (no external libraries)
- **Auto-Stretch**: Median-based statistical stretching (STF-compatible midtones transfer)
- **Bayer Debayering**: two paths for two jobs — super-pixel 2x2 block averaging for fast display (halves both axes), and **full-resolution VNG** (8-gradient, variable number of gradients) when the output feeds stacking or pixel-level inspection rather than a screen. All four patterns (RGGB, BGGR, GBRG, GRBG)
- **Preview Mode**: 2x2 binning for fast previews
- **SIMD Optimized**: SSE2/AVX2 (x86_64) and NEON (aarch64) with automatic detection
- **RGBA Output**: Optional RGBA pixel data for canvas/web display
- **In-Memory API**: Get raw pixel data without file I/O — ideal for GUI apps
- **Star Measurement**: the measurement behind Athenaeum's Analysis tab (`star_measure::analysis_mode`) — matched-filter detection, a Moffat fit per star, one beta per frame, FWHM/HFR/eccentricity/SNR medians, background and noise, and trail and defocus verdicts. `rustafits analyze` and `--annotate` run it. The older `ImageAnalyzer::analyze*` is deprecated since 1.8.0
- **Plate-Solving Blocks**: quad pattern matching, RANSAC outlier rejection, gnomonic projection, affine/SIP transform fitting and a WCS you can write to a FITS header — building blocks, no catalogue and no orchestration (see [docs/platesolving.md](docs/platesolving.md))
- **Fast Star Detection**: Lean single-pass detector that skips PSF fitting, SNR photometry, and trail detection — intended for pipelines that only need `(x, y, flux)` centroids (blind plate solving, quick previews). Runs in ~300–500 ms on a full-frame image vs. seconds for the legacy analyzer
- **Pure-Rust JPEG Encoder**: SIMD-accelerated JPEG encoding (NEON on aarch64, AVX2/SSE2 on x86_64, scalar fallback elsewhere) — baseline 4:2:0, no C toolchain, no cmake/nasm
- **Star Annotation**: Color-coded ellipse overlay showing PSF shape, elongation direction, and quality grading, drawn from the same measurement `analyze` prints (or from any star list you build)

## Supported Formats

| Format | Extensions | Data Types |
|--------|-----------|------------|
| FITS | `.fits`, `.fit` | 8/16/32-bit int, 32/64-bit float |
| XISF | `.xisf` | All sample formats, zlib/LZ4/Zstd compression |

## Installation

### Cargo (Recommended)

```bash
cargo install rustafits
```

**Build requirements:** none beyond the Rust toolchain. Every dependency is pure
Rust, so `cargo install` and cross-compilation work on x86_64 and aarch64 without
cmake, nasm, or a C compiler.

### From Source

```bash
git clone https://github.com/eg013ra1n/rustafits
cd rustafits
cargo build --release
sudo cp target/release/rustafits /usr/local/bin/
```

### Homebrew (macOS/Linux)

```bash
brew tap eg013ra1n/rustafits
brew install rustafits
```

## CLI Usage

```bash
# One file — the classic form still works
rustafits image.fits output.jpg
rustafits image.xisf output.png

# Output next to the input (light.jpg), or PNG with -f
rustafits light.fits
rustafits light.fits -f png

# Batch: many files or directories into one output directory
rustafits *.fits -o previews/
rustafits night1/ -r -o out/ -n      # recurse, mirror the tree, skip existing

# Pipe: JPEG on stdout
rustafits m31.xisf -o - > m31.jpg

# Fast preview, downscale, annotation, native-resolution OSC
rustafits large.fits -p                       # 2x2 binning
rustafits large.fits -d 4
rustafits image.fits annotated.jpg -a -v
rustafits osc.fits full.jpg --vng
```

Inputs are `.fits` `.fit` `.fts` `.xisf` files or directories. With exactly
two paths, the second is the output if it is `-`, ends in `.jpg`/`.jpeg`/`.png`,
or does not exist yet. A batch shows a progress bar and a summary on stderr; a
failing file is reported and the others carry on (exit code 1 if any failed).
Output names are checked before any work starts — an unsupported extension,
two inputs mapping to one output, or `-` with several inputs is an error.
Files are written through a temporary file and renamed, so an interrupted run
leaves no truncated images behind.

| Option | Default | |
|---|---|---|
| `-o, --output <PATH>` | next to input | Output file (one input; `-` = stdout) or directory (several inputs, existing dir, or trailing `/`) |
| `-f, --format jpg\|png` | jpg | Format for derived names and stdout; an explicit file name decides its own |
| `-r, --recursive` | off | Descend into subdirectories (mirrored under `-o DIR`) |
| `-n, --skip-existing` | off | Skip inputs whose output already exists |
| `-j, --jobs <N>` | min(4, cores); 2 with `-a` | Files converted concurrently |
| `--threads <N>` | all cores | Threads for pixel processing, shared by all jobs |
| `-q, --quality <Q>` | 95 | JPEG quality 1–100 |
| `--png-compression <L>` | fast | `none` \| `fast` \| `balanced` \| `high`. Astro frames barely compress at fast; none is ~2.5x faster, balanced ~2.5x smaller but ~30x slower |
| `-d, --downscale <N>` | 1 | Integer factor. Bayer/OSC: the super-pixel debayer already halves, so 2 = debayer only, 4 = debayer + 2x |
| `-p, --preview` | off | 2x2 binning for mono images |
| `--no-debayer` | off | Keep the raw Bayer mosaic |
| `--vng` | off | Debayer at native resolution (VNG) instead of 2x2 super-pixel |
| `-a, --annotate` | off | Overlay ellipses on the stars `analyze` measures (same numbers, same stars) |
| `--max-stars <N>` | 500 | Star budget for `-a`: the brightest good stars measured and drawn, 10–2000 |
| `-v, --verbose` | off | Per-file details and timings on stderr (`--log` is an alias) |
| `--quiet` | off | Errors only |

Batch throughput: 9 real 26 MP frames on 4 cores take 0.9 s in one call,
versus 2.5 s for a shell loop over single-file calls.

`convert` is the default command; `rustafits convert …` spells it out (for an
input file literally named `info`, `analyze` or `completions`).

### Headers: `rustafits info`

Reads only the header (a few KB per file), so it is instant on a whole night:

```bash
$ rustafits info tests/
FILE             SIZE       TYPE  BAYER    EXP  FILTER  GAIN   TEMP  DATE-OBS                 OBJECT
tests/mono.fits  6252x4176  u16   -       30.0  L         56   -9.5  2025-10-25T16:22:10.582  C/2025 A6 (Lemmon)
tests/osc.fits   6248x4176  u16   RGGB    15.0  -        100  -10.4  2025-10-25T15:53:47.995  C/2025 A6 (Lemmon)
tests/test.xisf  6224x4168  u16   RGGB   120.0  -        101  -10.0  2025-09-26T17:25:21.185  Ghost Nebula

rustafits info night1/ -r -k FOCALLEN -k AIRMASS   # extra keyword columns
rustafits info light.fits --all                    # every keyword
rustafits info night1/ --json | --csv
```

Columns come from the FITS cards N.I.N.A., SharpCap, ASIAIR etc. write, with
XISF property ids as fallback (`Instrument:ExposureTime`, …).

### Star metrics: `rustafits analyze`

Per-frame quality for culling subs. `analyze` measures a frame the way
Athenaeum's Analysis tab does (`star_measure::analysis_mode`): a matched filter
finds the stars at 6.0 σ, each gets a Moffat fit with beta fixed at 4, and the
brightest 500 good ones give the frame's medians. A file reads the same here as
in the tab. Only the flags you give change that.

```bash
$ rustafits analyze tests/mono.fits tests/osc.fits tests/test.xisf --sort fwhm
FILE             STARS  FWHM   HFR    ECC   SNR  WEIGHT    BG  NOISE  VERDICT
tests/mono.fits   4775  2.01  1.18  0.476  11.9   6.669   850   33.5  -
tests/osc.fits     590  2.50  2.49  0.424   9.7   6.170  1272   61.8  -
tests/test.xisf  24284  2.81  2.97  0.590  20.2   190.7   291   11.7  -

rustafits analyze night1/ --csv > night1.csv              # full precision
rustafits analyze night1/ --focal-length 448 --pixel-size 3.76   # + arcsec
rustafits analyze night1/ --detection-threshold 4 --psf-model auto
```

What the columns mean:

- `STARS` is the number of stars the detector found. The medians come from the
  brightest good stars among them, at most `--max-stars`.
- `FWHM` and `HFR` are in pixels (with `--focal-length` and `--pixel-size`, two
  more columns give arcseconds). `ECC` is 0 for round stars. `SNR` is the
  median star SNR, `WEIGHT` the SNR weight for stacking, `BG` and `NOISE` are
  in ADU.
- `VERDICT` is `-` for a usable frame. When the measurement judges a frame
  `trailed`, `defocused` or `no stars`, the word shows here, the medians show
  `-`, and `--sort` by a metric puts the frame after all the others.

Options:

| Option | Default | |
|---|---|---|
| `--detection-threshold <SIGMA>` | 6.0 | Detection threshold in σ above the noise, 2–20 |
| `--psf-model moffat4\|auto` | moffat4 | `moffat4` fixes beta at 4; `auto` fits the frame's own beta (the `median_beta` column) |
| `--max-stars <N>` | 500 | Star budget: the brightest good stars measured per frame, 10–2000 |
| `--focal-length <MM>`, `--pixel-size <UM>` | none | Both together add FWHM and HFR in arcseconds |
| `-s, --sort stars\|fwhm\|hfr\|ecc\|snr\|weight\|name` | input order | Best first; `--reverse` for worst first |
| `--json`, `--csv` | table | One object or row per frame, full precision |
| `-j, --jobs <N>` | min(2, cores, files) | Frames analyzed at once. Analysis already uses every core |
| `--threads <N>` | all cores | Threads for pixel work, shared by all jobs |
| `-r, --recursive`, `--quiet` | off | Descend into directories; errors only |

The CSV and JSON have 25 columns. The first 19 are the ones `analyze` had
before 1.8.0, in the same order, so a script that reads by position keeps
working. Six follow: `channels`, `stars_refused`, `frame_verdict`
(`trailed`, `defocused`, `no_stars` or empty), `measure_version`,
`detection_threshold` and `psf_model`. `median_beta` is empty unless you ask for
`--psf-model auto`.

A file the measurement cannot read is reported on stderr with its path, the
other files carry on, and the exit code is 1.

**Deprecated flags.** Scripts written for 1.7 still run. These flags are hidden
from `--help`, print a warning on stderr and go in 2.0: `--sigma` (now
`--detection-threshold`, for another detector: the default is 6.0, it was 5.0),
`--measure-cap` (no effect: the budget is `--max-stars`) and `--no-debayer`
(no effect: a mosaic is measured on its green sites).

### Shell completion

```bash
rustafits completions zsh  > ~/.zfunc/_rustafits
rustafits completions bash > ~/.local/share/bash-completion/completions/rustafits
rustafits completions fish > ~/.config/fish/completions/rustafits.fish
```

## Library Usage

Add to your `Cargo.toml`:

```toml
[dependencies]
rustafits = "1.8"
# Library only, without the CLI's clap/indicatif dependencies:
# rustafits = { version = "1.8", default-features = false }
```

### File output

```rust
use astroimage::ImageConverter;

ImageConverter::new()
    .with_preview_mode()
    .with_quality(90)
    .convert("input.fits", "output.jpg")?;
```

### In-memory processing

Get raw RGB pixel data without writing to disk — useful for GUI viewers, web backends, and Tauri apps:

```rust
use astroimage::{ImageConverter, ProcessedImage};

let image: ProcessedImage = ImageConverter::new()
    .with_downscale(2)
    .process("input.fits")?;

// image.data     - Vec<u8>, interleaved RGB or RGBA bytes
// image.width    - pixel width
// image.height   - pixel height
// image.channels - 3 (RGB) or 4 (RGBA)
// image.is_color - true if debayered/RGB, false if mono (gray replicated to RGB)
```

### Full-resolution debayer (VNG)

`ImageConverter` debayers with the super-pixel path, which averages each 2x2 CFA
block and so returns half the width and half the height. When you need the native
pixel grid — stacking input, or inspecting stars at 1:1 — call the VNG demosaic
directly:

```rust
use astroimage::processing::vng::vng_debayer_f32;
use astroimage::BayerPattern;

// `mosaic` is one f32 sample per pixel, width * height long.
let rgb: Vec<f32> = vng_debayer_f32(&mosaic, width, height, BayerPattern::Rggb);

// Planar RGB at the INPUT geometry: [R plane][G plane][B plane],
// each width * height long — 3 * width * height in total.
```

It is the more expensive of the two by a wide margin: on a 26 MP (6248x4176)
frame, super-pixel takes ~3.6 ms and returns 6 MP, while VNG takes ~0.7 s and
returns 26 MP — roughly 48x the cost per output pixel, for 4x the pixels, plus
~313 MB for the planar f32 result. Use it deliberately, not as the default for
previews or thumbnails.

The interior pass is parallelised over 64-row bands via rayon. The split is by
output row and every band reads the same immutable mosaic, so the result does not
depend on the thread count or on where a band boundary falls.

### Image analysis

Measure a frame the way Athenaeum's Analysis tab does. Three calls: read the
file as a raw frame, measure it, summarise the result.

```rust
use std::path::Path;
use astroimage::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
use astroimage::star_measure::raw::read_raw_frame;

// 1. Read: the file's own ADU values, in file row order. A Bayer mosaic stays a mosaic.
let mut frame = read_raw_frame(Path::new("light.fits"))?;
frame.drop_unmeasured_planes();   // frees the R and B planes of a debayered file

// 2. Measure. The defaults are the Analysis tab's: 6.0 σ, Moffat 4, 500 stars.
let opts = FrameOptions::default();
let pm = measure_raw_frame(&frame, &opts, None);   // None: rayon's global pool

// 3. Summarise: the frame figures the Analysis tab shows.
let s = FrameSummary::from_measurement(&pm, &opts);
println!("Stars: {}  FWHM: {:.2} px  HFR: {:.2} px  Ecc: {:.3}  SNR: {:.1}",
    s.stars_detected, s.median_fwhm, s.median_hfr, s.median_eccentricity, s.median_snr);
if let Some(v) = s.verdict {
    println!("Verdict: {v:?}");   // Trailed, Defocused or NoStars: the medians read 0.0
}
```

`FrameOptions` has three settings. It is `#[non_exhaustive]`, so you set the
fields on the default:

```rust
use astroimage::star_measure::PsfModel;

let mut opts = FrameOptions::default();
opts.detection_threshold = 4.0;      // σ above the noise; Athenaeum accepts 2 to 20
opts.psf_model = PsfModel::Auto;     // fit the frame's own beta instead of fixing it at 4
opts.max_stars = 1000;               // the brightest good stars measured
```

`FrameSummary` holds the frame figures. Sizes are in pixels, `background` and
`noise` in ADU:

| Field | Meaning |
|-------|---------|
| `stars_detected` | Stars the matched filter found (extended objects left out) |
| `median_fwhm`, `median_hfr`, `median_eccentricity`, `median_snr` | Medians over the measured stars; 0.0 when no fit stands in the statistic |
| `median_beta` | The frame's fitted beta under `PsfModel::Auto`, else `None` |
| `snr_weight`, `psf_signal`, `frame_snr` | SNR weight for stacking, median amplitude over noise, background over noise |
| `background`, `noise` | Sky level and noise σ |
| `trail_r_squared` | Rayleigh statistic of the star angles: 0 = no direction, 1 = all aligned |
| `stars_refused` | Stars the fit refused (the reasons are in the measurement) |
| `verdict` | `Some(Trailed \| Defocused \| NoStars)` or `None` |

The measurement `pm` also holds the star list (`pm.stars`), in file row order.
`frame.flip_vertical` says whether a display must flip it. A Bayer mosaic is
measured at its native size, on its green sites. `raw` and `analysis_mode` are
the stable entry of `star_measure`: `FrameOptions` and `FrameSummary` are
`#[non_exhaustive]`, so 1.x can add fields to them. The rest of `star_measure`
is public but outside that promise.

For many frames, run one measurement per file on your own workers and pass one
shared pool to each (`Some(&pool)`). `examples/analyze_batch.rs` does it.

### Fast star detection

When you only need rough `(x, y, flux)` centroids — e.g. for blind plate solving, quick
previews, or any quad-hash matching pipeline — `detect_fast` skips everything the legacy
analyzer does after the first detection pass. No PSF calibration, no Levenberg–Marquardt
Moffat fitting, no SNR photometry, no trail detection. On a full-frame (~26 Mpx) image this
drops star detection from ~5–6 s to ~300–500 ms in release builds. Centroid accuracy is
pass-1 (intensity-weighted, typically ~0.3–0.5 px) — fine for quad matching, insufficient
for stacking registration.

```rust
use astroimage::ImageAnalyzer;

let result = ImageAnalyzer::new()
    .with_detection_sigma(5.0)
    .with_max_stars(500)
    .detect_fast("light.fits")?;

println!("{} stars in {:.0} ms ({}x{})",
    result.stars.len(),
    result.timing.total_ms,
    result.width,
    result.height);

for s in result.stars.iter().take(10) {
    println!("  ({:.2}, {:.2}) flux={:.0}", s.x, s.y, s.flux);
}
```

All three entry-point variants mirror the legacy `analyze*` entry points:

| Method | Input | Use case |
|--------|-------|----------|
| `detect_fast(path)` | FITS / XISF file | Most common entry point |
| `detect_fast_data(&[f32], w, h, channels)` | Pre-loaded planar f32 | You already decoded pixels elsewhere |
| `detect_fast_raw(&ImageMetadata, &PixelData)` | Read-but-not-converted | Skip file I/O while keeping u16→f32 + debayer internal |

All three honor the same `ImageAnalyzer` builder settings as the legacy `analyze` —
`with_detection_sigma`, `with_max_stars`, `with_min_star_area`, `with_max_star_area`,
`with_saturation_fraction`, `without_debayer`, and `with_thread_pool`. Settings related
to the PSF fit (`with_measure_cap`, `with_fit_max_iter`, `with_fit_tolerance`,
`with_mrs_layers`, `with_trail_threshold`, `with_optics`) are simply ignored because
those stages don't run.

Pipeline stages in the fast path:

1. OSC green interpolation (for Bayer frames, same as the legacy `analyze`)
2. Luminance extraction if `channels == 3`
3. `estimate_background_mesh` (parallelized, MAD noise — MRS wavelet is never called)
4. Single matched-filter pass of `detect_stars` with fixed FWHM = 3.0 px
5. Sort by flux descending, truncate to `max_stars`, pack into `FastStar`

**When to use which**: if you need FWHM, eccentricity, SNR, HFR, per-star fit quality,
or trail detection, use `star_measure::analysis_mode` (see Image analysis). If you only
need positions and brightness ordering (quad matching, hash lookup, quick counts), use
`detect_fast`: it skips the PSF fits, SNR and verdicts and returns first-pass centroids.
On a 26 MP frame the two take comparable wall time (see Performance), so choose by what
you need, not by speed.

### Legacy analyzer (deprecated since 1.8.0, removed in 2.0)

`ImageAnalyzer::analyze`, `analyze_data`, `analyze_raw` and `analyze_batch` are
deprecated. The compiler says so and names `star_measure::analysis_mode` as the
replacement. They keep their old numbers in 1.x and go in 2.0. `detect_fast*`
and the builder methods are not deprecated (see
[Fast star detection](#fast-star-detection)).

#### Single frame

`analyze` detects stars, measures PSF shape and computes image quality metrics:

```rust
use astroimage::ImageAnalyzer;

let result = ImageAnalyzer::new()
    .with_max_stars(500)
    .with_optics(620.0, 3.76)  // focal length mm, pixel size µm → arcsec output
    .analyze("light.fits")?;

println!("Stars: {}  FWHM: {:.2} px ({:.1}\")  Ecc: {:.3}  Seeing: {:.1}\"",
    result.stars_detected, result.median_fwhm,
    result.median_fwhm_arcsec.unwrap_or(0.0),
    result.median_eccentricity,
    result.median_fwhm_arcsec.unwrap_or(0.0));

// Per-stage timing breakdown
let t = &result.stage_timing;
println!("Timing: bg={:.0}ms det={:.0}ms cal={:.0}ms meas={:.0}ms total={:.0}ms",
    t.background_ms, t.detection_pass1_ms, t.calibration_ms,
    t.measurement_ms, t.total_ms);
```

Default configuration uses a two-pass calibration pipeline: pass 1 fits free-beta Moffat
on bright calibration stars to derive the field PSF model (beta, FWHM). Pass 2 applies
fixed-beta Moffat to all detected stars with Gaussian and moments fallbacks. Background
uses parallelized mesh-grid with auto-tuned cell size and MAD noise estimation (MRS wavelet
available via `with_mrs_layers(4)`). OSC/Bayer images are green-interpolated before
detection and PSF fitting.

#### Batch

`analyze_batch` analyzes multiple images in parallel with progress reporting:

```rust
use astroimage::ImageAnalyzer;

let analyzer = ImageAnalyzer::new()
    .with_optics(620.0, 3.76);

let paths: Vec<&str> = vec!["frame001.fits", "frame002.fits", /* ... */];

let results = analyzer.analyze_batch(&paths, 4, |done, total, path| {
    println!("[{}/{}] {}", done, total, path.display());
});

for (path, result) in &results {
    match result {
        Ok(r) => println!("{}: FWHM={:.2}\" ecc={:.3}",
            path.display(),
            r.median_fwhm_arcsec.unwrap_or(0.0),
            r.median_eccentricity),
        Err(e) => eprintln!("{}: {}", path.display(), e),
    }
}
```

The `concurrency` parameter controls how many frames are analyzed simultaneously.
Results are returned in approximate completion order with their paths.

### Star annotation overlay

Draw color-coded ellipses showing PSF shape and quality. The stars come from a
`StarOverlay`: build one from an `analysis_mode` measurement, as the CLI does for
`--annotate`.

```rust
use std::path::Path;
use astroimage::{annotate_image_with, AnnotationConfig, ColorScheme, ImageConverter, StarOverlay};
use astroimage::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary};
use astroimage::star_measure::raw::read_raw_frame;

// Measure (see Image analysis), then render the same file.
let mut frame = read_raw_frame(Path::new("light.fits"))?;
frame.drop_unmeasured_planes();
let opts = FrameOptions::default();
let pm = measure_raw_frame(&frame, &opts, None);
let summary = FrameSummary::from_measurement(&pm, &opts);
let mut image = ImageConverter::new().process("light.fits")?;

// The star list, in the pixels the frame was measured on. A mosaic is measured at
// its native size and rendered at half size; the overlay scales for you.
let overlay = StarOverlay::from_measurement(frame.w, frame.h, &pm, &summary);

// Burn annotations with default settings (eccentricity color coding)
annotate_image_with(&mut image, &overlay, &AnnotationConfig::default());

// Or customize thresholds and color scheme instead
let config = AnnotationConfig {
    color_scheme: ColorScheme::Eccentricity,
    ecc_good: 0.5,   // ≤ 0.5 → green
    ecc_warn: 0.6,   // 0.51–0.6 → yellow, > 0.6 → red
    ..AnnotationConfig::default()
};
annotate_image_with(&mut image, &overlay, &config);

ImageConverter::save_processed(&image, "annotated.jpg", 95)?;
```

A `StarOverlay` is a plain struct, so you can fill one from your own star list:

| Field | Meaning |
|-------|---------|
| `stars` | `Vec<OverlayStar>`: `x`, `y` (source pixels, file row order), `fwhm_x` (major axis), `fwhm_y` (minor axis), `theta` (major axis angle, radians, counter-clockwise from +X), `eccentricity`, `fwhm` |
| `width`, `height` | Size of the pixel grid the stars were measured on |
| `median_fwhm` | The frame's median FWHM, for `ColorScheme::Fwhm` (0.0 colors every star green) |

Three API tiers for different integration needs:

| Function | Returns | Use Case |
|----------|---------|----------|
| `compute_overlay()` | `Vec<StarAnnotation>` | Raw geometry for custom rendering (Canvas2D, SwiftUI, SVG) |
| `create_overlay_layer()` | `Vec<u8>` (RGBA) | Transparent overlay for toggleable layer compositing |
| `annotate_image_with()` | modifies `ProcessedImage` | Burn-in for CLI or one-shot use |

Each has a legacy twin that takes an `AnalysisResult` from the deprecated
`ImageAnalyzer::analyze*`: `compute_annotations()`, `create_annotation_layer()` and
`annotate_image()`. They build `StarOverlay::from(&result)` and call the function above, so
both draw the same.

**`compute_overlay(overlay, width, height, flip_vertical, config)`** — Transforms star positions from the measured grid to output image coordinates (handling debayer scaling, downscale, and vertical flip), computes ellipse semi-axes from `fwhm_x`/`fwhm_y`, and assigns colors. Returns `Vec<StarAnnotation>` where each entry contains `x`, `y`, `semi_major`, `semi_minor`, `theta`, `eccentricity`, `fwhm`, and `color` — everything needed to draw the ellipse in any rendering system.

**`create_overlay_layer(overlay, width, height, flip_vertical, config)`** — Calls `compute_overlay()` internally, then rasterizes all ellipses and direction ticks onto a transparent RGBA buffer (same dimensions as the output image). Use as a compositable layer that can be toggled on/off without re-rendering the base image.

**`annotate_image_with(image, overlay, config)`** — Calls `compute_overlay()` internally, then draws directly onto the `ProcessedImage.data` buffer (RGB or RGBA). Reads `image.flip_vertical` automatically. Simplest path — one call, image modified in place.

**`ImageConverter::save_processed(image, path, quality)`** — Saves a `ProcessedImage` to disk as JPEG or PNG. Use after `annotate_image_with()` or any other post-processing on the pixel buffer. The format follows the extension: `.png` → PNG; `.jpg`, `.jpeg` or none → JPEG; anything else is an error. `save_processed_with(image, path, quality, PngCompression::None)` picks the PNG compression level (`ImageConverter::with_png_compression` does the same for `convert()`).

#### AnnotationConfig fields

| Field | Default | Description |
|-------|---------|-------------|
| `color_scheme` | `Eccentricity` | `Eccentricity` (tracking/optics), `Fwhm` (focus), or `Uniform` (all green) |
| `show_direction_tick` | `true` | Draw ticks along elongation axis (visible when ecc > 0.15) |
| `ellipse_scale` | `1.2` | Semi-axis as a multiple of FWHM. Was a hardcoded 2.5, which drew lassos wide enough to cover neighbouring stars on oversampled frames |
| `min_radius` | `6.0` | Minimum ellipse semi-axis in output pixels |
| `max_radius` | `60.0` | Maximum ellipse semi-axis in output pixels |
| `line_width` | `2` | Line thickness: `1` = 1px, `2` = 3px cross, `3` = 5px diamond |
| `ecc_good` | `0.5` | Eccentricity at or below this is green (good) |
| `ecc_warn` | `0.6` | Eccentricity between good and warn is yellow; above is red |
| `fwhm_good` | `1.3` | FWHM ratio (star/median) below this is green |
| `fwhm_warn` | `2.0` | FWHM ratio between good and warn is yellow; above is red |

See [Annotation Documentation](docs/annotation.md) for full API reference, integration examples, and coordinate transform details.

### ImageConverter builder methods

| Method | Description |
|--------|-------------|
| `with_downscale(n)` | Downscale by factor n (Bayer images: debayer counts as 2x, extra downscale applied for n > 2) |
| `with_quality(q)` | JPEG quality 1-100 |
| `without_debayer()` | Skip Bayer debayering |
| `with_preview_mode()` | 2x2 binning for fast previews |
| `with_rgba_output()` | Output RGBA instead of RGB (adds alpha=255 channel) |
| `with_thread_pool(pool)` | Use a custom rayon thread pool (see below) |

### ImageAnalyzer builder methods

The deprecated `analyze*` methods read all of these. `detect_fast*` reads the detection ones (it ignores the PSF-fit settings).

| Method | Description |
|--------|-------------|
| `with_detection_sigma(f32)` | Detection threshold in sigma above background (default 5.0) |
| `with_min_star_area(usize)` | Minimum star area in stamp (default 5 px) |
| `with_max_star_area(usize)` | Maximum star area in stamp (default 2000 px) |
| `with_saturation_fraction(f32)` | Reject stars above this fraction of 65535 (default 0.95) |
| `with_max_stars(usize)` | Keep only the brightest N stars (default 200) |
| `with_measure_cap(usize)` | Max stars to PSF-fit for statistics (default 500, 0 = all) |
| `with_mrs_layers(usize)` | Noise layers: 0 = fast MAD (default), 1-6 = MRS wavelet |
| `with_trail_threshold(f32)` | R² threshold for Rayleigh trail detection (default 0.5) |
| `with_optics(f64, f64)` | Focal length (mm) + pixel size (µm) → enables arcsec output |
| `without_debayer()` | Skip green-channel interpolation for OSC images |
| `with_centroid_refine(bool)` | Run the Moffat LM per detection to refine the pass-1 centroid and fill `sx`/`sy`/`fwhm` (default off, so `detect_fast` output is unchanged unless asked) |
| `with_fit_max_iter(usize)` | LM iteration ceiling per star |
| `with_fit_tolerance(f64)` | LM convergence tolerance |
| `with_fit_max_rejects(usize)` | Max outlier-rejection rounds in the fit |
| `with_thread_pool(pool)` | Use a custom rayon thread pool |

### AnalysisResult fields (legacy)

Returned by the deprecated `ImageAnalyzer::analyze*`.

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `width`, `height` | usize | pixels | Image dimensions |
| `background` | f32 | ADU | Global background level |
| `noise` | f32 | ADU | Background noise sigma |
| `stars_detected` | usize | — | Total detections (before measure cap) |
| `stars` | Vec\<StarMetrics\> | — | Per-star metrics (brightest N) |
| `median_fwhm` | f32 | pixels | Median FWHM across measured stars |
| `median_fwhm_arcsec` | Option\<f32\> | arcsec | Median FWHM (requires `with_optics`) |
| `median_eccentricity` | f32 | — | 0 = round, →1 = elongated |
| `median_hfr` | f32 | pixels | Median half-flux radius |
| `median_hfr_arcsec` | Option\<f32\> | arcsec | Median HFR (requires `with_optics`) |
| `median_snr` | f32 | — | Median per-star SNR |
| `plate_scale` | Option\<f32\> | arcsec/px | Plate scale (requires `with_optics`) |
| `trail_r_squared` | f32 | — | Rayleigh R̄² for directional coherence |
| `possibly_trailed` | bool | — | True if coherent trailing detected |
| `median_beta` | Option\<f32\> | — | Moffat β (None if Gaussian/moments) |
| `pass1_detections` | usize | — | Pass 1 detection count (before calibration) |
| `calibrated_fwhm` | f32 | pixels | Calibrated field FWHM from Moffat pass |
| `stars_measured` | usize | — | Stars that survived PSF fitting |
| `moffat_count` | usize | — | Moffat fits among measured stars |
| `gaussian_count` | usize | — | Gaussian fits among measured stars |
| `stage_timing` | StageTiming | ms | Per-stage timing breakdown |

### StarMetrics fields (legacy)

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `x`, `y` | f32 | pixels | Subpixel centroid position |
| `fwhm_x`, `fwhm_y` | f32 | pixels | FWHM along major/minor axis |
| `fwhm` | f32 | pixels | Geometric mean FWHM |
| `fwhm_arcsec` | Option\<f32\> | arcsec | FWHM (requires `with_optics`) |
| `eccentricity` | f32 | — | 0 = round, →1 = elongated |
| `theta` | f32 | radians | Position angle of major axis |
| `hfr` | f32 | pixels | Half-flux radius |
| `hfr_arcsec` | Option\<f32\> | arcsec | HFR (requires `with_optics`) |
| `snr` | f32 | — | Per-star aperture photometry SNR |
| `peak`, `flux` | f32 | ADU | Peak and total flux |
| `beta` | Option\<f32\> | — | Moffat β parameter |
| `fit_method` | FitMethod | — | FreeMoffat, FixedMoffat, Gaussian, or Moments |
| `fit_residual` | f32 | — | Normalized fit quality (lower = better) |

### StageTiming fields (legacy)

| Field | Type | Description |
|-------|------|-------------|
| `background_ms` | f64 | Background mesh + noise estimation |
| `detection_pass1_ms` | f64 | Pass 1 star detection |
| `calibration_ms` | f64 | Free-beta Moffat calibration |
| `detection_pass2_ms` | f64 | Pass 2 detection with refined kernel |
| `measurement_ms` | f64 | PSF measurement on measured stars |
| `snr_ms` | f64 | Per-star SNR computation |
| `statistics_ms` | f64 | Statistics aggregation |
| `total_ms` | f64 | Total pipeline wall time |

### FastAnalysisResult fields

Returned by `ImageAnalyzer::detect_fast*`. Lean result type — only what blind plate
solving and similar pipelines need.

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `width`, `height` | usize | pixels | Image dimensions (after debayer if applicable) |
| `stars` | Vec\<FastStar\> | — | Stars sorted by flux descending, capped at `max_stars` |
| `background` | f32 | ADU | Global background level from the mesh estimator |
| `noise` | f32 | ADU | Background noise sigma (MAD, no MRS wavelet) |
| `timing` | FastDetectTiming | ms | Per-stage timing breakdown |

### FastStar fields

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `x`, `y` | f32 | pixels | Subpixel centroid — pass-1 intensity-weighted, or PSF-refined when `with_centroid_refine(true)` |
| `raw_x`, `raw_y` | f32 | pixels | Always the unrefined pass-1 centroid; equals `x`/`y` when refinement is off or the per-star fit was rejected |
| `peak` | f32 | ADU | Background-subtracted peak value |
| `flux` | f32 | ADU | Background-subtracted total flux |
| `snr` | f32 | — | Aperture-photometry SNR, `flux / sqrt(flux + pi r^2 sigma^2)`. Separates compact sources from extended structure that carries high flux at low SNR |
| `sx`, `sy` | f32 | pixels | Per-axis Gaussian sigma from the Moffat fit (`FWHM / 2.3548`). `0.0` when refinement is off |
| `fwhm` | f32 | pixels | Mean fitted FWHM, `0.5 * (FWHM_x + FWHM_y)`. `0.0` when refinement is off |

### FastDetectTiming fields

| Field | Type | Description |
|-------|------|-------------|
| `read_ms` | f64 | File I/O (0.0 for `detect_fast_data` / `detect_fast_raw`) |
| `prep_ms` | f64 | f32 conversion + OSC green interpolation + luminance extraction |
| `background_ms` | f64 | Mesh-grid background and noise estimation |
| `detection_ms` | f64 | Single-pass matched-filter detection |
| `total_ms` | f64 | Wall clock from entry to return |

### Multi-image concurrent processing

By default, all parallel work (debayering, stretch, binning, byte conversion) runs on rayon's global thread pool. This works well for single-image processing, but when processing multiple images concurrently from separate threads, they all compete for the same pool — causing thread oversubscription and degraded throughput.

Use `with_thread_pool()` to route all parallel work to a dedicated or shared pool:

```rust
use std::sync::Arc;
use astroimage::{ImageConverter, ThreadPoolBuilder};

// Create a shared pool once at startup
let pool = Arc::new(
    ThreadPoolBuilder::new()
        .num_threads(num_cpus::get())
        .build()
        .unwrap()
);

// Process multiple images concurrently
let handles: Vec<_> = paths.iter().map(|path| {
    let pool = Arc::clone(&pool);
    let path = path.clone();
    std::thread::spawn(move || {
        ImageConverter::new()
            .with_thread_pool(pool)
            .process(&path)
    })
}).collect();

let results: Vec<_> = handles.into_iter()
    .map(|h| h.join().unwrap())
    .collect();
```

**Recommendations by concurrency level:**

| Concurrent images | Strategy |
|-------------------|----------|
| 1-3 | Default global pool is fine |
| 4-8 | Shared pool via `with_thread_pool()` with `num_cpus` threads |
| 8+ | Shared pool + limit concurrency with a semaphore or channel |

**Memory budget:** Each full-resolution image (e.g. 4096x3072 16-bit) uses ~150 MB peak. For 10 concurrent images, budget ~1.5 GB. Use `with_preview_mode()` or `with_downscale()` to reduce memory usage.

## Performance

Benchmarks on Apple M4 (6252x4176 16-bit images), wall time of one CLI call, warm cache:

| Mode | Time |
|------|------|
| Mono FITS → JPEG | ~107ms |
| OSC FITS → JPEG | ~67ms |
| XISF → JPEG | ~106ms |
| Mono FITS + annotate | ~480ms |
| `rustafits analyze`, one frame | ~260-390ms |

### SIMD Acceleration

SIMD is used across the processing pipeline with automatic runtime dispatch:

| Operation | SSE2 | AVX2 | NEON |
|-----------|------|------|------|
| Stretch | 4 px/iter | 8 px/iter | 4 px/iter |
| Binning | yes | yes | yes |
| u16 to f32 | yes | yes | yes |
| Gray to RGB | SSSE3 pshufb | AVX2 pshufb | yes |
| Debayer (f32) | yes | — | yes |
| JPEG encode | — | AVX2 DCT | NEON DCT |

## Architecture

```
rustafits/
├── src/
│   ├── lib.rs              # Library entry + public API
│   ├── types.rs            # Core types (PixelData, ProcessedImage, etc.)
│   ├── annotate.rs         # Star annotation overlay (3-tier API over a StarOverlay)
│   ├── converter.rs        # ImageConverter builder
│   ├── pipeline.rs         # Processing pipeline
│   ├── output.rs           # JPEG/PNG file output
│   ├── bin/rustafits.rs    # CLI tool
│   ├── formats/
│   │   ├── mod.rs          # Format dispatch
│   │   ├── fits.rs         # FITS reader
│   │   └── xisf.rs         # XISF reader (zlib/LZ4/Zstd)
│   ├── star_measure/       # Athenaeum's star measurement (analysis_mode, raw: the Analysis numbers)
│   ├── analysis/
│   │   ├── mod.rs            # Legacy analyzer (analyze* deprecated) + detect_fast
│   │   ├── background.rs     # Background estimation (global, mesh-grid, MRS wavelet)
│   │   ├── convolution.rs    # Separable convolution + B3-spline smoothing
│   │   ├── detection.rs      # Star detection (DAOFIND + proximity blend rejection)
│   │   ├── fitting.rs        # LM Gaussian & Moffat PSF fitting (free/fixed beta)
│   │   ├── metrics.rs        # FWHM, eccentricity, HFR measurement
│   │   └── snr.rs            # Per-star and image-wide SNR
│   ├── processing/
│   │   ├── mod.rs           # Processing module
│   │   ├── stretch.rs       # Auto-stretch (SIMD)
│   │   ├── debayer.rs       # Super-pixel debayer + green interpolation (SIMD)
│   │   ├── vng.rs           # Full-resolution 8-gradient VNG demosaic (rayon)
│   │   ├── binning.rs       # 2x2 binning (SIMD)
│   │   ├── downscale.rs     # Integer downscaling
│   │   └── color.rs         # Color conversions (SIMD)
│   └── platesolving/        # Quad pattern matching, RANSAC, WCS, gnomonic
│       │                    # projection, proper motion, transform fitting.
│       │                    # See docs/platesolving.md
│       └── ...
```

**Dependencies**: anyhow, flate2 (rust_backend), lz4_flex, ruzstd, image (PNG only), libjpeg-turbo-rs (pure-Rust JPEG encoder), quick-xml, base64, rayon, nalgebra, tracing

## Troubleshooting

**Slow conversion**: Use `--preview` for mono images or `--downscale 2`

**Black/white output**: the auto-stretch is computed from the frame's median and MAD, so a blank, fully saturated or heavily clipped frame stretches to flat black or white. (`--log` shows the conversion settings, not the stretch parameters.)

**Downscale + Bayer/OSC**: The super-pixel debayer already halves resolution (2x). A `--downscale 2` on a Bayer image produces debayer-only output with no extra downscale. Use `--downscale 4` or higher for additional reduction beyond debayering.

## References

- [FITS Standard](https://fits.gsfc.nasa.gov/)
- [XISF Specification](https://pixinsight.com/xisf/)
- Stetson, P.B. (1987) — DAOFIND star detection algorithm
- SExtractor — Background estimation methodology

## License

Apache-2.0
