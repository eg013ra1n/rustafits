//! Conversion to JPEG/PNG — the default command.
//!
//! One file (`rustafits in.fits out.jpg`, the classic form) or many
//! (`rustafits *.fits -o previews/ -j 4`): inputs are planned up front (output
//! names, collisions, unsupported extensions), then converted by `--jobs`
//! workers sharing one rayon pool. A failing file is reported and the rest
//! carry on; the exit code is 1 if any file failed.

use std::collections::HashMap;
use std::io::{IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use anyhow::{bail, Context, Result};
use astroimage::star_measure::analysis_mode::FrameOptions;
use astroimage::{
    annotate_image_with, encode_jpeg, encode_png, AnnotationConfig, ImageConverter, OutputFormat,
    PngCompression, ProcessedImage, StarOverlay, ThreadPool,
};
use clap::{Args, ValueEnum};

use crate::batch::{self, file_name, has_extension, Input, INPUT_EXTENSIONS};
use crate::measure::{budget_parser, measure_file, report_file_error, warn_ignored, MEASURE_CAP_WHY};

/// Files converted at once when `--jobs` is not given. Every file already uses
/// all cores for its pixel work; extra jobs mainly overlap the single-threaded
/// JPEG encoder and file I/O. Measured on 4 cores, 9 real 26 MP frames: plain
/// conversion 2.1 s at 1 job, 0.84 s at 4 (peak RSS 260 → 650 MB).
///
/// With --annotate, and for `analyze` (it takes `DEFAULT_MAX_JOBS_ANALYSIS`),
/// the measurement saturates the cores by itself, so the cap is lower.
/// Measured with the 1.8.0 measurement (`star_measure::analysis_mode`) on 10
/// cores, release build, four copies of a 26 MP mono frame: `convert -a` 2.5 s
/// and 0.86 GB peak RSS at 1 job, 1.9 s and 1.6 GB at 2, 1.8 s and 2.1 GB at 3
/// or 4; `analyze` 2.0 s and 0.86 GB at 1, 1.7 s and 1.5 GB at 2. The jobs
/// beyond 2 buy under 10 % and cost 0.5 GB, so the default stays 2: its peak
/// is under 2.25 GB, the 1.5 GB the pre-1.8.0 analyzer used at 2 jobs plus half.
const DEFAULT_MAX_JOBS: usize = 4;
pub const DEFAULT_MAX_JOBS_ANALYSIS: usize = 2;

#[derive(Args, Debug, Clone)]
pub struct ConvertArgs {
    /// Input FITS/XISF files or directories (.fits .fit .fts .xisf). With
    /// exactly two paths, the second is the output if it is `-`, ends in
    /// .jpg/.jpeg/.png, or does not exist (the classic `in out` form).
    #[arg(required = true, value_name = "INPUT")]
    paths: Vec<PathBuf>,

    /// Output file (one input; `-` = stdout) or directory (several inputs, an
    /// existing directory, or a path ending in `/`). Default: next to each input
    #[arg(short, long, value_name = "PATH")]
    output: Option<PathBuf>,

    /// Output format for derived names and stdout [default: jpg]. An explicit
    /// output file name decides its own format
    #[arg(short, long, value_enum)]
    format: Option<Format>,

    /// Descend into subdirectories of directory inputs (mirrored under -o DIR)
    #[arg(short, long)]
    recursive: bool,

    /// Skip inputs whose output file already exists
    #[arg(short = 'n', long)]
    skip_existing: bool,

    /// Files converted concurrently [default: min(4, cores, files); 2 with --annotate]
    #[arg(short, long, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    jobs: Option<u32>,

    /// Threads for pixel processing, shared by all jobs [default: all cores]
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    threads: Option<u32>,

    /// JPEG quality 1-100
    #[arg(short, long, default_value_t = 95, value_parser = clap::value_parser!(u8).range(1..=100))]
    quality: u8,

    /// PNG compression (astro frames barely compress at `fast`; `none` is
    /// ~2.5x faster, `balanced` ~2.5x smaller but ~30x slower)
    #[arg(long, value_enum, default_value_t = PngLevel::Fast, value_name = "LEVEL")]
    png_compression: PngLevel,

    /// Integer downscale factor. Bayer images: the 2x2 debayer already counts
    /// as 2, so 2 = debayer only, 4 = debayer + 2x
    #[arg(short, long, default_value_t = 1, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    downscale: u32,

    /// 2x2 binning for mono images (faster, half size)
    #[arg(short, long)]
    preview: bool,

    /// Keep the raw Bayer mosaic (no debayering)
    #[arg(long)]
    no_debayer: bool,

    /// Debayer at native resolution (VNG) instead of 2x2 super-pixel
    #[arg(long, conflicts_with = "no_debayer")]
    vng: bool,

    /// Overlay star detection ellipses on the output
    #[arg(short, long)]
    annotate: bool,

    /// Star budget for --annotate: the brightest good stars drawn, 10-2000
    #[arg(long, default_value_t = 500, value_name = "N", value_parser = budget_parser)]
    max_stars: usize,

    // Flag of the pre-1.8.0 analyzer, kept so old scripts still run (gone in 2.0).
    /// Ignored since 1.8.0: the star budget is --max-stars
    #[arg(long, value_name = "N", hide = true)]
    measure_cap: Option<usize>,

    /// Per-file details and timings, on stderr (`--log` is an alias)
    #[arg(short, long, alias = "log")]
    verbose: bool,

    /// Errors only: no progress bar, no summary
    #[arg(long, conflicts_with = "verbose")]
    quiet: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum Format {
    Jpg,
    Png,
}

impl Format {
    fn extension(self) -> &'static str {
        match self {
            Format::Jpg => "jpg",
            Format::Png => "png",
        }
    }

    fn output_format(self) -> OutputFormat {
        match self {
            Format::Jpg => OutputFormat::Jpeg,
            Format::Png => OutputFormat::Png,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum PngLevel {
    None,
    Fast,
    Balanced,
    High,
}

impl From<PngLevel> for PngCompression {
    fn from(l: PngLevel) -> Self {
        match l {
            PngLevel::None => PngCompression::None,
            PngLevel::Fast => PngCompression::Fast,
            PngLevel::Balanced => PngCompression::Balanced,
            PngLevel::High => PngCompression::High,
        }
    }
}

#[derive(Debug)]
enum Target {
    File(PathBuf),
    Stdout,
}

impl Target {
    fn display(&self) -> String {
        match self {
            Target::File(p) => p.display().to_string(),
            Target::Stdout => "<stdout>".to_string(),
        }
    }
}

#[derive(Debug)]
struct Job {
    input: PathBuf,
    target: Target,
    format: OutputFormat,
}

#[derive(Default)]
struct Tally {
    converted: usize,
    skipped: usize,
    failed: usize,
}

/// `Ok(false)`: the run happened but at least one file failed.
pub fn run(args: ConvertArgs) -> Result<bool> {
    // An old script keeps running; the flag is said once, even with --quiet.
    if args.measure_cap.is_some() {
        warn_ignored("--measure-cap", MEASURE_CAP_WHY);
    }
    let jobs = plan(&args)?;
    let pool = batch::thread_pool(args.threads)?;
    let default_jobs = if args.annotate { DEFAULT_MAX_JOBS_ANALYSIS } else { DEFAULT_MAX_JOBS };
    let workers = batch::workers(args.jobs, default_jobs, jobs.len());
    let bar = batch::progress_bar(jobs.len(), args.quiet);

    let started = Instant::now();
    let tally = Mutex::new(Tally::default());
    batch::run_parallel(
        &jobs,
        workers,
        &bar,
        |job| file_name(&job.input),
        |job| convert_one(job, &args, &pool),
        |job, outcome| {
            let mut t = tally.lock().unwrap();
            match outcome {
                Ok(Outcome::Converted(report)) => {
                    t.converted += 1;
                    if args.verbose {
                        eprintln!("{}", report.line(job));
                    }
                }
                Ok(Outcome::Skipped) => {
                    t.skipped += 1;
                    if args.verbose {
                        eprintln!("skip {} (exists)", job.target.display());
                    }
                }
                Err(e) => {
                    t.failed += 1;
                    report_file_error(&job.input, e);
                }
            }
        },
    );

    let t = tally.into_inner().unwrap();
    if !args.quiet && (jobs.len() > 1 || args.verbose) {
        let mut summary = format!("{} converted", t.converted);
        if t.skipped > 0 {
            summary += &format!(", {} skipped", t.skipped);
        }
        if t.failed > 0 {
            summary += &format!(", {} failed", t.failed);
        }
        eprintln!("{summary} in {:.1} s", started.elapsed().as_secs_f64());
    }
    Ok(t.failed == 0)
}

// ---------------------------------------------------------------- planning

/// Resolve inputs and output paths, and reject everything that can be caught
/// before any pixel work: missing inputs, unsupported output extensions,
/// two inputs mapping to one output, `-` with several inputs.
fn plan(args: &ConvertArgs) -> Result<Vec<Job>> {
    let (input_paths, positional_output) = split_positionals(&args.paths);
    let output = match (positional_output, &args.output) {
        (Some(_), Some(_)) => bail!("output given twice: as the second path and with -o"),
        (p, o) => p.or(o.clone()),
    };
    let mut inputs = batch::collect_inputs(&input_paths, args.recursive)?;

    let format = args.format.unwrap_or(Format::Jpg);
    let derived = |base: &Path| base.with_extension(format.extension());

    let jobs: Vec<Job> = match output {
        Some(o) if o.as_os_str() == "-" => {
            if inputs.len() != 1 {
                bail!("`-` (stdout) takes exactly one input, got {}", inputs.len());
            }
            if std::io::stdout().is_terminal() {
                bail!("refusing to write image data to a terminal; redirect stdout or use -o FILE");
            }
            vec![Job { input: inputs.remove(0).path, target: Target::Stdout, format: format.output_format() }]
        }
        Some(o) if inputs.len() > 1 || o.is_dir() || ends_with_separator(&o) => inputs
            .into_iter()
            .map(|Input { path, rel }| Job {
                target: Target::File(o.join(derived(&rel))),
                format: format.output_format(),
                input: path,
            })
            .collect(),
        Some(o) => {
            let out_format = OutputFormat::from_path(&o)?;
            if let Some(f) = args.format {
                if f.output_format() != out_format {
                    bail!("--format {} conflicts with the output name {}", f.extension(), o.display());
                }
            }
            vec![Job { input: inputs.remove(0).path, target: Target::File(o), format: out_format }]
        }
        None => inputs
            .into_iter()
            .map(|Input { path, .. }| Job {
                target: Target::File(derived(&path)),
                format: format.output_format(),
                input: path,
            })
            .collect(),
    };

    let mut seen: HashMap<&Path, &Path> = HashMap::new();
    for job in &jobs {
        if let Target::File(out) = &job.target {
            if let Some(first) = seen.insert(out, &job.input) {
                bail!(
                    "{} and {} would both be written to {}",
                    first.display(),
                    job.input.display(),
                    out.display()
                );
            }
        }
    }
    Ok(jobs)
}

/// The classic `rustafits in.fits out.jpg`: with exactly two paths, the second
/// is the output when it is `-`, has an output extension, or does not exist
/// (and is not named like an input). Otherwise every path is an input.
fn split_positionals(paths: &[PathBuf]) -> (Vec<PathBuf>, Option<PathBuf>) {
    if let [input, second] = paths {
        let is_output = second.as_os_str() == "-"
            || has_extension(second, &["jpg", "jpeg", "png"])
            || (!second.exists() && !has_extension(second, INPUT_EXTENSIONS));
        if is_output {
            return (vec![input.clone()], Some(second.clone()));
        }
    }
    (paths.to_vec(), None)
}

fn ends_with_separator(path: &Path) -> bool {
    path.as_os_str()
        .to_str()
        .is_some_and(|s| s.ends_with('/') || s.ends_with(std::path::MAIN_SEPARATOR))
}

// -------------------------------------------------------------- conversion

enum Outcome {
    Converted(Report),
    Skipped,
}

struct Report {
    width: usize,
    height: usize,
    /// The render's read of the file. With `--annotate` the measurement reads the file too, before it, and that read
    /// is inside `analysis`.
    read: Duration,
    /// `--annotate`: the measurement's time (its own read of the file included), the stars detected, the median FWHM.
    analysis: Option<(Duration, usize, f64)>,
    process: Duration,
    encode: Duration,
    write: Duration,
}

impl Report {
    fn line(&self, job: &Job) -> String {
        let ms = |d: Duration| d.as_secs_f64() * 1e3;
        let mut s = format!(
            "{} -> {}  {}x{}  read {:.0} ms",
            job.input.display(),
            job.target.display(),
            self.width,
            self.height,
            ms(self.read)
        );
        if let Some((d, stars, fwhm)) = self.analysis {
            s += &format!(" | analysis {:.0} ms ({stars} stars, FWHM {fwhm:.2})", ms(d));
        }
        s += &format!(
            " | process {:.0} ms | encode {:.0} ms | write {:.0} ms",
            ms(self.process),
            ms(self.encode),
            ms(self.write)
        );
        s
    }
}

fn convert_one(job: &Job, args: &ConvertArgs, pool: &Arc<ThreadPool>) -> Result<Outcome> {
    if args.skip_existing {
        if let Target::File(out) = &job.target {
            if out.exists() {
                return Ok(Outcome::Skipped);
            }
        }
    }

    // The measurement comes first and returns without the planes, so the render's read below starts after they are
    // freed and the two copies of a 26 MP frame never overlap in memory. (The file is read twice: the measurement's
    // reader states the ADU facts, the render's does not; the measurement dominates the time.)
    let measured = if args.annotate {
        let t = Instant::now();
        let mut opts = FrameOptions::default();
        opts.max_stars = args.max_stars;
        let m = measure_file(&job.input, &opts, pool)?;
        Some((m, t.elapsed()))
    } else {
        None
    };

    let t = Instant::now();
    let (meta, pixels) = ImageConverter::read_raw(&job.input)?;
    let read = t.elapsed();

    // The star positions are in file row order and the render flips by the file's own flag: the two readers must
    // agree, or the rings would be mirrored.
    if let Some((m, _)) = &measured {
        if m.flip_vertical != meta.flip_vertical {
            bail!(
                "cannot place the stars: the measurement reads the rows {} and the render {}",
                if m.flip_vertical { "bottom-up" } else { "top-down" },
                if meta.flip_vertical { "bottom-up" } else { "top-down" },
            );
        }
    }

    let mut converter = ImageConverter::new()
        .with_downscale(args.downscale as usize)
        .with_thread_pool(pool.clone());
    if args.no_debayer {
        converter = converter.without_debayer();
    }
    if args.vng {
        converter = converter.with_vng_debayer();
    }
    if args.preview {
        converter = converter.with_preview_mode();
    }

    let t = Instant::now();
    let mut image = converter.process_data(meta, pixels)?;
    if let Some((m, _)) = &measured {
        let overlay = StarOverlay::from_measurement(m.frame_w, m.frame_h, &m.pm, &m.summary);
        annotate_image_with(&mut image, &overlay, &AnnotationConfig::default());
    }
    let process = t.elapsed();

    let t = Instant::now();
    let bytes = encode(&image, job.format, args)?;
    let encode = t.elapsed();

    let t = Instant::now();
    match &job.target {
        Target::File(out) => write_atomically(out, &bytes)?,
        Target::Stdout => {
            let mut stdout = std::io::stdout().lock();
            stdout.write_all(&bytes).and_then(|_| stdout.flush()).context("failed to write to stdout")?;
        }
    }
    let write = t.elapsed();

    Ok(Outcome::Converted(Report {
        width: image.width,
        height: image.height,
        read,
        analysis: measured.map(|(m, d)| (d, m.summary.stars_detected, m.summary.median_fwhm)),
        process,
        encode,
        write,
    }))
}

fn encode(image: &ProcessedImage, format: OutputFormat, args: &ConvertArgs) -> Result<Vec<u8>> {
    let channels = image.channels as usize;
    match format {
        OutputFormat::Jpeg => encode_jpeg(&image.data, image.width, image.height, channels, args.quality),
        OutputFormat::Png => {
            encode_png(&image.data, image.width, image.height, channels, args.png_compression.into())
        }
    }
}

/// Write via a temporary file in the same directory and rename it into place,
/// so an interrupted run never leaves a truncated image that `--skip-existing`
/// would then take for done.
fn write_atomically(out: &Path, bytes: &[u8]) -> Result<()> {
    let dir = out.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
    std::fs::create_dir_all(dir).with_context(|| format!("cannot create directory {}", dir.display()))?;
    let tmp = dir.join(format!(".{}.{}.tmp", file_name(out), std::process::id()));
    let result = std::fs::write(&tmp, bytes)
        .and_then(|_| std::fs::rename(&tmp, out))
        .with_context(|| format!("cannot write {}", out.display()));
    if result.is_err() {
        std::fs::remove_file(&tmp).ok();
    }
    result
}
