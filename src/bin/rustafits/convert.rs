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
use astroimage::{
    annotate_image, encode_jpeg, encode_png, AnnotationConfig, ImageAnalyzer, ImageConverter,
    OutputFormat, PngCompression, ProcessedImage, ThreadPool,
};
use clap::{Args, ValueEnum};

use crate::batch::{self, file_name, has_extension, Input, INPUT_EXTENSIONS};

/// Files converted at once when `--jobs` is not given. Every file already uses
/// all cores for its pixel work; extra jobs mainly overlap the single-threaded
/// JPEG encoder and file I/O. Measured on 4 cores, 9 real 26 MP frames: plain
/// conversion 2.1 s at 1 job, 0.84 s at 4 (peak RSS 260 → 650 MB); with
/// --annotate the analysis saturates the cores by itself, 11.5 s at 1 job,
/// 8.6 s at 2 or 4, but 1.5 GB at 2 and 2.2 GB at 4 — hence a lower cap.
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

    /// Max stars kept for annotation
    #[arg(long, default_value_t = 500, value_name = "N", value_parser = clap::value_parser!(u32).range(1..))]
    max_stars: u32,

    /// Max stars PSF-fitted for annotation (0 = all)
    #[arg(long, default_value_t = 2000, value_name = "N")]
    measure_cap: usize,

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
                    eprintln!("error: {}: {e:#}", job.input.display());
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
    read: Duration,
    analysis: Option<(Duration, usize, f32)>,
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

    let t = Instant::now();
    let (meta, pixels) = ImageConverter::read_raw(&job.input)?;
    let read = t.elapsed();

    let analysis = if args.annotate {
        let t = Instant::now();
        let result = ImageAnalyzer::new()
            .with_max_stars(args.max_stars as usize)
            .with_measure_cap(args.measure_cap)
            .with_thread_pool(pool.clone())
            .analyze_raw(&meta, &pixels)
            .context("analysis failed")?;
        Some((result, t.elapsed()))
    } else {
        None
    };

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
    if let Some((result, _)) = &analysis {
        annotate_image(&mut image, result, &AnnotationConfig::default());
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
        analysis: analysis.map(|(r, d)| (d, r.stars_detected, r.median_fwhm)),
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
