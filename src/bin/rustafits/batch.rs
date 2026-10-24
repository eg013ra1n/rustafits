//! What every subcommand shares: turning paths into input files, the thread
//! pool, the per-file worker loop and the progress bar.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use anyhow::{bail, Context, Result};
use astroimage::{ThreadPool, ThreadPoolBuilder};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};

pub const INPUT_EXTENSIONS: &[&str] = &["fits", "fit", "fts", "xisf"];

/// One input file, with its path relative to the directory it was found in
/// (or just its file name when given directly) for mirroring under `-o DIR`.
pub struct Input {
    pub path: PathBuf,
    pub rel: PathBuf,
}

/// Files as given, directories expanded (`recursive` descends; entries sorted).
/// A missing path is an error; an empty directory a warning; nothing at all an error.
pub fn collect_inputs(paths: &[PathBuf], recursive: bool) -> Result<Vec<Input>> {
    let mut inputs = Vec::new();
    for path in paths {
        if path.is_dir() {
            let before = inputs.len();
            collect_dir(path, path, recursive, &mut inputs)?;
            if inputs.len() == before {
                eprintln!("warning: no FITS/XISF files in {}", path.display());
            }
        } else if path.is_file() {
            inputs.push(Input { path: path.clone(), rel: PathBuf::from(path.file_name().unwrap_or_default()) });
        } else {
            bail!("input not found: {}", path.display());
        }
    }
    if inputs.is_empty() {
        bail!("no FITS/XISF input files");
    }
    Ok(inputs)
}

fn collect_dir(root: &Path, dir: &Path, recursive: bool, out: &mut Vec<Input>) -> Result<()> {
    let mut entries: Vec<PathBuf> = std::fs::read_dir(dir)
        .and_then(|rd| rd.map(|e| e.map(|e| e.path())).collect::<std::io::Result<_>>())
        .with_context(|| format!("cannot read directory {}", dir.display()))?;
    entries.sort();
    for path in entries {
        if path.is_dir() {
            if recursive {
                collect_dir(root, &path, recursive, out)?;
            }
        } else if has_extension(&path, INPUT_EXTENSIONS) {
            let rel = path.strip_prefix(root).unwrap_or(&path).to_path_buf();
            out.push(Input { path, rel });
        }
    }
    Ok(())
}

pub fn has_extension(path: &Path, exts: &[&str]) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| exts.iter().any(|x| e.eq_ignore_ascii_case(x)))
}

pub fn file_name(path: &Path) -> String {
    path.file_name().map_or_else(|| path.display().to_string(), |n| n.to_string_lossy().into_owned())
}

/// The rayon pool all jobs share (`threads` None = all cores).
pub fn thread_pool(threads: Option<u32>) -> Result<Arc<ThreadPool>> {
    Ok(Arc::new(
        ThreadPoolBuilder::new()
            .num_threads(threads.map_or(0, |n| n as usize))
            .build()
            .context("failed to create thread pool")?,
    ))
}

/// Concurrent files: `requested`, else `default_cap` capped by the cores;
/// never more than there are items, never zero.
pub fn workers(requested: Option<u32>, default_cap: usize, items: usize) -> usize {
    let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
    requested.map_or(default_cap.min(cores), |n| n as usize).min(items).max(1)
}

/// A progress bar on stderr for batches (hidden for one item, with `quiet`,
/// or when stderr is not a terminal).
pub fn progress_bar(len: usize, quiet: bool) -> ProgressBar {
    if len <= 1 || quiet {
        return ProgressBar::hidden();
    }
    let bar = ProgressBar::with_draw_target(Some(len as u64), ProgressDrawTarget::stderr());
    bar.set_style(
        ProgressStyle::with_template("{elapsed_precise} [{bar:30}] {pos}/{len} eta {eta}  {msg}")
            .expect("valid template")
            .progress_chars("=> "),
    );
    bar
}

/// Run `work` over `items` on `workers` threads (each item once, in order of
/// pickup). `done` runs as each item finishes, one at a time, with the bar
/// suspended so its output never tears the bar. A panic in `work` becomes an
/// error for that item. Results come back in input order.
pub fn run_parallel<T: Sync, R: Send>(
    items: &[T],
    workers: usize,
    bar: &ProgressBar,
    label: impl Fn(&T) -> String + Sync,
    work: impl Fn(&T) -> Result<R> + Sync,
    done: impl Fn(&T, &Result<R>) + Sync,
) -> Vec<Result<R>> {
    let next = AtomicUsize::new(0);
    let slots: Vec<Mutex<Option<Result<R>>>> = items.iter().map(|_| Mutex::new(None)).collect();
    let report = Mutex::new(());

    std::thread::scope(|s| {
        for _ in 0..workers {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                let Some(item) = items.get(i) else { break };
                bar.set_message(label(item));
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| work(item)))
                    .unwrap_or_else(|_| Err(anyhow::anyhow!("internal error (panic)")));
                {
                    let _serial = report.lock().unwrap();
                    bar.suspend(|| done(item, &result));
                    bar.inc(1);
                }
                *slots[i].lock().unwrap() = Some(result);
            });
        }
    });
    bar.finish_and_clear();

    slots
        .into_iter()
        .map(|slot| slot.into_inner().unwrap().expect("every item is processed"))
        .collect()
}
