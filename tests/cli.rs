//! End-to-end tests of the `rustafits` binary on small synthetic FITS files
//! written by the tests themselves (no fixtures needed).
#![cfg(feature = "cli")]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use astroimage::{ImageAnalyzer, ImageConverter};

const BIN: &str = env!("CARGO_BIN_EXE_rustafits");

/// A fresh, empty scratch directory per test.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rustafits-cli-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Minimal 16-bit FITS (BZERO 32768) with a gradient, a few "stars" and an
/// optional Bayer pattern — enough structure for a non-flat stretch.
fn write_fits(path: &Path, w: usize, h: usize, bayer: Option<&str>) {
    let mut extra = Vec::new();
    if let Some(p) = bayer {
        extra.push(format!("BAYERPAT= '{p}'"));
    }
    write_fits_with(path, w, h, &extra, |x, y| {
        let noise = ((x * 7919 + y * 104_729) % 97) as f64;
        let star = if (x % 17 == 5) && (y % 13 == 4) { 20_000.0 } else { 0.0 };
        1500.0 + (x + y) as f64 * 3.0 + noise + star
    });
}

/// A star field the analyzer measures: round Gaussian stars of width `sigma`
/// on a 24 px grid over a noisy background, plus header keywords.
fn write_star_field(path: &Path, sigma: f64, extra: &[String]) {
    write_fits_with(path, 320, 240, extra, |x, y| {
        let noise = (((x * 7919 + y * 104_729) % 61) as f64 - 30.0) * 1.5;
        let (cx, cy) = (((x / 24) * 24 + 12) as f64, ((y / 24) * 24 + 12) as f64);
        let amp = 4000.0 + ((x / 24 * 7 + y / 24 * 13) % 9) as f64 * 1500.0;
        let r2 = (x as f64 - cx).powi(2) + (y as f64 - cy).powi(2);
        1000.0 + noise + amp * (-r2 / (2.0 * sigma * sigma)).exp()
    });
}

fn write_fits_with(path: &Path, w: usize, h: usize, extra: &[String], pixel: impl Fn(usize, usize) -> f64) {
    let mut cards = vec![
        "SIMPLE  =                    T".to_string(),
        "BITPIX  =                   16".to_string(),
        "NAXIS   =                    2".to_string(),
        format!("NAXIS1  = {w:>20}"),
        format!("NAXIS2  = {h:>20}"),
        "BZERO   =                32768".to_string(),
        "BSCALE  =                    1".to_string(),
    ];
    cards.extend(extra.iter().cloned());
    cards.push("END".to_string());
    let header: String = cards.iter().map(|c| format!("{c:<80}")).collect();
    let mut bytes = format!("{header:<width$}", width = header.len().next_multiple_of(2880)).into_bytes();
    for y in 0..h {
        for x in 0..w {
            let v = pixel(x, y).round().clamp(0.0, 65_535.0) as i64;
            bytes.extend_from_slice(&((v - 32_768) as i16).to_be_bytes());
        }
    }
    bytes.resize(bytes.len().next_multiple_of(2880), 0);
    std::fs::write(path, bytes).unwrap();
}

fn run(dir: &Path, args: &[&str]) -> Output {
    Command::new(BIN).args(args).current_dir(dir).output().unwrap()
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).into_owned()
}

fn jpeg_size(path: &Path) -> (u16, u16) {
    let mut d = jpeg_decoder::Decoder::new(std::io::BufReader::new(std::fs::File::open(path).unwrap()));
    d.decode().unwrap();
    let info = d.info().unwrap();
    (info.width, info.height)
}

#[test]
fn classic_form_writes_what_the_library_writes() {
    let dir = scratch("classic");
    write_fits(&dir.join("m.fits"), 64, 48, None);
    let o = run(&dir, &["m.fits", "out.jpg"]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(o.stdout.is_empty() && o.stderr.is_empty(), "a single conversion stays silent");

    ImageConverter::new().convert(dir.join("m.fits"), dir.join("lib.jpg")).unwrap();
    assert_eq!(std::fs::read(dir.join("out.jpg")).unwrap(), std::fs::read(dir.join("lib.jpg")).unwrap());

    // Extensionless output (classic form) is still a JPEG.
    assert!(run(&dir, &["m.fits", "noext"]).status.success());
    assert_eq!(&std::fs::read(dir.join("noext")).unwrap()[..2], &[0xFF, 0xD8]);
}

#[test]
fn default_output_lands_next_to_the_input() {
    let dir = scratch("default");
    write_fits(&dir.join("m.fits"), 64, 48, None);
    assert!(run(&dir, &["m.fits"]).status.success());
    assert!(dir.join("m.jpg").is_file());

    assert!(run(&dir, &["m.fits", "-f", "png"]).status.success());
    assert_eq!(&std::fs::read(dir.join("m.png")).unwrap()[..8], b"\x89PNG\r\n\x1a\n");
}

#[test]
fn two_existing_inputs_are_a_batch_not_input_and_output() {
    let dir = scratch("two-inputs");
    write_fits(&dir.join("a.fits"), 64, 48, None);
    write_fits(&dir.join("b.fits"), 64, 48, None);
    let o = run(&dir, &["a.fits", "b.fits", "--quiet"]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(dir.join("a.jpg").is_file() && dir.join("b.jpg").is_file());
}

#[test]
fn directories_batch_and_recursion_mirrors_the_tree() {
    let dir = scratch("dirs");
    std::fs::create_dir_all(dir.join("in/sub")).unwrap();
    write_fits(&dir.join("in/a.fits"), 64, 48, None);
    write_fits(&dir.join("in/sub/b.FIT"), 64, 48, None);
    std::fs::write(dir.join("in/notes.txt"), "not an image").unwrap();

    let o = run(&dir, &["in", "-o", "flat/", "--quiet"]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(dir.join("flat/a.jpg").is_file());
    assert!(!dir.join("flat/sub").exists(), "no recursion without -r");

    let o = run(&dir, &["in", "-r", "-o", "tree"]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(dir.join("tree/a.jpg").is_file() && dir.join("tree/sub/b.jpg").is_file());
    assert!(stderr(&o).contains("2 converted"), "{}", stderr(&o));
    let leftovers: Vec<_> = std::fs::read_dir(dir.join("tree")).unwrap().flatten()
        .filter(|e| e.file_name().to_string_lossy().ends_with(".tmp"))
        .collect();
    assert!(leftovers.is_empty(), "temp files left behind");
}

#[test]
fn skip_existing_leaves_outputs_untouched() {
    let dir = scratch("skip");
    write_fits(&dir.join("a.fits"), 64, 48, None);
    write_fits(&dir.join("b.fits"), 64, 48, None);
    std::fs::create_dir_all(dir.join("out")).unwrap();
    std::fs::write(dir.join("out/a.jpg"), "sentinel").unwrap();

    let o = run(&dir, &["a.fits", "b.fits", "-o", "out", "-n"]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert_eq!(std::fs::read_to_string(dir.join("out/a.jpg")).unwrap(), "sentinel");
    assert!(dir.join("out/b.jpg").is_file());
    assert!(stderr(&o).contains("1 converted, 1 skipped"), "{}", stderr(&o));
}

#[test]
fn stdout_carries_the_same_bytes_as_a_file() {
    let dir = scratch("stdout");
    write_fits(&dir.join("m.fits"), 64, 48, None);
    let o = run(&dir, &["m.fits", "-o", "-"]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(run(&dir, &["m.fits", "-o", "f.jpg"]).status.success());
    assert_eq!(o.stdout, std::fs::read(dir.join("f.jpg")).unwrap());

    let o = run(&dir, &["m.fits", "-o", "-", "-f", "png"]);
    assert_eq!(&o.stdout[..8], b"\x89PNG\r\n\x1a\n");
}

#[test]
fn bad_output_names_fail_before_any_work() {
    let dir = scratch("bad-names");
    write_fits(&dir.join("m.fits"), 64, 48, None);
    std::fs::create_dir_all(dir.join("d1")).unwrap();
    std::fs::create_dir_all(dir.join("d2")).unwrap();
    write_fits(&dir.join("d1/x.fits"), 64, 48, None);
    write_fits(&dir.join("d2/x.fits"), 64, 48, None);

    let o = run(&dir, &["m.fits", "out.tif"]);
    assert!(!o.status.success() && stderr(&o).contains("unsupported output format"));
    assert!(!dir.join("out.tif").exists());

    let o = run(&dir, &["m.fits", "-o", "x.png", "-f", "jpg"]);
    assert!(!o.status.success() && stderr(&o).contains("conflicts"));

    let o = run(&dir, &["d1/x.fits", "d2/x.fits", "-o", "out/"]);
    assert!(!o.status.success() && stderr(&o).contains("would both be written"));
    assert!(!dir.join("out").exists(), "nothing may be written");

    let o = run(&dir, &["m.fits", "d1/x.fits", "-o", "-"]);
    assert!(!o.status.success() && stderr(&o).contains("exactly one input"));

    let o = run(&dir, &["missing.fits"]);
    assert!(!o.status.success() && stderr(&o).contains("input not found"));
}

#[test]
fn a_failing_file_does_not_stop_the_batch() {
    let dir = scratch("partial");
    write_fits(&dir.join("good1.fits"), 64, 48, None);
    write_fits(&dir.join("good2.fits"), 64, 48, None);
    std::fs::write(dir.join("broken.fits"), b"SIMPLE  =                    T").unwrap();

    let o = run(&dir, &["good1.fits", "broken.fits", "good2.fits", "-o", "out", "-j", "1"]);
    assert_eq!(o.status.code(), Some(1), "any failure → exit 1");
    assert!(dir.join("out/good1.jpg").is_file() && dir.join("out/good2.jpg").is_file());
    assert!(!dir.join("out/broken.jpg").exists());
    let err = stderr(&o);
    assert!(err.contains("broken.fits") && err.contains("2 converted, 1 failed"), "{err}");
}

#[test]
fn debayer_flags_set_the_output_size() {
    let dir = scratch("debayer");
    write_fits(&dir.join("osc.fits"), 64, 48, Some("RGGB"));
    for (args, size) in [
        (vec!["osc.fits", "sp.jpg"], (32, 24)),
        (vec!["osc.fits", "vng.jpg", "--vng"], (64, 48)),
        (vec!["osc.fits", "raw.jpg", "--no-debayer"], (64, 48)),
        (vec!["osc.fits", "ds.jpg", "-d", "4"], (16, 12)),
    ] {
        let o = run(&dir, &args);
        assert!(o.status.success(), "{args:?}: {}", stderr(&o));
        assert_eq!(jpeg_size(&dir.join(args[1])), size, "{args:?}");
    }
    let o = run(&dir, &["osc.fits", "--vng", "--no-debayer"]);
    assert_eq!(o.status.code(), Some(2), "clap usage error");
}

#[test]
fn verbose_reports_each_file_on_stderr() {
    let dir = scratch("verbose");
    write_fits(&dir.join("m.fits"), 64, 48, None);
    for flag in ["-v", "--log"] {
        let o = run(&dir, &["m.fits", "out.jpg", flag]);
        assert!(o.status.success());
        assert!(o.stdout.is_empty(), "logs never go to stdout");
        let err = stderr(&o);
        assert!(err.contains("m.fits -> out.jpg") && err.contains("64x48"), "{flag}: {err}");
    }
}

fn stdout(o: &Output) -> String {
    String::from_utf8_lossy(&o.stdout).into_owned()
}

#[test]
fn info_summarises_headers_without_pixels() {
    let dir = scratch("info");
    let extra = |exp: &str, filter: &str| {
        vec![
            format!("EXPTIME = {exp:>20} / [s]"),
            format!("FILTER  = '{filter}'"),
            "FOCALLEN=                448.0".to_string(),
            "COMMENT  synthetic".to_string(),
        ]
    };
    write_star_field(&dir.join("a.fits"), 1.5, &extra("30.0", "Ha"));
    write_star_field(&dir.join("b.fits"), 1.5, &extra("120.0", "OIII"));

    let o = run(&dir, &["info", "a.fits", "b.fits", "-k", "FOCALLEN"]);
    assert!(o.status.success(), "{}", stderr(&o));
    let out = stdout(&o);
    let lines: Vec<&str> = out.lines().collect();
    assert_eq!(lines.len(), 3, "{out}");
    assert!(lines[0].starts_with("FILE") && lines[0].contains("EXP") && lines[0].ends_with("FOCALLEN"));
    assert!(lines[1].contains("a.fits") && lines[1].contains("320x240") && lines[1].contains("30.0") && lines[1].contains("Ha"));
    assert!(lines[2].contains("120.0") && lines[2].contains("OIII") && lines[2].contains("448.0"));

    let o = run(&dir, &["info", "a.fits", "--csv"]);
    assert_eq!(stdout(&o).lines().nth(1).unwrap(), "a.fits,320x240,u16,-,30.0,Ha,-,-,-,-");

    let o = run(&dir, &["info", "a.fits", "--json"]);
    let json = stdout(&o);
    assert!(json.contains("\"width\": 320") && json.contains("\"bayer\": null"), "{json}");
    assert!(json.contains("{\"key\": \"FILTER\", \"value\": \"Ha\", \"comment\": \"\"}"), "{json}");

    let o = run(&dir, &["info", "a.fits", "--all"]);
    let all = stdout(&o);
    assert!(all.contains("EXPTIME  = 30.0  / [s]") && all.contains("COMMENT  = synthetic"), "{all}");
}

/// Parse `analyze --csv` into (header, rows).
fn analyze_csv(dir: &Path, args: &[&str]) -> (Vec<String>, Vec<Vec<String>>) {
    let mut full = vec!["analyze", "--csv", "--quiet"];
    full.extend_from_slice(args);
    let o = run(dir, &full);
    assert!(o.status.success(), "{}", stderr(&o));
    let out = stdout(&o);
    let mut lines = out.lines().map(|l| l.split(',').map(str::to_string).collect::<Vec<_>>());
    let header = lines.next().unwrap();
    (header, lines.collect())
}

#[test]
fn analyze_reports_exactly_what_the_library_measures() {
    let dir = scratch("analyze");
    write_star_field(&dir.join("sharp.fits"), 1.2, &[]);
    write_star_field(&dir.join("soft.fits"), 2.2, &[]);

    let (header, rows) = analyze_csv(&dir, &["sharp.fits", "soft.fits"]);
    let col = |name: &str| header.iter().position(|h| h == name).unwrap();
    for (row, file) in rows.iter().zip(["sharp.fits", "soft.fits"]) {
        let lib = ImageAnalyzer::new().analyze(dir.join(file)).unwrap();
        assert!(lib.stars_detected > 20, "{file}: the field must be measurable");
        assert_eq!(row[col("file")], file);
        assert_eq!(row[col("stars_detected")], lib.stars_detected.to_string());
        assert_eq!(row[col("median_fwhm")], lib.median_fwhm.to_string());
        assert_eq!(row[col("median_hfr")], lib.median_hfr.to_string());
        assert_eq!(row[col("median_eccentricity")], lib.median_eccentricity.to_string());
        assert_eq!(row[col("snr_weight")], lib.snr_weight.to_string());
    }

    // Best-first by FWHM puts the sharp field first whatever the input order.
    let (_, rows) = analyze_csv(&dir, &["soft.fits", "sharp.fits", "--sort", "fwhm"]);
    assert_eq!((rows[0][0].as_str(), rows[1][0].as_str()), ("sharp.fits", "soft.fits"));
    let (_, rows) = analyze_csv(&dir, &["sharp.fits", "soft.fits", "--sort", "fwhm", "--reverse"]);
    assert_eq!(rows[0][0], "soft.fits");

    // Table and JSON carry the same frames.
    let o = run(&dir, &["analyze", "sharp.fits", "soft.fits", "--quiet"]);
    let table = stdout(&o);
    assert!(table.starts_with("FILE") && table.contains("sharp.fits") && table.contains("soft.fits"), "{table}");
    let o = run(&dir, &["analyze", "sharp.fits", "--json"]);
    assert!(stdout(&o).contains("\"file\": \"sharp.fits\""));
}

#[test]
fn analyze_reports_a_broken_file_and_keeps_the_rest() {
    let dir = scratch("analyze-broken");
    write_star_field(&dir.join("good.fits"), 1.5, &[]);
    std::fs::write(dir.join("broken.fits"), b"SIMPLE  =                    T").unwrap();
    let o = run(&dir, &["analyze", "good.fits", "broken.fits", "--csv", "--quiet"]);
    assert_eq!(o.status.code(), Some(1));
    assert!(stderr(&o).contains("broken.fits"));
    assert_eq!(stdout(&o).lines().count(), 2, "header + the good frame");
}

#[test]
fn convert_subcommand_and_completions() {
    let dir = scratch("subcommands");
    write_fits(&dir.join("m.fits"), 64, 48, None);
    assert!(run(&dir, &["convert", "m.fits", "a.jpg"]).status.success());
    assert!(run(&dir, &["m.fits", "b.jpg"]).status.success());
    assert_eq!(std::fs::read(dir.join("a.jpg")).unwrap(), std::fs::read(dir.join("b.jpg")).unwrap());

    for shell in ["bash", "zsh", "fish"] {
        let o = run(&dir, &["completions", shell]);
        assert!(o.status.success(), "{shell}");
        let script = stdout(&o);
        assert!(script.contains("rustafits") && script.contains("analyze") && script.contains("png-compression"), "{shell}");
    }
}
