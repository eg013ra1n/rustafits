//! End-to-end tests of the `rustafits` binary on small synthetic FITS files
//! written by the tests themselves (no fixtures needed).
#![cfg(feature = "cli")]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use astroimage::star_measure::analysis_mode::{measure_raw_frame, FrameOptions, FrameSummary, Verdict};
use astroimage::star_measure::raw::{read_raw_frame, RawFrame, RawPlanes};
use astroimage::star_measure::{PlaneMeasurement, PsfModel, MEASURE_VERSION};
use astroimage::ImageConverter;

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

/// The CSV header: the 19 columns the CLI always had, in their order, then the six `analysis_mode` adds.
const CSV_COLUMNS: [&str; 25] = [
    "file", "width", "height", "stars_detected", "stars_measured", "median_fwhm", "median_hfr",
    "median_eccentricity", "median_snr", "median_beta", "snr_weight", "psf_signal", "frame_snr",
    "background", "noise", "trail_r_squared", "possibly_trailed", "median_fwhm_arcsec", "median_hfr_arcsec",
    "channels", "stars_refused", "frame_verdict", "measure_version", "detection_threshold", "psf_model",
];

/// What the library says about `path`: `read_raw_frame` + `measure_raw_frame` + `FrameSummary::from_measurement`.
/// No result reads a thread count (G-6), so `None` for the pool equals the CLI's.
fn lib_measure(path: &Path, opts: &FrameOptions) -> (RawFrame, PlaneMeasurement, FrameSummary) {
    let frame = read_raw_frame(path).unwrap();
    let pm = measure_raw_frame(&frame, opts, None);
    let summary = FrameSummary::from_measurement(&pm, opts);
    (frame, pm, summary)
}

/// `FrameOptions` from the three CLI knobs (it is `#[non_exhaustive]`: no struct literal outside the crate).
fn frame_options(detection_threshold: f64, psf_model: PsfModel, max_stars: usize) -> FrameOptions {
    let mut opts = FrameOptions::default();
    opts.detection_threshold = detection_threshold;
    opts.psf_model = psf_model;
    opts.max_stars = max_stars;
    opts
}

/// The CSV row the library's numbers make (no optics: the two arcsec cells are empty), cell by cell as
/// `f64::to_string()`. `file` is the first cell as the CLI prints it.
fn lib_row(file: &str, path: &Path, opts: &FrameOptions, psf_model: &str) -> Vec<String> {
    let (frame, pm, s) = lib_measure(path, opts);
    let verdict = match s.verdict {
        Some(Verdict::Trailed) => "trailed",
        Some(Verdict::Defocused) => "defocused",
        Some(Verdict::NoStars) => "no_stars",
        None => "",
    };
    vec![
        file.to_string(),
        frame.w.to_string(),
        frame.h.to_string(),
        s.stars_detected.to_string(),
        pm.stars.len().to_string(),
        s.median_fwhm.to_string(),
        s.median_hfr.to_string(),
        s.median_eccentricity.to_string(),
        s.median_snr.to_string(),
        s.median_beta.map_or(String::new(), |b| b.to_string()),
        s.snr_weight.to_string(),
        s.psf_signal.to_string(),
        s.frame_snr.to_string(),
        s.background.to_string(),
        s.noise.to_string(),
        s.trail_r_squared.to_string(),
        pm.trailed.to_string(),
        String::new(),
        String::new(),
        frame.channels.to_string(),
        s.stars_refused.to_string(),
        verdict.to_string(),
        MEASURE_VERSION.to_string(),
        opts.detection_threshold.to_string(),
        psf_model.to_string(),
    ]
}

/// Compare a CLI row with the library's, cell by cell, naming the column that differs.
fn assert_row_eq(got: &[String], want: &[String], what: &str) {
    assert_eq!(got.len(), want.len(), "{what}: {got:?}");
    for ((g, w), name) in got.iter().zip(want).zip(CSV_COLUMNS) {
        assert_eq!(g, w, "{what}: column {name}");
    }
}

#[test]
fn analyze_reports_exactly_what_the_library_measures() {
    let dir = scratch("analyze");
    write_star_field(&dir.join("sharp.fits"), 1.2, &[]);
    write_star_field(&dir.join("soft.fits"), 2.2, &[]);

    let (header, rows) = analyze_csv(&dir, &["sharp.fits", "soft.fits"]);
    assert_eq!(header, CSV_COLUMNS);
    let opts = FrameOptions::default();
    for (row, file) in rows.iter().zip(["sharp.fits", "soft.fits"]) {
        let (_, pm, s) = lib_measure(&dir.join(file), &opts);
        assert!(pm.stars.len() >= 20 && s.verdict.is_none(), "{file}: the field must be measurable ({s:?})");
        assert_row_eq(row, &lib_row(file, &dir.join(file), &opts, "moffat4"), file);
        // Named, so a drift in the defaults shows here and not only in the row comparison.
        let col = |name: &str| header.iter().position(|h| h == name).unwrap();
        assert_eq!(row[col("median_beta")], "");
        assert_eq!(row[col("frame_verdict")], "");
        assert_eq!(row[col("measure_version")], "8");
        assert_eq!(row[col("psf_model")], "moffat4");
    }

    // Best-first by FWHM puts the sharp field first whatever the input order.
    let (_, rows) = analyze_csv(&dir, &["soft.fits", "sharp.fits", "--sort", "fwhm"]);
    assert_eq!((rows[0][0].as_str(), rows[1][0].as_str()), ("sharp.fits", "soft.fits"));
    let (_, rows) = analyze_csv(&dir, &["sharp.fits", "soft.fits", "--sort", "fwhm", "--reverse"]);
    assert_eq!(rows[0][0], "soft.fits");

    // Table and JSON carry the same frames; the table ends in VERDICT.
    let o = run(&dir, &["analyze", "sharp.fits", "soft.fits", "--quiet"]);
    let table = stdout(&o);
    assert!(table.starts_with("FILE") && table.contains("sharp.fits") && table.contains("soft.fits"), "{table}");
    assert!(table.lines().next().unwrap().ends_with("VERDICT"), "{table}");
    assert!(!table.contains("TRAIL"), "VERDICT replaces TRAIL: {table}");
    let o = run(&dir, &["analyze", "sharp.fits", "--json"]);
    let json = stdout(&o);
    assert!(json.contains("\"file\": \"sharp.fits\""), "{json}");
    assert!(json.contains("\"frame_verdict\": null") && json.contains("\"psf_model\": \"moffat4\""), "{json}");
    assert!(json.contains("\"measure_version\": 8") && json.contains("\"detection_threshold\": 6,"), "{json}");
}

#[test]
fn analyze_options_reach_the_measurement() {
    let dir = scratch("analyze-options");
    write_star_field(&dir.join("sharp.fits"), 1.2, &[]);
    let path = dir.join("sharp.fits");

    // --psf-model auto: the frame's own beta, whatever the library says (Some only for a fitted source).
    let auto = frame_options(6.0, PsfModel::Auto, 500);
    let (_, rows) = analyze_csv(&dir, &["sharp.fits", "--psf-model", "auto"]);
    assert_row_eq(&rows[0], &lib_row("sharp.fits", &path, &auto, "auto"), "--psf-model auto");
    let (_, _, s) = lib_measure(&path, &auto);
    eprintln!("median_beta under auto on sharp.fits: {:?}", s.median_beta);

    // --detection-threshold: the library's detection at 12.0.
    let high = frame_options(12.0, PsfModel::Moffat4, 500);
    let (_, rows) = analyze_csv(&dir, &["sharp.fits", "--detection-threshold", "12"]);
    assert_row_eq(&rows[0], &lib_row("sharp.fits", &path, &high, "moffat4"), "--detection-threshold 12");
    assert_eq!(rows[0][CSV_COLUMNS.iter().position(|c| *c == "detection_threshold").unwrap()], "12");

    // --max-stars: the budget. The default keeps more than 20 on this field, so 20 bites.
    let (_, rows) = analyze_csv(&dir, &["sharp.fits"]);
    let measured = CSV_COLUMNS.iter().position(|c| *c == "stars_measured").unwrap();
    assert!(rows[0][measured].parse::<usize>().unwrap() > 20, "{:?}", rows[0]);
    let small = frame_options(6.0, PsfModel::Moffat4, 20);
    let (_, rows) = analyze_csv(&dir, &["sharp.fits", "--max-stars", "20"]);
    assert_row_eq(&rows[0], &lib_row("sharp.fits", &path, &small, "moffat4"), "--max-stars 20");
    assert!(rows[0][measured].parse::<usize>().unwrap() <= 20, "{:?}", rows[0]);

    // Optics add the arcsec cells: px x pixel size / focal length x 206.265.
    let (_, rows) = analyze_csv(&dir, &["sharp.fits", "--focal-length", "448", "--pixel-size", "3.76"]);
    let (_, _, s) = lib_measure(&path, &FrameOptions::default());
    let scale = 3.76 / 448.0 * 206.265;
    let cell = |name: &str| rows[0][CSV_COLUMNS.iter().position(|c| *c == name).unwrap()].parse::<f64>().unwrap();
    assert!((cell("median_fwhm_arcsec") - s.median_fwhm * scale).abs() < 1e-9);
    assert!((cell("median_hfr_arcsec") - s.median_hfr * scale).abs() < 1e-9);
}

#[test]
fn analyze_measures_a_mosaic_on_its_green_sites() {
    let dir = scratch("analyze-mosaic");
    write_star_field(&dir.join("mosaic.fits"), 1.2, &["BAYERPAT= 'RGGB'".to_string()]);
    let path = dir.join("mosaic.fits");
    let opts = FrameOptions::default();
    let (frame, pm, s) = lib_measure(&path, &opts);
    assert!(matches!(frame.planes, RawPlanes::Cfa { .. }), "the file reads as a mosaic");
    assert!(pm.stars.len() >= 20 && s.verdict.is_none(), "the mosaic must be measurable ({s:?})");

    let (_, rows) = analyze_csv(&dir, &["mosaic.fits"]);
    assert_row_eq(&rows[0], &lib_row("mosaic.fits", &path, &opts, "moffat4"), "mosaic");
    // The native size, one channel: not the half-size super-pixel render.
    let cell = |name: &str| rows[0][CSV_COLUMNS.iter().position(|c| *c == name).unwrap()].as_str();
    assert_eq!((cell("width"), cell("height"), cell("channels")), ("320", "240", "1"));
}

/// A frame of noise and no star: pseudo-random pixels from a hash of the position.
fn write_noise_frame(path: &Path) {
    write_fits_with(path, 320, 240, &[], |x, y| {
        let mut h = (x as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ (y as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
        h = (h ^ (h >> 29)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        h ^= h >> 32;
        1000.0 + ((h % 61) as f64 - 30.0) * 1.5
    });
}

#[test]
fn a_frame_with_a_verdict_prints_dashes_and_sorts_last() {
    let dir = scratch("analyze-verdict");
    write_noise_frame(&dir.join("flat.fits"));
    write_star_field(&dir.join("sharp.fits"), 1.2, &[]);
    let opts = FrameOptions::default();
    let (_, _, s) = lib_measure(&dir.join("flat.fits"), &opts);
    assert_eq!(s.verdict, Some(Verdict::NoStars), "the noise frame has no star ({s:?})");

    let (_, rows) = analyze_csv(&dir, &["flat.fits", "sharp.fits"]);
    assert_row_eq(&rows[0], &lib_row("flat.fits", &dir.join("flat.fits"), &opts, "moffat4"), "flat");
    let col = |name: &str| CSV_COLUMNS.iter().position(|c| *c == name).unwrap();
    assert_eq!(rows[0][col("frame_verdict")], "no_stars");
    assert_eq!(rows[0][col("median_fwhm")], "0", "the CSV keeps what Athenaeum stores");
    assert_eq!(rows[1][col("frame_verdict")], "");

    // Table: dashes for the medians of a frame without a statistic, the verdict in words.
    let o = run(&dir, &["analyze", "flat.fits", "sharp.fits", "--quiet"]);
    assert!(o.status.success(), "{}", stderr(&o));
    let table = stdout(&o);
    let lines: Vec<Vec<&str>> = table.lines().map(|l| l.split_whitespace().collect()).collect();
    assert_eq!(lines[0], ["FILE", "STARS", "FWHM", "HFR", "ECC", "SNR", "WEIGHT", "BG", "NOISE", "VERDICT"]);
    let flat = &lines[1];
    assert_eq!(flat[0], "flat.fits");
    assert_eq!(flat[2..=6], ["-", "-", "-", "-", "-"], "FWHM HFR ECC SNR WEIGHT: {table}");
    assert_eq!(flat[9..].join(" "), "no stars", "{table}");
    let sharp = &lines[2];
    assert_eq!(sharp[0], "sharp.fits");
    assert!(sharp[2..=6].iter().all(|c| *c != "-"), "{table}");
    assert_eq!(sharp[9..], ["-"], "no verdict: {table}");

    // With optics the arcsec columns dash out too.
    let o = run(&dir, &["analyze", "flat.fits", "sharp.fits", "--quiet", "--focal-length", "448", "--pixel-size", "3.76"]);
    let table = stdout(&o);
    let flat: Vec<&str> = table.lines().nth(1).unwrap().split_whitespace().collect();
    assert_eq!(flat[2..=8], ["-"; 7], "FWHM HFR FWHM\" HFR\" ECC SNR WEIGHT: {table}");

    // --sort fwhm: a frame with a verdict is never ranked best, with either input order.
    for order in [["flat.fits", "sharp.fits"], ["sharp.fits", "flat.fits"]] {
        let (_, rows) = analyze_csv(&dir, &[order[0], order[1], "--sort", "fwhm"]);
        assert_eq!((rows[0][0].as_str(), rows[1][0].as_str()), ("sharp.fits", "flat.fits"), "{order:?}");
    }
}

#[test]
fn old_flags_still_run_with_a_warning() {
    let dir = scratch("analyze-old-flags");
    write_star_field(&dir.join("sharp.fits"), 1.2, &[]);

    let o = run(&dir, &["analyze", "sharp.fits", "--sigma", "8", "--measure-cap", "0", "--no-debayer", "--csv", "--quiet"]);
    assert_eq!(o.status.code(), Some(0), "{}", stderr(&o));
    let new = run(&dir, &["analyze", "sharp.fits", "--detection-threshold", "8", "--csv", "--quiet"]);
    assert!(new.status.success(), "{}", stderr(&new));
    assert_eq!(stdout(&o), stdout(&new), "--sigma 8 is --detection-threshold 8");
    assert!(new.stderr.is_empty(), "the new spelling is silent: {}", stderr(&new));

    let err = stderr(&o);
    let lines: Vec<&str> = err.lines().collect();
    assert_eq!(
        lines,
        [
            "warning: --sigma is --detection-threshold since 1.8.0: another detector, default 6.0",
            "warning: --measure-cap has no effect since 1.8.0: the star budget is --max-stars",
            "warning: --no-debayer has no effect since 1.8.0: a mosaic is measured on its green sites",
        ],
        "{err}"
    );

    // The old and the new spelling of the same knob are one knob.
    let both = run(&dir, &["analyze", "sharp.fits", "--sigma", "8", "--detection-threshold", "8", "--quiet"]);
    assert_eq!(both.status.code(), Some(2), "clap usage error: {}", stderr(&both));

    // The hidden flags stay out of --help.
    let help = stdout(&run(&dir, &["analyze", "--help"]));
    assert!(help.contains("--detection-threshold") && help.contains("--psf-model") && help.contains("--max-stars"), "{help}");
    assert!(!help.contains("--sigma") && !help.contains("--measure-cap") && !help.contains("--no-debayer"), "{help}");
}

#[test]
fn out_of_range_options_are_refused() {
    let dir = scratch("analyze-ranges");
    write_star_field(&dir.join("sharp.fits"), 1.2, &[]);
    for (args, message) in [
        (["--max-stars", "5"], "must be between 10 and 2000"),
        (["--max-stars", "2001"], "must be between 10 and 2000"),
        (["--detection-threshold", "1.5"], "must be between 2 and 20"),
        (["--detection-threshold", "25"], "must be between 2 and 20"),
        (["--sigma", "25"], "must be between 2 and 20"),
    ] {
        let o = run(&dir, &["analyze", "sharp.fits", "--csv", "--quiet", args[0], args[1]]);
        assert_eq!(o.status.code(), Some(2), "{args:?}: {}", stderr(&o));
        assert!(stderr(&o).contains(message), "{args:?}: {}", stderr(&o));
        assert!(o.stdout.is_empty(), "{args:?}: nothing on stdout");
    }
}

/// G5-f, the claim of port phase 5: on the real frames the CLI prints `analysis_mode`'s numbers, every column at
/// full precision. Real 26 MP frames are slow in a test run: `cargo test --release --test cli -- --ignored`.
#[test]
#[ignore = "26 MP frames: run on request, in release mode"]
fn analyze_matches_the_library_on_the_test_frames() {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let opts = FrameOptions::default();
    let mut checked = 0;
    for name in ["mono.fits", "osc.fits", "cocoon.fits", "test.xisf"] {
        let path = tests.join(name);
        if !path.exists() {
            // cocoon.fits is not tracked (the other three are); say so rather than pass silently.
            eprintln!("Skipping: tests/{name} not found");
            continue;
        }
        let (_, rows) = analyze_csv(&tests, &[name]);
        assert_row_eq(&rows[0], &lib_row(name, &path, &opts, "moffat4"), name);
        checked += 1;
    }
    assert!(checked >= 3, "mono.fits, osc.fits and test.xisf are tracked: {checked} frames checked");
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
