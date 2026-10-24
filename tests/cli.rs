//! End-to-end tests of the `rustafits` binary on small synthetic FITS files
//! written by the tests themselves (no fixtures needed).
#![cfg(feature = "cli")]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

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
    let mut cards = vec![
        "SIMPLE  =                    T".to_string(),
        "BITPIX  =                   16".to_string(),
        "NAXIS   =                    2".to_string(),
        format!("NAXIS1  = {w:>20}"),
        format!("NAXIS2  = {h:>20}"),
        "BZERO   =                32768".to_string(),
        "BSCALE  =                    1".to_string(),
    ];
    if let Some(p) = bayer {
        cards.push(format!("BAYERPAT= '{p}'"));
    }
    cards.push("END".to_string());
    let mut header: String = cards.iter().map(|c| format!("{c:<80}")).collect();
    while header.len() % 2880 != 0 {
        header.push(' ');
    }
    let mut bytes = header.into_bytes();
    for y in 0..h {
        for x in 0..w {
            let noise = ((x * 7919 + y * 104_729) % 97) as i64;
            let star = if (x % 17 == 5) && (y % 13 == 4) { 20_000 } else { 0 };
            let v = (1500 + (x + y) as i64 * 3 + noise + star).min(65_535);
            bytes.extend_from_slice(&((v - 32_768) as i16).to_be_bytes());
        }
    }
    while bytes.len() % 2880 != 0 {
        bytes.push(0);
    }
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
