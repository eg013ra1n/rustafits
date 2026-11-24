//! Old header rules against new, over files (starloom phase 2, Task 2, ruling P2-R23).
//!
//! ```text
//! cargo run --release --features test-support --example header_rules_diff -- \
//!     [--all] [--write-pinned <dir>] [--list <file>]... FILES_OR_DIRS...
//! ```
//!
//! Every file is read through the crate's older copy of Athenaeum's header rules (`formats::stored_keys`, as
//! `star_measure::raw`, `formats::fits` and `formats::xisf` apply them today; `astroimage::old_header_rules`) and
//! through the new `astroimage::header` (Athenaeum's rules at ff7e0ca3). For `BITPIX`, `BZERO`, `BSCALE`, `DATAMAX`,
//! `EGAIN`, `XBAYROFF`, `YBAYROFF` and the folded Bayer pattern (`BAYER`), one line per key the two read differently:
//!
//! ```text
//! <file>\t<key>\told=<v>\tnew=<v>
//! ```
//!
//! and last `files: N, differing files: K, differing keys: M`. A value is `none` when absent or unreadable, an
//! integer as written, a number in `f64`'s round-trip form. The new values: `HeaderKeyMap::number` for the numbers,
//! `HeaderKeyMap::integer` for FITS `BITPIX` and for the offsets (absent is 0, as the readers fold), over
//! `header::read_header_key_map`; the new pattern is the reader's own pattern refolded with the new offsets. A file
//! one side reads and the other cannot is a `read` line; a file neither reads is an `unreadable` line, not counted
//! as differing.
//!
//! Two checks of the new module itself, printed before the summary and not counted as differences:
//! `text-map:<key>` lines (`bytes=` against `text=`) where `header_key_map(header_text(f))` answers otherwise than
//! the map built from the bytes, and `facts:<key>` lines (`map=` against `facts=`) where `read_frame_facts`' scanner
//! line reads an offset otherwise than `HeaderKeyMap::integer`.
//!
//! `--all` prints every key of every file, with a fifth column `same` or `differs`. A directory stands for its
//! `.fit`, `.fits`, `.fts`, `.xisf`, `.fz` and `.gz` files, recursively, sorted by path; `--list <file>` adds one path
//! per line. `--write-pinned <dir>` first writes every input the crate's header-rule tests pin as a file into
//! `<dir>` (and the rule-list probes into `<dir>/probes`), then exits.

use std::path::{Path, PathBuf};

use astroimage::header::{header_key_map, header_text, read_frame_facts, read_header_key_map, HeaderFormat};
use astroimage::old_header_rules::{self as old, KEYS};

fn usage() -> ! {
    eprintln!(
        "usage: header_rules_diff [--all] [--write-pinned <dir>] [--list <file>]... FILES_OR_DIRS..."
    );
    std::process::exit(2);
}

fn main() {
    let mut args = std::env::args().skip(1);
    let mut all = false;
    let mut inputs: Vec<PathBuf> = Vec::new();
    while let Some(a) = args.next() {
        match a.as_str() {
            "--all" => all = true,
            "--write-pinned" => {
                let dir = PathBuf::from(args.next().unwrap_or_else(|| usage()));
                let n = write_pinned(&dir).unwrap_or_else(|e| {
                    eprintln!("write-pinned: {e:#}");
                    std::process::exit(1)
                });
                eprintln!("wrote {n} files into {}", dir.display());
                return;
            }
            "--list" => {
                let list = args.next().unwrap_or_else(|| usage());
                let text = std::fs::read_to_string(&list).unwrap_or_else(|e| {
                    eprintln!("{list}: {e}");
                    std::process::exit(1)
                });
                inputs.extend(text.lines().filter(|l| !l.trim().is_empty()).map(PathBuf::from));
            }
            _ if a.starts_with("--") => usage(),
            _ => inputs.push(PathBuf::from(a)),
        }
    }
    if inputs.is_empty() {
        usage();
    }
    let mut files = Vec::new();
    for p in inputs {
        if p.is_dir() {
            let mut found = Vec::new();
            collect(&p, &mut found);
            found.sort();
            files.extend(found);
        } else {
            files.push(p);
        }
    }

    let mut totals = Totals::default();
    for f in &files {
        compare(f, all, &mut totals);
    }
    println!(
        "checks: unreadable by both: {}, text-map differing keys: {}, facts-vs-map differing keys: {}",
        totals.unreadable, totals.text_map, totals.facts
    );
    println!("files: {}, differing files: {}, differing keys: {}", files.len(), totals.files, totals.keys);
}

#[derive(Default)]
struct Totals {
    files: usize,
    keys: usize,
    unreadable: usize,
    text_map: usize,
    facts: usize,
}

fn collect(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        eprintln!("{}: cannot read the directory", dir.display());
        return;
    };
    for e in rd.flatten() {
        let p = e.path();
        if p.is_dir() {
            collect(&p, out);
        } else if p.extension().and_then(|x| x.to_str()).is_some_and(|x| {
            ["fit", "fits", "fts", "xisf", "fz", "gz"].iter().any(|k| x.eq_ignore_ascii_case(k))
        }) {
            out.push(p);
        }
    }
}

fn show(v: &Option<String>) -> &str {
    v.as_deref().unwrap_or("none")
}

fn one_line(e: &anyhow::Error) -> String {
    format!("{e:#}").replace(['\n', '\t'], " ")
}

fn odd(n: i64) -> bool {
    n.rem_euclid(2) == 1
}

/// The new rules' value of each of `KEYS`, and the new offsets.
fn new_values(map: &astroimage::header::HeaderKeyMap) -> (Vec<Option<String>>, (i64, i64)) {
    let n = |k: &str| old::number(map.number(k));
    let x = map.integer("XBAYROFF").unwrap_or(0);
    let y = map.integer("YBAYROFF").unwrap_or(0);
    let bitpix = match map.format() {
        HeaderFormat::Fits => old::int(map.integer("BITPIX")),
        HeaderFormat::Xisf => n("BITPIX"),
    };
    (
        vec![bitpix, n("BZERO"), n("BSCALE"), n("DATAMAX"), n("EGAIN"), old::int(Some(x)), old::int(Some(y))],
        (x, y),
    )
}

fn compare(f: &Path, all: bool, t: &mut Totals) {
    let file = f.display();
    let old_read = old::read(f);
    let new_map = read_header_key_map(f);
    let (o, m) = match (old_read, new_map) {
        (Ok(o), Ok(m)) => (o, m),
        (Err(a), Err(b)) => {
            t.unreadable += 1;
            println!("{file}\tunreadable\told=error: {}\tnew=error: {}", one_line(&a), one_line(&b));
            return;
        }
        (a, b) => {
            t.files += 1;
            t.keys += 1;
            let side = |r: Result<(), String>| r.map_or_else(|e| format!("error: {e}"), |_| "read".to_string());
            println!(
                "{file}\tread\told={}\tnew={}",
                side(a.map(|_| ()).map_err(|e| one_line(&e))),
                side(b.map(|_| ()).map_err(|e| one_line(&e)))
            );
            return;
        }
    };

    let (new_vals, (nx, ny)) = new_values(&m);
    let new_bayer = old::shift_bayer_pattern(o.bayer, odd(o.offsets.0) ^ odd(nx), odd(o.offsets.1) ^ odd(ny));
    let mut rows: Vec<(String, String, String)> = KEYS
        .iter()
        .zip(o.values.iter().zip(&new_vals))
        .map(|(k, (a, b))| (k.to_string(), show(a).to_string(), show(b).to_string()))
        .collect();
    rows.push(("BAYER".into(), format!("{:?}", o.bayer), format!("{new_bayer:?}")));
    if o.xisf != (m.format() == HeaderFormat::Xisf) {
        rows.push(("FORMAT".into(), if o.xisf { "xisf" } else { "fits" }.into(), format!("{:?}", m.format())));
    }
    let mut differing = 0;
    for (k, a, b) in &rows {
        let same = a == b;
        if !same {
            differing += 1;
        }
        if all {
            println!("{file}\t{k}\told={a}\tnew={b}\t{}", if same { "same" } else { "differs" });
        } else if !same {
            println!("{file}\t{k}\told={a}\tnew={b}");
        }
    }
    if differing > 0 {
        t.files += 1;
        t.keys += differing;
    }

    // The module's own checks: the stored text answers as the bytes do; the facts read the offsets as the map does.
    match header_text(f) {
        Ok((format, text)) => {
            let tm = header_key_map(format, &text);
            let (text_vals, _) = new_values(&tm);
            for ((k, a), b) in KEYS.iter().zip(&new_vals).zip(&text_vals) {
                if a != b {
                    t.text_map += 1;
                    println!("{file}\ttext-map:{k}\tbytes={}\ttext={}", show(a), show(b));
                }
            }
        }
        Err(e) => {
            t.text_map += 1;
            println!("{file}\ttext-map:read\tbytes=read\ttext=error: {}", one_line(&e));
        }
    }
    match read_frame_facts(f) {
        Ok(facts) => {
            for (k, fact) in [("XBAYROFF", facts.xbayroff), ("YBAYROFF", facts.ybayroff)] {
                let (a, b) = (old::int(m.integer(k)), old::int(fact));
                if a != b {
                    t.facts += 1;
                    println!("{file}\tfacts:{k}\tmap={}\tfacts={}", show(&a), show(&b));
                }
            }
        }
        Err(e) => {
            t.facts += 1;
            println!("{file}\tfacts:read\tmap=read\tfacts=error: {}", one_line(&e));
        }
    }
}

// ----- the pinned inputs, as files -----

/// An 80-column `KEY = value` card, the value right-aligned to column 30.
fn card(key: &str, value: &str) -> String {
    format!("{:<80}", format!("{key:<8}= {value:>20}"))
}

/// A 16-bit FITS image of `w` x `h` zero samples whose header is `cards` (each padded to 80 columns), then `END`.
fn write_fits(path: &Path, cards: &[String], w: usize, h: usize) -> std::io::Result<()> {
    let mut bytes: Vec<u8> = cards.iter().map(|c| format!("{c:<80}")).collect::<String>().into_bytes();
    bytes.extend(format!("{:<80}", "END").bytes());
    bytes.resize(bytes.len().div_ceil(2880) * 2880, b' ');
    bytes.extend(std::iter::repeat_n(0u8, w * h * 2));
    bytes.resize(bytes.len().div_ceil(2880) * 2880, 0);
    std::fs::write(path, bytes)
}

/// The structural cards of a 2 x 2 16-bit image, then `extra` as written.
fn fits_2x2(extra: &[&str]) -> Vec<String> {
    let mut c = vec![card("SIMPLE", "T"), card("BITPIX", "16"), card("NAXIS", "2"), card("NAXIS1", "2"), card("NAXIS2", "2")];
    c.extend(extra.iter().map(|s| s.to_string()));
    c
}

/// A 4 x 2 unsigned 16-bit top-down RGGB mosaic (`formats::fits`' `write_cfa_fits_raw`), then `extra` as written.
fn fits_cfa(extra: &[&str]) -> Vec<String> {
    let mut c = vec![
        card("SIMPLE", "T"),
        card("BITPIX", "16"),
        card("NAXIS", "2"),
        card("NAXIS1", "4"),
        card("NAXIS2", "2"),
        card("BZERO", "32768"),
        card("BAYERPAT", "'RGGB'"),
        card("ROWORDER", "'TOP-DOWN'"),
    ];
    c.extend(extra.iter().map(|s| s.to_string()));
    c
}

/// An XISF file: `prefix`, the XML declaration and an `<Image>` holding `children`, `suffix`, then the samples.
fn write_xisf(
    path: &Path,
    geometry: &str,
    sample_format: &str,
    children: &str,
    prefix: &str,
    suffix: &str,
) -> std::io::Result<()> {
    let parts: Vec<usize> = geometry.split(':').map(|p| p.parse().unwrap()).collect();
    let bps = if sample_format == "Float32" { 4 } else { 2 };
    let size = parts.iter().product::<usize>() * bps;
    let xml = |pos: usize| {
        format!(
            "{prefix}<?xml version=\"1.0\" encoding=\"UTF-8\"?><xisf version=\"1.0\">\
<Image geometry=\"{geometry}\" sampleFormat=\"{sample_format}\" colorSpace=\"Gray\" byteOrder=\"little\" \
location=\"attachment:{pos:08}:{size}\">{children}</Image></xisf>{suffix}"
        )
    };
    let len = xml(0).len();
    let mut bytes = b"XISF0100".to_vec();
    bytes.extend_from_slice(&(len as u32).to_le_bytes());
    bytes.extend_from_slice(&[0u8; 4]);
    bytes.extend_from_slice(xml(16 + len).as_bytes());
    bytes.extend(std::iter::repeat_n(0u8, size));
    std::fs::write(path, bytes)
}

fn kw(name: &str, value: &str) -> String {
    format!(r#"<FITSKeyword name="{name}" value="{value}" comment=""/>"#)
}

enum Pin {
    /// FITS cards.
    Fits(Vec<String>, usize, usize),
    /// XISF: geometry, sample format, children, prefix, suffix.
    Xisf(&'static str, &'static str, String, String, String),
}

fn x(geometry: &'static str, sample: &'static str, children: String) -> Pin {
    Pin::Xisf(geometry, sample, children, String::new(), String::new())
}

/// A 4 x 2 UInt16 mosaic with `BAYERPAT = 'RGGB'` first (`formats::xisf`'s `write_cfa_xisf`).
fn x_cfa(children: String) -> Pin {
    x("4:2:1", "UInt16", format!("{}{children}", kw("BAYERPAT", "'RGGB'")))
}

fn pinned() -> Vec<(&'static str, Pin)> {
    let f2 = |extra: &[&str]| Pin::Fits(fits_2x2(extra), 2, 2);
    let fc = |extra: &[&str]| Pin::Fits(fits_cfa(extra), 4, 2);
    let xf = |children: String| x("2:2:1", "Float32", children);
    let mut p: Vec<(&'static str, Pin)> = vec![
        // formats::stored_keys tests: xisf_rule_pins_athenaeums_results (XBAYROFF), on an RGGB mosaic.
        ("sk-x01-xbayroff-empty-then-1", x_cfa(format!("{}{}", kw("XBAYROFF", ""), kw("XBAYROFF", "1")))),
        ("sk-x02-xbayroff-quotepair-then-1", x_cfa(format!("{}{}", kw("XBAYROFF", "''"), kw("XBAYROFF", "1")))),
        ("sk-x03-xbayroff-novalue-then-1", x_cfa(format!(r#"<FITSKeyword name="XBAYROFF"/>{}"#, kw("XBAYROFF", "1")))),
        ("sk-x04-xbayroff-1-slash-off", x_cfa(kw("XBAYROFF", "1 / off"))),
        ("sk-x05-xbayroff-quoted1-slash-off", x_cfa(kw("XBAYROFF", "'1' / off"))),
        ("sk-x06-xbayroff-unclosed-quote", x_cfa(kw("XBAYROFF", "'1"))),
        ("sk-x07-xbayroff-blank-before-name", x_cfa(kw(" XBAYROFF", "1"))),
        ("sk-x08-xbayroff-lowercase-name", x_cfa(kw("xbayroff", "1"))),
        ("sk-x09-xbayroff-char-reference", x_cfa(kw("XBAYROFF", "&#49;"))),
        ("sk-x10-xbayroff-0-then-1", x_cfa(format!("{}{}", kw("XBAYROFF", "0"), kw("XBAYROFF", "1")))),
        ("sk-x11-xbayroff-start-end-pair", x_cfa(r#"<FITSKeyword name="XBAYROFF" value="1"></FITSKeyword>"#.into())),
        ("sk-x12-xbayroff-property", x_cfa(r#"<Property id="XBAYROFF" type="Int32" value="1"/>"#.into())),
        ("sk-x13-xbayroff-lone-quote", x_cfa(kw("XBAYROFF", "'"))),
        ("sk-x14-xbayroff-quote1quote0quote", x_cfa(kw("XBAYROFF", "'1'0'"))),
        (
            "sk-x15-xbayroff-before-the-xml-cut",
            Pin::Xisf("4:2:1", "UInt16", kw("BAYERPAT", "'RGGB'"), kw("XBAYROFF", "1"), String::new()),
        ),
        // xisf_num_and_offset_read_as_athenaeum_reads_them, on the keys the readers read.
        ("sk-x16-xbayroff-0p6", x_cfa(kw("XBAYROFF", "0.6"))),
        ("sk-x17-xbayroff-minus1", x_cfa(kw("XBAYROFF", "-1"))),
        ("sk-x18-xbayroff-nan", x_cfa(kw("XBAYROFF", "nan"))),
        ("sk-x19-egain-quoted-1p5", xf(kw("EGAIN", "'1.5'"))),
        // an_xisf_without_keywords_reads_no_keys, a_keyword_after_the_closing_tag_is_not_read,
        // the_parse_abort_log_line_carries_the_path_when_given, a_keyword_value_holding_an_element_is_read_as_xml_first,
        // dangling_ampersand_in_xisf_text_still_yields_keywords (with EGAIN as the keyword after the `&`).
        ("sk-x20-no-keywords", xf(String::new())),
        (
            "sk-x21-datamax-after-the-closing-tag",
            Pin::Xisf("2:2:1", "Float32", kw("EGAIN", "1.5"), "junk".into(), kw("DATAMAX", "9")),
        ),
        ("sk-x22-egain-after-a-comment-cut", xf(format!("<!-- </xisf> -->{}", kw("EGAIN", "2")))),
        ("sk-x23-xbayroff-in-a-value-element", x_cfa(kw("Z", "<FITSKeyword name='XBAYROFF' value='1'/>"))),
        (
            "sk-x24-egain-after-a-dangling-amp",
            xf(format!(r#"<Property id="note">Dark & Flat calibration</Property>{}"#, kw("EGAIN", "1.5"))),
        ),
        // fits_rule_pins_athenaeums_results (EGAIN).
        ("sk-f01-egain-no-equals-in-column-9", f2(&["EGAIN    1.5"])),
        ("sk-f02-egain-equals-no-space", f2(&["EGAIN   =1.5"])),
        ("sk-f03-egain-valueless-then-1p5", f2(&["EGAIN", "EGAIN   = 1.5"])),
        ("sk-f04-egain-hierarch", f2(&["HIERARCH EGAIN = 1.5"])),
        ("sk-f05-egain-blank-before-key", f2(&[" EGAIN  = 1.5"])),
        ("sk-f06-egain-string-amp", f2(&["EGAIN   = '1.5&'"])),
        ("sk-f07-egain-continue", f2(&["EGAIN   = '1.&'", "CONTINUE  '5'"])),
        ("sk-f08-egain-d-exponent-comment", f2(&["EGAIN   = 1.5D0 / e-/ADU"])),
        ("sk-f09-egain-repeated", f2(&["EGAIN   = 1.5", "EGAIN   = 2.5"])),
        ("sk-f10-egain-quoted", f2(&["EGAIN   = '1.5'"])),
        (
            "sk-f11-bitpix-16p0",
            Pin::Fits(
                vec![card("SIMPLE", "T"), "BITPIX  = 16.0".into(), card("NAXIS", "2"), card("NAXIS1", "2"), card("NAXIS2", "2")],
                2,
                2,
            ),
        ),
        // star_measure::raw: fits_facts_read_as_athenaeum_reads_them_on_unusual_cards.
        ("rf-f01-no-cards", f2(&[])),
        ("rf-f02-egain-no-equals", f2(&["EGAIN    1.5"])),
        ("rf-f03-datamax-no-equals", f2(&["DATAMAX   30000"])),
        ("rf-f04-egain-valueless-then-1p5", f2(&["EGAIN", "EGAIN   = 1.5"])),
        ("rf-f05-egain-hierarch", f2(&["HIERARCH EGAIN = 1.5"])),
        ("rf-f06-datamax-hierarch", f2(&["HIERARCH DATAMAX = 30000"])),
        ("rf-f07-egain-blank-before-key", f2(&[" EGAIN  = 1.5"])),
        ("rf-f08-egain-string-amp", f2(&["EGAIN   = '1.5&'"])),
        ("rf-f09-egain-continue", f2(&["EGAIN   = '1.&'", "CONTINUE  '5'"])),
        ("rf-f10-egain-repeated", f2(&["EGAIN   = 1.5", "EGAIN   = 2.5"])),
        // star_measure::raw: xisf_facts_read_as_athenaeum_reads_them_on_unusual_values.
        ("rx-01-egain-empty-then-1p5", xf(format!("{}{}", kw("EGAIN", ""), kw("EGAIN", "1.5")))),
        ("rx-02-egain-quotepair-then-1p5", xf(format!("{}{}", kw("EGAIN", "''"), kw("EGAIN", "1.5")))),
        ("rx-03-egain-1p5-slash-comment", xf(kw("EGAIN", "1.5 / e-/ADU"))),
        ("rx-04-egain-quoted-slash-comment", xf(kw("EGAIN", "'1.5' / e-/ADU"))),
        ("rx-05-egain-unclosed-quote", xf(kw("EGAIN", "'1.5"))),
        ("rx-06-egain-lowercase-blank-name", xf(kw(" egain", "1.5"))),
        ("rx-07-egain-char-reference", xf(kw("EGAIN", "1&#46;5"))),
        ("rx-08-egain-property", xf(r#"<Property id="EGAIN" type="Float32" value="1.5"/>"#.into())),
        ("rx-09-egain-start-end-pair", xf(r#"<FITSKeyword name="EGAIN" value="1.5"></FITSKeyword>"#.into())),
        ("rx-10-datamax-empty-then-0p9", xf(format!("{}{}", kw("DATAMAX", ""), kw("DATAMAX", "0.9")))),
        ("rx-11-datamax-slash-comment", xf(kw("DATAMAX", "0.9 / saturation"))),
        ("rx-12-datamax-unclosed-quote", xf(kw("DATAMAX", "'0.9"))),
        // star_measure::raw: an_xisf_mosaic_reads_unusual_offset_values_as_athenaeum_reads_them.
        ("rm-01-xbayroff-empty-then-1", x_cfa(format!("{}{}", kw("XBAYROFF", ""), kw("XBAYROFF", "1")))),
        ("rm-02-xbayroff-1-slash-off", x_cfa(kw("XBAYROFF", "1 / off"))),
        ("rm-03-xbayroff-unclosed-quote", x_cfa(kw("XBAYROFF", "'1"))),
    ];
    // formats::xisf: xisf_bayer_offset_values_read_as_athenaeums_reads_them.
    for (name, value) in [
        ("xv-01-xbayroff-1", "1"),
        ("xv-02-xbayroff-1p0", "1.0"),
        ("xv-03-xbayroff-quoted-1", "'1'"),
        ("xv-04-xbayroff-minus1", "-1"),
        ("xv-05-xbayroff-3", "3"),
        ("xv-06-xbayroff-0p6", "0.6"),
        ("xv-07-xbayroff-2", "2"),
        ("xv-08-xbayroff-x", "x"),
        ("xv-09-xbayroff-empty", ""),
        ("xv-10-xbayroff-nan", "nan"),
    ] {
        p.push((name, x_cfa(kw("XBAYROFF", value))));
    }
    p.extend([
        ("xv-11-no-offsets", x_cfa(String::new())),
        ("xv-12-ybayroff-0p6", x_cfa(kw("YBAYROFF", "0.6"))),
        ("xv-13-xbayroff-lowercase-name", x_cfa(kw("xbayroff", "1"))),
        ("xv-14-xbayroff-1-then-0", x_cfa(format!("{}{}", kw("XBAYROFF", "1"), kw("XBAYROFF", "0")))),
        ("xv-15-xbayroff-0-then-1", x_cfa(format!("{}{}", kw("XBAYROFF", "0"), kw("XBAYROFF", "1")))),
        // formats::xisf: xisf_bayer_offsets_follow_athenaeums_rule_on_unusual_values.
        ("xu-01-xbayroff-empty-then-1", x_cfa(format!("{}{}", kw("XBAYROFF", ""), kw("XBAYROFF", "1")))),
        ("xu-02-xbayroff-quotepair-then-1", x_cfa(format!("{}{}", kw("XBAYROFF", "''"), kw("XBAYROFF", "1")))),
        ("xu-03-xbayroff-novalue-then-1", x_cfa(format!(r#"<FITSKeyword name="XBAYROFF"/>{}"#, kw("XBAYROFF", "1")))),
        ("xu-04-xbayroff-1-slash-off", x_cfa(kw("XBAYROFF", "1 / off"))),
        ("xu-05-xbayroff-quoted1-slash-off", x_cfa(kw("XBAYROFF", "'1' / off"))),
        ("xu-06-ybayroff-1-slash-off", x_cfa(kw("YBAYROFF", "1 / off"))),
        ("xu-07-xbayroff-unclosed-quote", x_cfa(kw("XBAYROFF", "'1"))),
        ("xu-08-xbayroff-lone-quote", x_cfa(kw("XBAYROFF", "'"))),
        ("xu-09-xbayroff-blank-before-name", x_cfa(kw(" XBAYROFF", "1"))),
        ("xu-10-xbayroff-char-reference", x_cfa(kw("XBAYROFF", "&#49;"))),
        ("xu-11-xbayroff-start-end-pair", x_cfa(r#"<FITSKeyword name="XBAYROFF" value="1"></FITSKeyword>"#.into())),
        // formats::fits: quoted_bayer_offsets_fold_like_the_catalog.
        ("ff-01-xbayroff-quoted-1", fc(&["XBAYROFF= '1'"])),
        ("ff-02-ybayroff-quoted-1", fc(&["YBAYROFF= '1'"])),
        ("ff-03-both-quoted-1", fc(&["XBAYROFF= '1'", "YBAYROFF= '1'"])),
        ("ff-04-padded-and-commented", fc(&["XBAYROFF= '1       ' / padded and commented", "YBAYROFF= '0'"])),
        ("ff-05-quoted-2-and-3", fc(&["XBAYROFF= '2'", "YBAYROFF= '3'"])),
        // formats::fits: bayer_offsets_follow_the_catalog_rule_on_unusual_cards.
        ("ff-06-plain-commented", fc(&["XBAYROFF= 1 / Bayer pattern X axis offset", "YBAYROFF= 0 / Bayer pattern Y axis offset"])),
        ("ff-07-xbayroff-0-then-1", fc(&["XBAYROFF= 0", "XBAYROFF= 1"])),
        ("ff-08-xbayroff-1-then-0", fc(&["XBAYROFF= 1", "XBAYROFF= 0"])),
        ("ff-09-ybayroff-0-then-1-xbayroff-quoted", fc(&["YBAYROFF= 0", "YBAYROFF= 1", "XBAYROFF= '1'"])),
        ("ff-10-xbayroff-lowercase", fc(&["xbayroff= 1"])),
        ("ff-11-xbayroff-lowercase-blank", fc(&["xbayroff = 1"])),
        ("ff-12-ybayroff-mixed-case", fc(&["Ybayroff= 1"])),
        ("ff-13-xbayroff-1p0", fc(&["XBAYROFF= 1.0"])),
        ("ff-14-fractions", fc(&["XBAYROFF= 0.6", "YBAYROFF= 1.9"])),
        ("ff-15-xbayroff-0p1e1", fc(&["XBAYROFF= 0.1E1"])),
        ("ff-16-xbayroff-1d0", fc(&["XBAYROFF= 1D0"])),
        ("ff-17-signs", fc(&["XBAYROFF= -1", "YBAYROFF= +1"])),
        ("ff-18-commented-and-minus", fc(&["XBAYROFF= 1 / comment", "YBAYROFF= -1"])),
    ]);
    p
}

/// Inputs that no test pins, one per rule the step-1 list says changed, to check that list on files.
fn probes() -> Vec<(&'static str, Pin)> {
    let f2 = |extra: &[&str]| Pin::Fits(fits_2x2(extra), 2, 2);
    let fc = |extra: &[&str]| Pin::Fits(fits_cfa(extra), 4, 2);
    let xf = |children: String| x("2:2:1", "Float32", children);
    vec![
        ("pr-f01-hierarch-lowercase-egain", f2(&["HIERARCH egain = 1.5"])),
        ("pr-f02-hierarch-lowercase-then-plain", f2(&["HIERARCH egain = 1.5", "EGAIN   = 2.5"])),
        ("pr-f03-hierarch-lowercase-xbayroff", fc(&["HIERARCH xbayroff = 1"])),
        ("pr-x01-egain-repeated", xf(format!("{}{}", kw("EGAIN", "1.5"), kw("EGAIN", "2.5")))),
        ("pr-x02-egain-then-empty", xf(format!("{}{}", kw("EGAIN", "1.5"), kw("EGAIN", "")))),
        ("pr-x03-egain-d-exponent", xf(kw("EGAIN", "1.5D0"))),
        ("pr-x04-egain-trailing-quote", xf(kw("EGAIN", "1.5'"))),
        ("pr-x05-xbayroff-1p5", x_cfa(kw("XBAYROFF", "1.5"))),
        ("pr-x06-xbayroff-2p5", x_cfa(kw("XBAYROFF", "2.5"))),
        ("pr-x07-xbayroff-minus-0p6", x_cfa(kw("XBAYROFF", "-0.6"))),
        ("pr-x08-xbayroff-inf", x_cfa(kw("XBAYROFF", "inf"))),
        ("pr-x09-xbayroff-1-then-empty", x_cfa(format!("{}{}", kw("XBAYROFF", "1"), kw("XBAYROFF", "")))),
        ("pr-x10-xbayroff-case-variants", x_cfa(format!("{}{}", kw("XBAYROFF", "1"), kw("xbayroff", "0")))),
        ("pr-x11-xbayroff-1d0", x_cfa(kw("XBAYROFF", "1D0"))),
        ("pr-x12-datamax-quoted-padded", xf(kw("DATAMAX", "'  0.9  '"))),
    ]
}

fn write_pins(dir: &Path, pins: Vec<(&'static str, Pin)>) -> anyhow::Result<usize> {
    std::fs::create_dir_all(dir)?;
    let n = pins.len();
    for (name, pin) in pins {
        match pin {
            Pin::Fits(cards, w, h) => write_fits(&dir.join(format!("{name}.fits")), &cards, w, h)?,
            Pin::Xisf(g, s, children, prefix, suffix) => {
                write_xisf(&dir.join(format!("{name}.xisf")), g, s, &children, &prefix, &suffix)?
            }
        }
    }
    Ok(n)
}

fn write_pinned(dir: &Path) -> anyhow::Result<usize> {
    let a = write_pins(dir, pinned())?;
    let b = write_pins(&dir.join("probes"), probes())?;
    Ok(a + b)
}
