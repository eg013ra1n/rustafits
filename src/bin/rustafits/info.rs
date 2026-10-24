//! `rustafits info` — header summary of many files without reading pixels.

use std::path::PathBuf;

use anyhow::Result;
use astroimage::{BayerPattern, FileFormat, ImageConverter, ImageHeader};
use clap::Args;

use crate::batch;
use crate::report::{json_str, print_csv, print_table};

#[derive(Args, Debug)]
pub struct InfoArgs {
    /// Input FITS/XISF files or directories
    #[arg(required = true, value_name = "INPUT")]
    paths: Vec<PathBuf>,

    /// Descend into subdirectories of directory inputs
    #[arg(short, long)]
    recursive: bool,

    /// Extra columns from header keywords (repeatable), e.g. -k FOCALLEN -k AIRMASS
    #[arg(short = 'k', long = "keyword", value_name = "KEY")]
    keywords: Vec<String>,

    /// Every keyword of every file instead of the summary table
    #[arg(long, conflicts_with_all = ["json", "csv"])]
    all: bool,

    /// JSON: one object per file with the layout and every keyword
    #[arg(long, conflicts_with = "csv")]
    json: bool,

    /// CSV of the summary columns (plus -k columns)
    #[arg(long)]
    csv: bool,
}

/// A summary column: header name and the keywords tried in order (FITS cards
/// as N.I.N.A./SharpCap/ASIAIR write them, then the XISF property id).
const COLUMNS: &[(&str, &[&str])] = &[
    ("EXP", &["EXPTIME", "EXPOSURE", "Instrument:ExposureTime"]),
    ("FILTER", &["FILTER", "Instrument:Filter:Name"]),
    ("GAIN", &["GAIN", "Instrument:Camera:Gain"]),
    ("TEMP", &["CCD-TEMP", "CCD_TEMP", "Instrument:Sensor:Temperature"]),
    ("DATE-OBS", &["DATE-OBS", "Observation:Time:Start"]),
    ("OBJECT", &["OBJECT", "Observation:Object:Name"]),
];

pub fn run(args: InfoArgs) -> Result<bool> {
    let inputs = batch::collect_inputs(&args.paths, args.recursive)?;
    let mut ok = true;
    let mut headers = Vec::new();
    for input in &inputs {
        match ImageConverter::read_header(&input.path) {
            Ok(h) => headers.push((input.path.display().to_string(), h)),
            Err(e) => {
                ok = false;
                eprintln!("error: {}: {e:#}", input.path.display());
            }
        }
    }

    if args.json {
        print_json(&headers);
    } else if args.all {
        print_all(&headers);
    } else {
        let mut cols: Vec<&str> = vec!["FILE", "SIZE", "TYPE", "BAYER"];
        cols.extend(COLUMNS.iter().map(|(name, _)| *name));
        cols.extend(args.keywords.iter().map(String::as_str));
        let rows: Vec<Vec<String>> = headers
            .iter()
            .map(|(file, h)| {
                let mut row = vec![file.clone(), size(h), h.sample_type.to_string(), bayer(h.bayer_pattern).to_string()];
                row.extend(COLUMNS.iter().map(|(_, keys)| first_of(h, keys)));
                row.extend(args.keywords.iter().map(|k| first_of(h, &[k.as_str()])));
                row
            })
            .collect();
        if args.csv {
            print_csv(&cols, &rows);
        } else {
            let right: Vec<bool> = cols.iter().map(|c| matches!(*c, "EXP" | "GAIN" | "TEMP")).collect();
            print_table(&cols, &rows, &right);
        }
    }
    Ok(ok)
}

fn first_of(h: &ImageHeader, keys: &[&str]) -> String {
    keys.iter().find_map(|k| h.get(k)).filter(|v| !v.is_empty()).unwrap_or("-").to_string()
}

fn size(h: &ImageHeader) -> String {
    if h.channels > 1 {
        format!("{}x{}x{}", h.width, h.height, h.channels)
    } else {
        format!("{}x{}", h.width, h.height)
    }
}

fn bayer(p: BayerPattern) -> &'static str {
    match p {
        BayerPattern::None => "-",
        BayerPattern::Rggb => "RGGB",
        BayerPattern::Bggr => "BGGR",
        BayerPattern::Gbrg => "GBRG",
        BayerPattern::Grbg => "GRBG",
    }
}

fn format_name(f: FileFormat) -> &'static str {
    match f {
        FileFormat::Fits => "FITS",
        FileFormat::Xisf => "XISF",
    }
}

fn print_all(headers: &[(String, ImageHeader)]) {
    for (i, (file, h)) in headers.iter().enumerate() {
        if i > 0 {
            println!();
        }
        println!(
            "== {file}  ({} {} {}, Bayer {}, {})",
            format_name(h.format),
            size(h),
            h.sample_type,
            bayer(h.bayer_pattern),
            if h.flip_vertical { "bottom-up" } else { "top-down" }
        );
        let width = h.keywords.iter().map(|k| k.key.chars().count()).max().unwrap_or(0);
        for k in &h.keywords {
            if k.comment.is_empty() {
                println!("{:<width$} = {}", k.key, k.value);
            } else {
                println!("{:<width$} = {}  / {}", k.key, k.value, k.comment);
            }
        }
    }
}

fn print_json(headers: &[(String, ImageHeader)]) {
    println!("[");
    for (i, (file, h)) in headers.iter().enumerate() {
        let keywords: Vec<String> = h
            .keywords
            .iter()
            .map(|k| {
                format!(
                    "{{\"key\": {}, \"value\": {}, \"comment\": {}}}",
                    json_str(&k.key),
                    json_str(&k.value),
                    json_str(&k.comment)
                )
            })
            .collect();
        let bayer = match h.bayer_pattern {
            BayerPattern::None => "null".to_string(),
            p => json_str(bayer(p)),
        };
        println!(
            "  {{\"file\": {}, \"format\": {}, \"width\": {}, \"height\": {}, \"channels\": {}, \
             \"sample_type\": {}, \"bayer\": {bayer}, \"flip_vertical\": {}, \"keywords\": [{}]}}{}",
            json_str(file),
            json_str(format_name(h.format)),
            h.width,
            h.height,
            h.channels,
            json_str(h.sample_type),
            h.flip_vertical,
            keywords.join(", "),
            if i + 1 < headers.len() { "," } else { "" }
        );
    }
    println!("]");
}
