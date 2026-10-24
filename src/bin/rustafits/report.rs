//! Plain-text tables, JSON and CSV for `info` and `analyze` (no serde: the
//! records are flat and the escaping rules are short).

/// Columns padded to their widest cell; `right[i]` right-aligns column i.
pub fn print_table(headers: &[&str], rows: &[Vec<String>], right: &[bool]) {
    let mut widths: Vec<usize> = headers.iter().map(|h| h.chars().count()).collect();
    for row in rows {
        for (w, cell) in widths.iter_mut().zip(row) {
            *w = (*w).max(cell.chars().count());
        }
    }
    let line = |cells: Vec<&str>| {
        let mut out = String::new();
        for (i, cell) in cells.iter().enumerate() {
            let pad = widths[i].saturating_sub(cell.chars().count());
            let last = i + 1 == cells.len();
            if right.get(i).copied().unwrap_or(false) {
                out.push_str(&" ".repeat(pad));
                out.push_str(cell);
            } else {
                out.push_str(cell);
                if !last {
                    out.push_str(&" ".repeat(pad));
                }
            }
            if !last {
                out.push_str("  ");
            }
        }
        println!("{out}");
    };
    line(headers.to_vec());
    for row in rows {
        line(row.iter().map(String::as_str).collect());
    }
}

pub fn json_str(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// A JSON number, or `null` for NaN/±inf (which JSON cannot carry).
pub fn json_num(v: f64) -> String {
    if v.is_finite() {
        format!("{v}")
    } else {
        "null".to_string()
    }
}

/// A CSV field, quoted when it holds a comma, quote or line break.
pub fn csv_field(s: &str) -> String {
    if s.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

pub fn print_csv(headers: &[&str], rows: &[Vec<String>]) {
    println!("{}", headers.iter().map(|h| csv_field(h)).collect::<Vec<_>>().join(","));
    for row in rows {
        println!("{}", row.iter().map(|c| csv_field(c)).collect::<Vec<_>>().join(","));
    }
}
