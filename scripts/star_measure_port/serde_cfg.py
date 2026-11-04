"""Shared transform: put serde derives / attributes / imports behind `feature = "serde"`, drop ts_rs::TS."""
import re

def transform(src: str) -> str:
    out = []
    for line in src.split('\n'):
        stripped = line.strip()
        indent = line[:len(line) - len(line.lstrip())]
        if stripped == 'use serde::{Deserialize, Serialize};':
            out.append(indent + '#[cfg(feature = "serde")]')
            out.append(line)
            continue
        m = re.fullmatch(r'#\[derive\((.*)\)\]', stripped)
        if m and ('Serialize' in m.group(1) or 'ts_rs::TS' in m.group(1)):
            items = [t.strip() for t in m.group(1).split(',')]
            has_serde = 'Serialize' in items or 'Deserialize' in items
            rest = [t for t in items if t not in ('Serialize', 'Deserialize', 'ts_rs::TS')]
            if rest:
                out.append(f'{indent}#[derive({", ".join(rest)})]')
            if has_serde:
                out.append(f'{indent}#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]')
            continue
        m = re.fullmatch(r'#\[serde\((.*)\)\]', stripped)
        if m:
            out.append(f'{indent}#[cfg_attr(feature = "serde", serde({m.group(1)}))]')
            continue
        out.append(line)
    return '\n'.join(out)
