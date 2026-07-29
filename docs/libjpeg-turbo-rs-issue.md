# 4:2:0 trailing-MCU padding diverges from C libjpeg-turbo (upstream issue #362)

**Status: RESOLVED.** Reported against `libjpeg-turbo-rs` 0.6.3 (also reproduced on
0.6.2). Upstream fixed it in **0.7.0** (commit `5d88fcf`) and ships a permanent
regression test, `tests/regression_issue_362_420_trailing_mcu.rs`. rustafits builds
on **0.8.0**, where the divergence is gone — verified locally on both architectures,
see "Verification on 0.8.0" below.

This file is kept for the record: it documents what we saw, what the real cause turned
out to be (not what we first thought), and how to re-check the behaviour if the pinned
dependency is ever bumped.

## Original report

At **4:2:0** only, `compress()` produced a different entropy-coded scan from C
libjpeg-turbo whenever `width % 16` fell in `1..=8` — that is, whenever the second
8×8 luma block column of the trailing MCU contains no real image samples at all.
Measured on `x86_64-unknown-linux-gnu` with AVX2, `simd` feature enabled (default),
against `turbojpeg` 1.5.1 (C libjpeg-turbo).

The decoded image was **bit-identical** either way, so this was a compatibility /
compression-ratio issue rather than a correctness bug. The cost was a consistent
**+0.7 % on the compressed size** for affected widths (synthetic noise; ~+0.02…0.12 %
on real photographic content).

4:4:4 and 4:2:2 matched C byte-for-byte at every width.

```
width sweep, height 320 — '=' means byte-identical to C libjpeg-turbo

  4:4:4  W=320..336  =================
  4:2:2  W=320..336  =================
  4:2:0  W=320..336  =########========
        w%16       01234567890123450
```

## Root cause — the report's diagnosis was wrong

We attributed the divergence to chroma: `expand_right_edge` running in RGB space
before `h2v2` downsampling rather than per-component after it. Upstream's fix commit
and regression test show that diagnosis is incorrect, and the difference matters for
anyone reading this later:

- The trailing MCU's second luma block column lies entirely outside the image and must
  be a **dummy** block — AC all zero, the previous block's DC copied in so the coded DC
  difference is zero (C `jccoefct.c`) — never the FDCT of replicated edge pixels.
- 4:2:0 chroma has exactly `ceil(width/16)` block columns, i.e. **never a padding
  column at all**, so no chroma-side padding strategy could have produced the window
  we measured.

The actual defect was in the **4:2:0 row fast path**: it guarded the last partial MCU
*row* but not the last partial MCU *column*, so it forward-transformed replicated edge
pixels where C emits a dummy block. `width % 16` in `1..=8` is exactly the residue set
that makes `ceil(width/8)` odd, which is why the reported window had that shape (same
defect as upstream issue #314).

The fix excludes a partial final MCU column from the fast path and lets the generic,
dummy-block-aware path encode it:

```rust
// libjpeg-turbo-rs 0.8.0, src/encode/pipeline.rs
let fast_cols: usize = if y_last_col_width == y_mcu_width { mcus_x } else { mcus_x - 1 };
…
generic_start_col = fast_cols;
```

0.6.3 had no such guard — it ran `for mcu_col in 0..mcus_x` unconditionally.

## Verification on 0.8.0

`repro.rs` (`libjpeg-turbo-rs-issue-repro.rs`, deps `libjpeg-turbo-rs` and
`turbojpeg = "1.5"` as the C reference) sweeps widths 320..=336 at height 320 for all
three subsamplings. Re-run 2026-07-29 on macOS (Apple Silicon), release build:

| build                                        | 4:4:4 | 4:2:2 | 4:2:0             |
| -------------------------------------------- | ----- | ----- | ----------------- |
| 0.8.0, `aarch64-apple-darwin` (NEON)         | `=`   | `=`   | `=` all widths    |
| 0.8.0, `x86_64-apple-darwin` (SSE2, Rosetta) | `=`   | `=`   | `=` all widths    |
| 0.6.3, `aarch64-apple-darwin` (NEON)         | `=`   | `=`   | `=` all widths    |
| 0.6.3, `x86_64-apple-darwin` (SSE2, Rosetta) | `=`   | `=`   | `#` **all** widths |

Two things worth keeping:

1. **aarch64 does not discriminate.** Both the broken and the fixed release are
   byte-clean on NEON, so a green test run on Apple Silicon proves nothing about this
   bug. Any re-check must include an x86_64 run.
2. **The x86_64 fault was wider than the filed report.** On this machine Rosetta
   exposes SSE2/SSE4.2 but **no AVX/AVX2**, so the AVX2 fast path never executes here.
   Under SSE2, 0.6.3 diverged from C at *every* width in the sweep (−13 % scan size),
   not only at `width % 16 ∈ 1..=8`. 0.8.0 is byte-exact with C at every width on both
   arches and both SIMD tiers.

The AVX2 code path itself was **not executed** in this verification (no AVX2-capable
execution environment locally: Rosetta and OrbStack's amd64 emulation both lack it).
It is covered by upstream CI and by the structure of the fix — with the `fast_cols`
guard the AVX2 path no longer touches the trailing partial column at all, and the
generic path that now owns that column is the one validated byte-exact above.

## Re-checking after a dependency bump

Two independent checks, in ascending cost:

1. `cargo test --test jpeg_encoder` and `cargo test --target x86_64-apple-darwin
   --test jpeg_encoder` — rustafits' own backend-agnostic contract tests.
2. Copy `libjpeg-turbo-rs-issue-repro.rs` into a scratch crate with the new version
   pinned and run it under **both** targets. Requires `cmake` + `nasm` for the
   `turbojpeg` reference only — rustafits itself needs neither.

Upstream's own `tests/regression_issue_362_420_trailing_mcu.rs` can also be run
directly from the vendored crate source; it byte-compares against stock `cjpeg`
(`brew install jpeg-turbo`) and additionally reads back coefficients to assert the
C dummy-block contract without needing any C toolchain.
