# 4:2:0 trailing-MCU padding diverges from C libjpeg-turbo when `width % 16` is in 1..=8

**Version:** `libjpeg-turbo-rs` 0.6.3 (also reproduces on 0.6.2), compared against
`turbojpeg` 1.5.1 (C libjpeg-turbo). Target: `x86_64-unknown-linux-gnu`, AVX2, `simd`
feature enabled (default).

## Summary

At **4:2:0** only, `compress()` produces a different entropy-coded scan from C
libjpeg-turbo whenever `width % 16` falls in `1..=8` — that is, whenever the second
8×8 luma block column of the trailing MCU contains no real image samples at all.

The decoded image is **bit-identical** either way, so this is a compatibility /
compression-ratio issue rather than a correctness bug. The cost is a consistent
**+0.7 % on the compressed size** for affected widths (measured on synthetic noise;
~+0.02…0.12 % on real photographic content).

4:4:4 and 4:2:2 match C byte-for-byte at every width, which points at the `h2v2`
(vertically subsampled) downsampling path specifically, not at horizontal
downsampling or at the FDCT.

## Reproduction

`repro.rs` is attached (single file, deps `libjpeg-turbo-rs = "0.6"` and
`turbojpeg = "1.5"`). Output on my machine:

```
width sweep, height 320 — '=' means byte-identical to C libjpeg-turbo

  4:4:4  W=320..336  =================
  4:2:2  W=320..336  =================
  4:2:0  W=320..336  =########========
        w%16       01234567890123450

the divergent output equals the C encoder fed a pre-padded image:

  W=322 (w%16= 2)  rust(W).scan == C(replicate-padded to 336).scan : true   [86132 B vs C native 85494 B, +0.746%]
  W=324 (w%16= 4)  rust(W).scan == C(replicate-padded to 336).scan : true   [86413 B vs C native 85791 B, +0.725%]
  W=326 (w%16= 6)  rust(W).scan == C(replicate-padded to 336).scan : true   [86718 B vs C native 86074 B, +0.748%]
  W=328 (w%16= 8)  rust(W).scan == C(replicate-padded to 336).scan : true   [86697 B vs C native 86030 B, +0.775%]
```

## Diagnosis

The second block of the reproduction pins the mechanism down. For an affected width
`W`, this holds exactly:

```
rust_encode(image, W)  ==  c_encode(replicate_last_column(image, W → ceil16(W)), ceil16(W))
```

byte-for-byte over the entropy-coded scan.

In other words, this crate behaves as if the right edge were expanded by column
replication **in RGB space, before color conversion and chroma downsampling**, out to
the full padded MCU width. C libjpeg-turbo does the expansion **per component, after
downsampling** (`expand_right_edge` runs on the already-downsampled component rows in
`compress_data`), so the synthetic chroma samples in the padding differ, and so do the
resulting quantized chroma coefficients.

When `width % 16` is in `9..=15` the trailing MCU's second luma block still contains
real samples, the two padding strategies coincide over the region that matters, and
the outputs agree — which is why the divergence window is exactly `1..=8`.

## Impact

Low, but worth flagging for two reasons:

1. **Bit-exactness with C libjpeg-turbo** is a reasonable expectation for a crate
   described as a reimplementation, and anyone diffing outputs or checksumming
   encoder results across a migration will hit this.
2. It leaves ~0.7 % of compressed size on the table for affected widths.

## Suggested fix

Move the right-edge expansion into the per-component path so it runs after chroma
downsampling, matching `expand_right_edge`'s position in libjpeg's `compress_data`,
rather than expanding the interleaved source image up front.

## Note

Not a complaint — the crate did exactly what I needed (dropping cmake+nasm from an
astrophotography tool's build, on both x86_64 and aarch64, at parity encode speed).
Filing this only because the divergence is cheap to reproduce and cheap to close.
