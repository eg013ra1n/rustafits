# fpack fixtures

Tile-compressed images written by the reference `fpack` (CFITSIO 4.6.4), for
`tests/fpack_fixtures.rs`. Synthetic data: a 64 × 64 unsigned 16-bit frame
(stored as `int16` with `BZERO = 32768`), background about 1000, five
Gaussian stars with flat saturated cores at 65535.

Copied from [OrbitalCommons/fits-test-cases](https://github.com/OrbitalCommons/fits-test-cases),
`tile-compression/`, commit `b13951b7fdecd78e913c958c7eb5da6a292b3f5f`, which
also holds each file's provenance note (`*.fits.md`), the generator
(`make_fpack_fixtures.sh`, `fpack_originals.c`) and the reference decodes
(`*.expected`, little-endian `float64`). Licence: Apache-2.0.

| File | Bytes | Source | Written by |
|------|------:|--------|------------|
| `fpack_rice_u16.fz` | 14,400 | `fpack_rice_u16.fits` | `fpack -r` (RICE_1, lossless) |
| `fpack_hc_u16_s4.fz` | 8,640 | `fpack_hc_u16_s4.fits` | `fpack -h -s 4` (HCOMPRESS_1, lossy) |
| `fpack_plio_u16.fz` | 17,280 | `fpack_plio_u16.fits` | `fpack -p` (PLIO_1) |
| `u16_stars.fits` | 11,520 | `fpack_rice_u16.expected` | the uncompressed original, stored plain |
| `fpack_hc_u16_s4.expected.fits` | 11,520 | `fpack_hc_u16_s4.expected` | the reference decode of the HCOMPRESS file, stored plain |

The `.fz` files are byte-for-byte copies. The two plain files hold the
`.expected` values as unsigned 16-bit samples (`BZERO = 32768`), every value
an exact integer in range. RICE_1 and PLIO_1 are lossless, so their
`.expected` files are identical and equal the original frame. The lossy
HCOMPRESS decode overshoots 65535 next to the saturated cores; the reference
clips those values to 65535 (213 pixels), and its largest error against the
original is 1558.
