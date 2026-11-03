# Plate Solving

`astroimage::platesolving` is a set of **building blocks**, not a solver. It
carries no star catalogue, does no catalogue I/O, and has no orchestration loop:
you bring detected image stars and catalogue stars, and this module gives you
quad matching, outlier rejection, transform fitting and a WCS you can write to a
FITS header. Athenaeum's solver, `solvemyastro`, does its own matching and uses
rustafits only to decode frames and detect stars (`ImageAnalyzer::detect_fast`).
Athenaeum itself uses this module's types (`WcsSolution`, `SipCoefficients`,
`SolveHints`, `ProperMotionCorrector`).

Everything here is `f64`. RA/Dec are degrees; tangent-plane coordinates (xi, eta)
are radians; pixel coordinates are pixels.

## The shape of a solve

```
detected stars ──► build_quads ──┐
                                 ├──► match_quads ──► fit_affine_from_centers
catalogue stars ─► build_quads ──┘         │                     │
   (projected to the tangent plane)        │                     ▼
                                           │            coarse AffineTransform
                                           ▼                     │
                                  per-star correspondences       │
                                           │                     │
                                           ▼                     ▼
                                    RansacFilter ──────► TransformFitter
                                                                 │
                                                                 ▼
                                                          WcsSolution
```

The two halves are deliberately separate: quad matching establishes *which* star
is which without knowing the plate scale, and the transform fit turns those
correspondences into geometry.

## Quad matching

A **quad** is one star plus its three nearest neighbours (`QUAD_SIZE = 4`). Its
six pairwise distances (`NUM_EDGES = C(4,2) = 6`) are sorted and normalised by the
longest, which leaves five scale-invariant ratios — the sixth is always 1.0 and is
not stored.

Matching is a brute-force comparison of those ratio vectors with an early exit on
the first mismatch. Two quads are similar when, for all five ratios, the smaller
divided by the larger is at least `1 - tolerance`:

```rust
let ratio = if a > b { b / a } else { a / b };
if ratio < 1.0 - tolerance { /* reject */ }
```

Each surviving pair yields a **center correspondence** — the quad's own star, not
one of the four in some order — so there is no star-ordering ambiguity to resolve.
`QuadMatch` also carries `scale_ratio` (catalogue longest distance over image
longest distance), which is the plate-scale estimate that falls out of the match
for free.

Three or more center correspondences are enough for
`fit_affine_from_centers` to produce a coarse `AffineTransform`
(`a1 b1 c1 / a2 b2 c2`, i.e. `x' = a1 x + b1 y + c1`).

### Densifying the pool on sparse fields

`build_quads` makes exactly one quad per star: star + 3 nearest neighbours. On a
sparse field that pool can be too small to find a match at all.

`build_quads_multi(stars, max_stars, group_size)` forms a local group of
`[center star, its group_size - 1 nearest neighbours]` and emits **every 4-subset
of that group**, so the pool grows combinatorially: `group_size = 5` gives
`C(5,4) = 5` quads per star, `group_size = 6` gives `C(6,4) = 15`. `build_quads`
is preserved exactly as the `group_size = 4` case, so switching to `_multi` with
`group_size = 4` changes nothing.

Two behaviours to know:

- `group_size` is clamped to **6** (the group buffer is fixed at 6), and raised to
  `QUAD_SIZE` if smaller — asking for 8 silently gives you 6.
- Quads are deduplicated by their index set, so subsets shared between
  neighbouring centre stars are emitted once.

The cost is quadratic in pool size on both sides — `match_quads` is a nested loop —
so densify deliberately.

### Which stars to feed it

Use the fast detector (`detect_fast`, backed by `detect_stars_adaptive`; see
[detection.md](detection.md)). Two of its outputs matter here:

- The **adaptive ladder** is what makes long-focal-length and nebula-swamped
  fields yield enough well-centroided stars to solve at all.
- **`FastStar::snr`** is the right ranking key for the quad pool. Extended
  structure — galaxy and nebula knots — carries huge flux spread over a large
  aperture, so it has low SNR and would otherwise dominate a flux-ranked pool and
  fill it with things that are not point sources.

`with_centroid_refine(true)` additionally refines each centroid with the Moffat
LM and keeps the pass-1 value in `raw_x`/`raw_y`, so a solver can fit from both
and keep whichever gives the better residual.

## Projection

`GnomonicProjection` is the TAN projection, used to bring catalogue RA/Dec onto
the flat plane the affine fit lives in:

| Function | Direction |
| ---- | ---- |
| `sky_to_tangent(ra, dec, ra0, dec0) -> (xi, eta)` | degrees in, radians out |
| `tangent_to_sky(xi, eta, ra0, dec0) -> (ra, dec)` | radians in, degrees out |

`(ra0, dec0)` is the tangent point — the projection centre, normally the field
centre hint.

## Outlier rejection

`RansacFilter::filter` takes matched pairs (image pixels ↔ catalogue tangent
plane) and finds the largest consistent set by random sampling plus affine
fitting. Defaults (`RansacConfig`):

| Field | Default | Meaning |
| ---- | ---- | ---- |
| `threshold_px` | `2.5` | Inlier residual ceiling, pixels |
| `max_iterations` | `100` | Sampling rounds |
| `min_inliers` | `6` | Below this the result is rejected outright |

## Transform fitting

`TransformFitter` turns surviving correspondences into geometry:

- `fit_wcs(...)` — least squares for the CD matrix:
  `xi_deg = CD1_1 * u + CD1_2 * v`, `eta_deg = CD2_1 * u + CD2_2 * v`, where
  `u = x - crpix_x` and `v = y - crpix_y`.
- SIP distortion — fit the linear CD matrix first, then the polynomial
  correction on top (`FitModel::Sip { order }`, coefficients in
  `SipCoefficients { order, coeffs: [[f64; 6]; 6] }`).
- `fit_pixel(...)` — a pixel-to-pixel `PixelTransform` with no sky coordinates at
  all, for registration rather than solving.

## The solution

`WcsSolution` holds `crpix`, `crval`, the `cd` matrix and optional forward and
reverse SIP coefficients, and answers the questions a caller actually asks:

| Method | Returns |
| ---- | ---- |
| `pixel_to_sky(x, y)` | `(ra, dec)` in degrees |
| `sky_to_pixel(ra, dec)` | `(x, y)` in pixels |
| `pixel_scale_arcsec()` | arcsec per pixel |
| `field_rotation_deg()` | field rotation |
| `to_fits_headers()` | `Vec<(String, String)>` ready to write into a header |

`ProperMotionCorrector::propagate` moves catalogue positions between epochs
before projection, which matters as soon as the catalogue epoch and the frame
epoch are years apart.

## Configuration

`PatternMatcherConfig` defaults:

| Field | Default | Meaning |
| ---- | ---- | ---- |
| `max_stars` | `100` | Stars considered when building quads |
| `hash_tolerance` | `0.01` | Ratio tolerance in `match_quads` |
| `scale_hint` | `None` | Expected arcsec/px — "reject quads whose scale differs >20%" |
| `multi_probe` | `true` | "Check adjacent hash bins (+-1 in each dimension)" |

**`PatternMatcherConfig` is a carrier, not a knob this crate reads.** Only
`hash_tolerance` and `max_stars` have counterparts here, and even those are passed
to `match_quads`/`build_quads` explicitly as arguments rather than read from the
struct. `scale_hint` and `multi_probe` are referenced nowhere in
`src/platesolving/` — and `multi_probe`'s own description mentions hash bins,
which this matcher does not have at all: `match_quads` is the brute-force ratio
comparison described above. Treat both fields as a contract with the consumer
(`solvemyastro`), not as documentation of behaviour implemented here.

`SolveHints` carries what the caller already knows — `ra`, `dec`, `fov_deg`,
`rotation`, `pixel_scale_arcsec`, all `Option` — and narrowing them is the
difference between a targeted solve and a blind one.

## Types at a glance

| Type | Fields |
| ---- | ---- |
| `ImageStar` | `x`, `y`, `flux` |
| `CatalogStar` | `ra`, `dec`, `mag` |
| `ProjectedStar` | `xi`, `eta`, `mag`, `ra`, `dec` |
| `StarMatch` | `image_idx`, `catalog_idx`, `residual_px` |
| `Quad` | `star_indices[4]`, `ratios[5]`, `center`, `longest_dist` |
| `QuadMatch` | `image_center`, `catalog_center`, `scale_ratio` |
| `AffineTransform` | `a1`, `b1`, `c1`, `a2`, `b2`, `c2` |
| `FitModel` | `Affine` / `Projective` / `Sip { order }` |
