Name:           rustafits
Version:        0.8.0
Release:        1%{?dist}
Summary:        High-performance FITS/XISF to JPEG converter with auto-stretch

License:        Apache-2.0
URL:            https://github.com/eg013ra1n/rustafits
Source0:        %{url}/archive/v%{version}/%{name}-%{version}.tar.gz

BuildRequires:  rust-packaging >= 21
BuildRequires:  cargo
BuildRequires:  rust

%description
High-performance FITS/XISF to JPEG converter for astronomical images with
PixInsight STF-compatible auto-stretch and Bayer debayering.

Features:
- FITS and XISF format support (including compressed XISF)
- Auto-Stretch: Median-based statistical stretching (robust to outliers)
- Bayer Debayering: Super-pixel 2x2 block averaging (RGGB, BGGR, GBRG, GRBG)
- Preview Mode: 2x2 binning for mono images (4x fewer pixels, ~90ms processing)
- Fast: 60-100ms for 4096x4096 images (Quickselect + SIMD + compiler optimizations)
- Multi-Platform SIMD: ARM NEON + x86_64 SSE2 with automatic selection
- Pure Rust: No C dependencies, fully self-contained static binary

%prep
%autosetup -n %{name}-%{version}
%cargo_prep

%build
%cargo_build

%install
%cargo_install

%check
%cargo_test

%files
%license LICENSE
%doc README.md
%{_bindir}/%{name}

%changelog
* Mon Mar 24 2026 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.8.0-1
- Add detection-stage morphological filters (sharpness, concentration index, edge margin)
- Add moments-based pre-screening gate before LM fitting
- Add spatial grid selection for balanced star measurement
- Lower default measure_cap from 2000 to 500

* Thu Mar 12 2026 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.7.5-1
- Fix duplicate star detections: replace Voronoi deblending with multi-peak skip
- Fix NMS boundary off-by-one in star detection
- Add post-detection centroid dedup safety net

* Tue Mar 11 2026 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.7.3-1
- Add fit-residual-weighted statistics for improved PixInsight agreement
- Add two-stage trail detection (Rayleigh + PSF-fit eccentricity)
- Add trail-aware statistics (bypass ecc filter on trailed frames)
- Add LmResult struct with fit_residual from LM solvers
- Fix Rayleigh test: add R² floor, proper even-length median, raise min stars to 20
- Update all documentation for new pipeline features

* Sun Mar 08 2026 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.6.4-1
- Add MRS wavelet noise estimation (B3-spline à trous transform)
- Add fixed-beta Moffat PSF fitting (7 free parameters)
- Add configurable distortion pre-filter (with_max_distortion)
- Add compare subcommand to debug CLI for bulk metric comparison
- Update documentation for all new algorithms and settings

* Fri Mar 06 2026 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.6.0-1
- Add 2D elliptical Moffat PSF fitting with automatic Gaussian fallback
- Add separable matched-filter convolution for star detection
- Add iterative source-masked background estimation
- Fix annotation ellipse orientation when optimizer converges with swapped axes
- Fix annotation theta not negated on vertical flip (FITS Y-up images)
- Canonicalize PSF fit output: fwhm_x >= fwhm_y, theta along major axis
- Improve background mesh estimation with bilinear interpolation

* Sun Mar 01 2026 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.5.5-1
- Add per-star PSF position angle (theta) to StarMetrics
- Fix analysis pipeline panics on real-world FITS images
- Add dual-path Rayleigh trail rejection for under- and oversampled stars
- Make trail detection advisory: expose R² and possibly_trailed on AnalysisResult
- Add green-channel interpolation and green-pixel-only PSF fitting for OSC
- Add star annotation overlay with 3-tier API and CLI support
- Add peak-based deblending for crowded-field star detection
- Add single-read path for annotated image generation
- Fix deblending regression on extended objects

* Thu Feb 26 2026 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.4.5-1
- Add RGBA output support (with_rgba_output, channels field on ProcessedImage)
- SSSE3 gray-to-RGB on x86_64 (replaces scalar fallback)
- Parallel downscale via rayon
- NEON-accelerated gray-to-RGBA
- Generalized vertical flip (RGB and RGBA)
- RGBA-aware PNG/JPEG output (auto-strips alpha for JPEG)

* Wed Feb 26 2026 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.4.3-1
- Fix downscale ordering: debayer before downscale for Bayer/OSC images
- Bayer-aware downscale: super-pixel debayer counts as 2x reduction

* Mon Feb 10 2025 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.4.2-1
- Add configurable rayon thread pool for multi-image concurrent processing
- Re-export rayon ThreadPool and ThreadPoolBuilder from library API

* Sun Feb 09 2025 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.3.0-1
- Rewrite all C code to pure Rust
- Remove system library dependencies (zlib, lz4, zstd)
- Use pure Rust compression crates (flate2, lz4_flex, ruzstd)
- Use image crate for JPEG/PNG output (replaces stb_image_write)
- Use quick-xml for XISF parsing (replaces fragile strstr-based parser)
- Fully self-contained static binary with zero external native dependencies

* Tue Dec 17 2024 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.2.1-1
- Fix portable builds for package managers (Homebrew, RPM)
- Add RUSTAFITS_PORTABLE environment variable for portable binaries

* Wed Nov 06 2024 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.2.0-1
- Add XISF format support with LZ4, zstd, and zlib compression
- SIMD optimizations for XISF processing

* Wed Nov 06 2024 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.1.0-1
- Initial package
- High-performance FITS to JPEG converter
- PixInsight STF-compatible auto-stretch
- Super-pixel Bayer debayering
- Multi-platform SIMD optimization (NEON/SSE2)
