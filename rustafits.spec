Name:           rustafits
Version:        0.3.0
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
QuickFits/PixInsight-compatible auto-stretch and Bayer debayering.

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
- QuickFits/PixInsight-compatible auto-stretch
- Super-pixel Bayer debayering
- Multi-platform SIMD optimization (NEON/SSE2)
