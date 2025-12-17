Name:           rustafits
Version:        0.2.1
Release:        1%{?dist}
Summary:        High-performance FITS/XISF to JPEG converter with auto-stretch

License:        GPLv3
URL:            https://github.com/eg013ra1n/rustafits
Source0:        %{url}/archive/v%{version}/%{name}-%{version}.tar.gz

BuildRequires:  rust-packaging >= 21
BuildRequires:  cargo
BuildRequires:  rust
BuildRequires:  pkgconfig
BuildRequires:  lz4-devel
BuildRequires:  libzstd-devel
BuildRequires:  zlib-devel

Requires:       lz4
Requires:       libzstd
Requires:       zlib

%description
High-performance FITS/XISF to JPEG converter for astronomical images with
QuickFits/PixInsight-compatible auto-stretch and Bayer debayering.

Features:
- FITS and XISF format support (including compressed XISF)
- Auto-Stretch: Median-based statistical stretching (robust to outliers)
- Bayer Debayering: Super-pixel 2×2 block averaging (RGGB, BGGR, GBRG, GRBG)
- Preview Mode: 2×2 binning for mono images (4× fewer pixels, ~90ms processing)
- Fast: 60-100ms for 4096×4096 images (Quickselect + SIMD + compiler optimizations)
- Multi-Platform SIMD: ARM NEON + x86_64 SSE2 with automatic selection

%prep
%autosetup -n %{name}-%{version}
%cargo_prep

%build
export RUSTAFITS_PORTABLE=1
%cargo_build

%install
export RUSTAFITS_PORTABLE=1
%cargo_install

%check
%cargo_test

%files
%license LICENSE
%doc README.md
%{_bindir}/%{name}

%changelog
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
