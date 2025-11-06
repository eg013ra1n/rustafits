Name:           rustafits
Version:        0.1.0
Release:        1%{?dist}
Summary:        High-performance FITS to JPEG converter with auto-stretch

License:        GPLv3
URL:            https://github.com/eg013ra1n/rustafits
Source0:        %{url}/archive/v%{version}/%{name}-%{version}.tar.gz

BuildRequires:  rust-packaging >= 21
BuildRequires:  cargo
BuildRequires:  rust
BuildRequires:  cfitsio-devel
BuildRequires:  pkgconfig

Requires:       cfitsio

%description
High-performance FITS to JPEG converter for astronomical images with
QuickFits/PixInsight-compatible auto-stretch and Bayer debayering.

Features:
- Auto-Stretch: Median-based statistical stretching (robust to outliers)
- Bayer Debayering: Super-pixel 2×2 block averaging (RGGB, BGGR, GBRG, GRBG)
- Preview Mode: 2×2 binning for mono images (4× fewer pixels, ~90ms processing)
- Fast: 60-100ms for 4096×4096 images (Quickselect + SIMD + compiler optimizations)
- Multi-Platform SIMD: ARM NEON + x86_64 SSE2 with automatic selection

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
%doc README.md PERFORMANCE.md
%{_bindir}/%{name}

%changelog
* Wed Nov 06 2025 Vilen Sharifov <vilen.sharifov@gmail.com> - 0.1.0-1
- Initial package
- High-performance FITS to JPEG converter
- QuickFits/PixInsight-compatible auto-stretch
- Super-pixel Bayer debayering
- Multi-platform SIMD optimization (NEON/SSE2)
