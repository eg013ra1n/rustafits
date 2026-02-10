# Maintainer: Vilen Sharifov <vilen.sharifov@gmail.com>
pkgname=rustafits
pkgver=0.1.0
pkgrel=1
pkgdesc="High-performance FITS to JPEG converter with auto-stretch and Bayer debayering"
arch=('x86_64' 'aarch64')
url="https://github.com/eg013ra1n/rustafits"
license=('Apache-2.0')
depends=('cfitsio')
makedepends=('rust' 'cargo' 'pkg-config')
source=("$pkgname-$pkgver.tar.gz::https://github.com/eg013ra1n/$pkgname/archive/v$pkgver.tar.gz")
sha256sums=('SKIP')  # Update with actual checksum after first release

prepare() {
    cd "$pkgname-$pkgver"
    cargo fetch --locked --target "$CARCH-unknown-linux-gnu"
}

build() {
    cd "$pkgname-$pkgver"
    export RUSTUP_TOOLCHAIN=stable
    export CARGO_TARGET_DIR=target
    cargo build --release --locked
}

check() {
    cd "$pkgname-$pkgver"
    cargo test --release --locked
}

package() {
    cd "$pkgname-$pkgver"
    install -Dm755 "target/release/$pkgname" "$pkgdir/usr/bin/$pkgname"
    install -Dm644 README.md "$pkgdir/usr/share/doc/$pkgname/README.md"
    install -Dm644 PERFORMANCE.md "$pkgdir/usr/share/doc/$pkgname/PERFORMANCE.md"
}
