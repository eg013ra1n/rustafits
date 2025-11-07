# Building macOS Binaries Manually

Since GitHub Actions charges for macOS runners, macOS binaries need to be built and uploaded manually.

## Prerequisites

```bash
brew install cfitsio
```

## Build for Your Architecture

### Apple Silicon (M1/M2/M3/M4)

```bash
cargo build --release
cd target/release
tar czf rustafits-macos-aarch64.tar.gz rustafits
shasum -a 256 rustafits-macos-aarch64.tar.gz > rustafits-macos-aarch64.tar.gz.sha256
```

### Intel Mac (x86_64)

If you're on Apple Silicon but want to build for Intel:

```bash
rustup target add x86_64-apple-darwin
cargo build --release --target x86_64-apple-darwin
cd target/x86_64-apple-darwin/release
tar czf rustafits-macos-x86_64.tar.gz rustafits
shasum -a 256 rustafits-macos-x86_64.tar.gz > rustafits-macos-x86_64.tar.gz.sha256
```

If you're on Intel Mac:

```bash
cargo build --release
cd target/release
tar czf rustafits-macos-x86_64.tar.gz rustafits
shasum -a 256 rustafits-macos-x86_64.tar.gz > rustafits-macos-x86_64.tar.gz.sha256
```

## Upload to GitHub Release

1. Go to: https://github.com/eg013ra1n/rustafits/releases
2. Edit the release
3. Drag and drop:
   - `rustafits-macos-aarch64.tar.gz`
   - `rustafits-macos-aarch64.tar.gz.sha256`
   - `rustafits-macos-x86_64.tar.gz` (if you built it)
   - `rustafits-macos-x86_64.tar.gz.sha256` (if you built it)
4. Update release notes to mention macOS binaries

## Alternative: Use GitHub CLI

```bash
# Build first (see above)
cd target/release  # or target/x86_64-apple-darwin/release

# Upload to release
gh release upload v0.1.0 rustafits-macos-aarch64.tar.gz rustafits-macos-aarch64.tar.gz.sha256
```

## For Future Releases

If you enable GitHub Actions billing for macOS runners, update `.github/workflows/release.yml` to include:

```yaml
- os: macos-latest
  target: aarch64-apple-darwin
  asset_name: rustafits-macos-aarch64
```
