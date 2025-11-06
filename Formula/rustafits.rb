class Rustafits < Formula
  desc "High-performance FITS to JPEG converter with auto-stretch"
  homepage "https://github.com/eg013ra1n/rustafits"
  url "https://github.com/eg013ra1n/rustafits/archive/v0.1.0.tar.gz"
  sha256 "CHECKSUM_HERE"
  license "GPL-3.0"

  depends_on "rust" => :build
  depends_on "pkg-config" => :build
  depends_on "cfitsio"

  def install
    system "cargo", "install", *std_cargo_args
  end

  test do
    assert_match "FITS to JPEG Converter", shell_output("#{bin}/rustafits --help", 1)
  end
end
