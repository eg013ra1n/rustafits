class Rustafits < Formula
  desc "High-performance FITS/XISF to JPEG/PNG converter with auto-stretch"
  homepage "https://github.com/eg013ra1n/rustafits"
  url "https://github.com/eg013ra1n/rustafits/archive/v0.4.5.tar.gz"
  sha256 "CHECKSUM_HERE"
  license "Apache-2.0"

  depends_on "rust" => :build

  def install
    system "cargo", "install", *std_cargo_args
  end

  test do
    assert_match "FITS/XISF to JPEG Converter", shell_output("#{bin}/rustafits --help", 1)
  end
end
