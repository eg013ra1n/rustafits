fn main() {
    // Compile C code
    let mut build = cc::Build::new();

    // Try to use Homebrew LLVM on macOS for OpenMP support
    #[cfg(target_os = "macos")]
    {
        if let Ok(llvm_path) = std::env::var("HOMEBREW_PREFIX") {
            let llvm_bin = format!("{}/opt/llvm/bin/clang", llvm_path);
            if std::path::Path::new(&llvm_bin).exists() {
                build.compiler(&llvm_bin);
                println!("cargo:rustc-link-search={}/opt/llvm/lib", llvm_path);
            }
        } else {
            // Fallback to common Homebrew paths
            for prefix in &["/opt/homebrew", "/usr/local"] {
                let llvm_bin = format!("{}/opt/llvm/bin/clang", prefix);
                if std::path::Path::new(&llvm_bin).exists() {
                    build.compiler(&llvm_bin);
                    println!("cargo:rustc-link-search={}/opt/llvm/lib", prefix);
                    break;
                }
            }
        }
    }

    // Get compression library paths from pkg-config (for XISF support)
    let zlib = pkg_config::Config::new().probe("zlib").ok();
    let lz4 = pkg_config::Config::new().probe("liblz4").ok();
    let zstd = pkg_config::Config::new().probe("libzstd").ok();

    build
        .file("c_src/fits_processor.c")
        .file("c_src/fits_reader.c")
        .file("c_src/jpeg_writer.c")
        .file("c_src/xisf_reader.c")
        .file("c_src/base64.c")
        .include("c_src")
        .opt_level(3)  // Maximum optimization
        .flag("-ffast-math")  // Fast math optimizations
        .flag("-ftree-vectorize")  // Enable auto-vectorization
        .flag("-funroll-loops")  // Unroll loops for better performance
        .flag("-fomit-frame-pointer");  // Free up a register

    // Use native CPU optimizations only for local builds (not package managers)
    // Set RUSTAFITS_PORTABLE=1 to build portable binaries (e.g., for Homebrew bottles)
    if std::env::var("RUSTAFITS_PORTABLE").is_err() {
        build.flag("-march=native");  // Use native CPU instructions
        build.flag("-mtune=native");  // Tune for native CPU
    }

    build.warnings(false);

    // Add compression library include paths
    if let Some(ref z) = zlib {
        for path in &z.include_paths {
            build.include(path);
        }
    }
    if let Some(ref l) = lz4 {
        for path in &l.include_paths {
            build.include(path);
        }
    }
    if let Some(ref z) = zstd {
        for path in &z.include_paths {
            build.include(path);
        }
    }

    build.compile("fits_processor");

    // Link math library
    println!("cargo:rustc-link-lib=m");

    // Link compression libraries for XISF support
    println!("cargo:rustc-link-lib=z");      // zlib
    println!("cargo:rustc-link-lib=lz4");    // LZ4
    println!("cargo:rustc-link-lib=zstd");   // Zstandard

    // Tell cargo to rerun if C sources change
    println!("cargo:rerun-if-changed=c_src/fits_processor.c");
    println!("cargo:rerun-if-changed=c_src/fits_processor.h");
    println!("cargo:rerun-if-changed=c_src/fits_reader.c");
    println!("cargo:rerun-if-changed=c_src/fits_reader.h");
    println!("cargo:rerun-if-changed=c_src/jpeg_writer.c");
    println!("cargo:rerun-if-changed=c_src/xisf_reader.c");
    println!("cargo:rerun-if-changed=c_src/xisf_reader.h");
    println!("cargo:rerun-if-changed=c_src/base64.c");
}
