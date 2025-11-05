fn main() {
    // Get cfitsio paths from pkg-config
    let cfitsio = pkg_config::Config::new()
        .probe("cfitsio")
        .expect("Failed to find cfitsio library");

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

    build
        .file("c_src/fits_processor.c")
        .file("c_src/jpeg_writer.c")
        .file("ffi/debayer.c")
        .file("ffi/stretch.c")
        .include("c_src")
        .include("ffi")
        .opt_level(3)  // Maximum optimization
        .flag("-march=native")  // Use native CPU instructions (enables NEON on ARM, SSE2 on x86)
        .flag("-ffast-math")  // Fast math optimizations
        .warnings(false);

    // Add cfitsio include paths
    for path in &cfitsio.include_paths {
        build.include(path);
    }

    build.compile("fits_processor");

    // Link CFITSIO
    println!("cargo:rustc-link-lib=cfitsio");
    println!("cargo:rustc-link-lib=m");  // Math library

    // Tell cargo to rerun if C sources change
    println!("cargo:rerun-if-changed=c_src/fits_processor.c");
    println!("cargo:rerun-if-changed=c_src/fits_processor.h");
    println!("cargo:rerun-if-changed=c_src/jpeg_writer.c");
    println!("cargo:rerun-if-changed=ffi/debayer.c");
    println!("cargo:rerun-if-changed=ffi/debayer.h");
    println!("cargo:rerun-if-changed=ffi/stretch.c");
    println!("cargo:rerun-if-changed=ffi/stretch.h");

    // Try to find cfitsio with pkg-config
    if let Ok(cfitsio) = pkg_config::Config::new().probe("cfitsio") {
        for path in cfitsio.include_paths {
            println!("cargo:include={}", path.display());
        }
    } else {
        // Fallback paths for macOS Homebrew
        println!("cargo:rustc-link-search=/opt/homebrew/lib");
        println!("cargo:rustc-link-search=/usr/local/lib");
    }
}
