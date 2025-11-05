fn main() {
    // Get cfitsio paths from pkg-config
    let cfitsio = pkg_config::Config::new()
        .probe("cfitsio")
        .expect("Failed to find cfitsio library");

    // Compile C code
    let mut build = cc::Build::new();
    build
        .file("c_src/fits_processor.c")
        .file("c_src/jpeg_writer.c")
        .file("ffi/debayer.c")
        .file("ffi/stretch.c")
        .include("c_src")
        .include("ffi")
        .opt_level(3)  // Maximum optimization
        .flag("-march=native")  // Use native CPU instructions
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
