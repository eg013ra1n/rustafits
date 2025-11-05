// Re-export the C-based implementation
mod lib_c;
pub use lib_c::FitsConverter;

// Re-export modules for library users
pub mod fits;
pub mod processor;
pub mod output;
pub mod debayer;
pub mod stretch;
pub mod downscale;
