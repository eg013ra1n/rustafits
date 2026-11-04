//! Geometry helpers [`crate::star_measure`] shares with Athenaeum, copied verbatim from Athenaeum's `geometry`
//! (at 0b6fb183). Crate-private: Athenaeum keeps its own copies for registration.

// The module uses `build` and `within`; the rest of the tree's API stays as in Athenaeum.
#[allow(dead_code, clippy::manual_is_multiple_of)]
pub mod kdtree;
#[cfg(test)]
#[allow(dead_code)]
pub mod ransac;
