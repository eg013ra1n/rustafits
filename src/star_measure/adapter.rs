//! The helpers Athenaeum's stacking adapter and Analysis layer call beyond this module's API (port plan D-3): `pub`
//! for it, hidden from the docs and outside semver. Each is reachable at its own path as well; this list is the
//! inventory. Methods in the same position: [`PlanePrep::seeds`](super::PlanePrep::seeds),
//! [`verdict_time`](super::PlanePrep::verdict_time), [`verdict_iterations`](super::PlanePrep::verdict_iterations),
//! [`is_storable`](super::PlanePrep::is_storable), [`NoiseCorrelation::rho`], and `PlanePrep::for_pin` (the
//! encoding-pin record; compiled for tests and under `test-support`).

pub use super::fit::{inside_half_max, mark_blend_neighbours};
pub use super::frame::finite_median;
pub use super::install_in_span;
pub use super::seeds::{envelope_floor_factor, read_noise_correlation, NoiseCorrelation};
