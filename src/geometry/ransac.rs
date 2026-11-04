//! The seeded generator of Athenaeum's `geometry::ransac` (at 0b6fb183), for the test fixtures that draw from it
//! (`crate::test_support`). The RANSAC itself stays in Athenaeum.

/// SplitMix64 — tiny, deterministic, good enough to draw samples.
pub struct SplitMix64(pub u64);

impl SplitMix64 {
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    #[inline]
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n.max(1) as u64) as usize
    }
}
