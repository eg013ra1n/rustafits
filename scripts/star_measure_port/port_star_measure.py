"""Copy Athenaeum's star_measure into rustafits (port plan phase 2). Edits: paths, serde gates, helper visibility,
three Athenaeum-bound tests removed, `adapter` / `testing` declared, the provenance paragraph and the clippy
allow-list added. Nothing else.

    python3 scripts/star_measure_port/port_star_measure.py <athenaeum>/crates/athenaeum-core/src/star_measure \
        src/star_measure

Regenerates every copied file in place; `adapter.rs` and `testing.rs` are rustafits' own and are not touched. Run it
from a checkout of the Athenaeum sha the module doc names (0b6fb183), or update that sha in the provenance lines.
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from serde_cfg import transform

SRC, DST = sys.argv[1], sys.argv[2]

PATHS = [
    ('crate::integration::stats', 'crate::stats'),
    ('astroimage::', 'crate::'),
    ('crate::stacking::solve::tests::captured_logs', 'crate::test_support::captured_logs'),
    ('crate::stacking::telemetry::capture_logs', 'crate::test_support::capture_logs'),
]

# (file, exact item head) made `pub` and hidden from the docs: what Athenaeum calls outside the module (plan D-3).
HELPERS = [
    ('mod.rs', 'pub(crate) fn install_in_span<'),
    ('mod.rs', 'pub(crate) fn seeds(&self)'),
    ('mod.rs', 'pub(crate) fn verdict_time(&self)'),
    ('mod.rs', 'pub(crate) fn is_storable(&self)'),
    ('mod.rs', 'pub(crate) fn verdict_iterations(&self)'),
    ('mod.rs', 'pub(crate) fn for_pin(starless: bool)'),
    ('frame.rs', 'pub(crate) fn finite_median('),
    ('fit.rs', 'pub(crate) fn inside_half_max('),
    ('fit.rs', 'pub(crate) fn mark_blend_neighbours('),
    ('seeds.rs', 'pub(crate) fn envelope_floor_factor('),
    ('seeds.rs', 'pub(crate) fn read_noise_correlation('),
    ('seeds.rs', 'pub(crate) struct NoiseCorrelation'),
    ('seeds.rs', 'pub(crate) fn rho(&self)'),
]

# `#[cfg(test)]` fixtures Athenaeum's tests call: also under `test-support`.
TEST_SUPPORT = [
    ('mod.rs', 'pub(crate) fn for_pin(starless: bool)'),
    ('seeds.rs', 'pub(crate) fn for_pin() -> NoiseCorrelation'),
    ('fit.rs', 'impl SeedPass {'),
]

# Tests that need Athenaeum's stacking adapter or integration engine: they stay there (plan 2.4, phase 3).
DROP_TESTS = [
    ('mod.rs', 'fn module_measure_equals_the_adapter_bit_for_bit()'),
    ('mod.rs', 'fn adapter_matched_filter_seeds_fit_exactly_as_the_module_matched_filter()'),
    ('robust.rs', 'fn integration_copies_agree_with_this_module()'),
]

def fn_span(L, head):
    """Line range [attrs .. closing brace] of the 4-space-indented test fn whose line contains `head`."""
    i = next(k for k, l in enumerate(L) if head in l)
    s = i
    while L[s - 1].strip().startswith('#[') or L[s - 1].strip().startswith('///'):
        s -= 1
    e = next(k for k in range(i + 1, len(L)) if L[k] == '    }')
    # swallow one blank line after
    if e + 1 < len(L) and L[e + 1].strip() == '':
        e += 1
    return s, e

for name in sorted(os.listdir(SRC)):
    text = open(os.path.join(SRC, name)).read()
    for a, b in PATHS:
        text = text.replace(a, b)
    text = transform(text)
    L = text.split('\n')
    for f, head in DROP_TESTS:
        if f == name:
            s, e = fn_span(L, head)
            del L[s:e + 1]
    # tests that (de)serialize run under `--features serde` only; diagnostics' test module holds nothing else
    if name == 'diagnostics.rs':
        k = L.index('mod tests {')
        assert L[k - 1] == '#[cfg(test)]'
        L[k - 1] = '#[cfg(all(test, feature = "serde"))]'
    i = 0
    while i < len(L):
        if name != 'diagnostics.rs' and re.search(r'\b(serde_json|postcard|xxhash_rust)::', L[i]):
            j = i
            while not re.match(r'\s{4}fn \w+', L[j]):
                j -= 1
            k = j
            while L[k - 1].strip().startswith('#['):
                k -= 1
            if L[k].strip() != '#[cfg(feature = "serde")]':
                L.insert(k, '    #[cfg(feature = "serde")]')
                i += 1
            # skip to the end of this fn
            i = next(m for m in range(i, len(L)) if L[m] == '    }')
        i += 1
    # only serde's `default = "not_measured"` calls it
    if name == 'fit.rs':
        k = next(k for k, l in enumerate(L) if l.startswith('fn not_measured()'))
        L.insert(k, '#[cfg(feature = "serde")]')
    # Athenaeum's own tests build the encoding-pin record (`PlanePrep::for_pin`, with its two inner fixtures): it is
    # compiled under `test-support` as well (D-3)
    for f, head in TEST_SUPPORT:
        if f != name:
            continue
        k = next(k for k, l in enumerate(L) if l.strip().startswith(head))
        j = next(j for j in range(k - 1, k - 4, -1) if L[j].strip() == '#[cfg(test)]')
        L[j] = L[j].replace('#[cfg(test)]', '#[cfg(any(test, feature = "test-support"))]')
        if head == 'impl SeedPass {':
            # its `for_test` is this crate's tests' alone
            indent = L[j][:len(L[j]) - len(L[j].lstrip())]
            L.insert(j + 1, indent + '#[cfg_attr(not(test), allow(dead_code))]')
    for f, head in HELPERS:
        if f != name:
            continue
        hits = [k for k, l in enumerate(L) if l.strip().startswith(head)]
        assert len(hits) == 1, (f, head, hits)
        k = hits[0]
        indent = L[k][:len(L[k]) - len(L[k].lstrip())]
        L[k] = L[k].replace('pub(crate) ', 'pub ', 1)
        L.insert(k, indent + '#[doc(hidden)]')
    if name == 'mod.rs':
        k = L.index('//! defect by defect against the synthetic truth lab.')
        L[k + 1:k + 1] = [
            '//!',
            '//! **In rustafits** (port plan `docs/superpowers/plans/2026-11-03-star-measure-port.md`, phase 2): copied from',
            "//! Athenaeum's `athenaeum-core` at 0b6fb183 with paths, serde gates (feature `serde`) and the visibility of the",
            "//! [`adapter`] helpers changed, nothing else. Its numbers are Athenaeum's stacking and Analysis numbers. The specs,",
            "//! plans, rulings and modules the docs cite (`stacking::…`, `analysis::…`, `R-…`, `PC-…`) are Athenaeum's.",
        ]
        k = L.index("//! plans, rulings and modules the docs cite (`stacking::…`, `analysis::…`, `R-…`, `PC-…`) are Athenaeum's.")
        L[k + 1:k + 1] = [
            '',
            "// Athenaeum's code as copied (port plan G-3): the clippy style lints it trips are allowed, never \"fixed\" in a way",
            '// that could move a bit.',
            '#![allow(',
            '    clippy::bool_assert_comparison,',
            '    clippy::doc_lazy_continuation,',
            '    clippy::manual_repeat_n,',
            '    clippy::needless_lifetimes,',
            '    clippy::needless_range_loop,',
            '    clippy::neg_cmp_op_on_partial_ord,',
            '    clippy::ptr_arg,',
            '    clippy::question_mark,',
            '    clippy::too_many_arguments,',
            '    clippy::type_complexity,',
            '    clippy::unnecessary_map_or',
            ')]',
        ]
        # the test's second assignment is a no-op (it equalizes the timings before comparing): harmless
        k = L.index('        direct.timings = stored.timings;')
        j = k
        while not L[j].startswith('    fn '):
            j -= 1
        L.insert(j, '    #[allow(clippy::almost_swapped)]')
        k = L.index('pub mod substeps;')
        L[k + 1:k + 1] = [
            '',
            '#[doc(hidden)]',
            'pub mod adapter;',
            '#[cfg(any(test, feature = "test-support"))]',
            'mod testing;',
        ]
    open(os.path.join(DST, name), 'w').write('\n'.join(L))
print('ok')
