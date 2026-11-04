#!/usr/bin/env bash
# Port plan phase 3: the gates that need the owner's Mac (aarch64), in one run.
#
#   1. Digest: `port_digest.rs` built against Athenaeum's `main` (its own star_measure) and against the phase-3
#      branch (rustafits' star_measure through the re-export) must write byte-identical files.
#   2. Suites with the exact aarch64 pins: `cargo test -p rustafits --features serde` and
#      `cargo test -p athenaeum-core --all-targets`, both on the phase-3 branch.
#
# The lab scorecard (athenaeum-lab `star_probe`, before = `main`, after = the branch, 0 moved) is the third gate; run
# it with the lab's own harness, this script only reminds you.
#
# Usage:
#   scripts/star_measure_port/mac_check.sh <athenaeum checkout> [work dir]
#
# The work dir (default ~/.athenaeum-acc/star-measure-port) holds two worktrees and one cargo target; give it a
# REGISTRY.tsv row and remove it when done (`git -C <athenaeum> worktree remove …`, then the directory).
#
# The branches come from the remote MAC_CHECK_REMOTE (default `origin`). In the owner's checkout `origin` is GitLab:
# push the branch there first, or set MAC_CHECK_REMOTE=github.
set -euo pipefail

if [ $# -lt 1 ]; then
    sed -n '2,20p' "$0"
    exit 2
fi
ATH=$(cd "$1" && pwd)
WORK=${2:-$HOME/.athenaeum-acc/star-measure-port}
HERE=$(cd "$(dirname "$0")" && pwd)
RUSTAFITS=$(cd "$HERE/../.." && pwd)
FRAMES=$RUSTAFITS/tests
BRANCH=star-measure-port-phase-3
REMOTE=${MAC_CHECK_REMOTE:-origin}

mkdir -p "$WORK"
export CARGO_TARGET_DIR=$WORK/target
git -C "$ATH" fetch -q "$REMOTE" main "$BRANCH"

for side in main "$BRANCH"; do
    wt=$WORK/$side
    if [ ! -d "$wt" ]; then
        git -C "$ATH" worktree add -q --detach "$wt" "$REMOTE/$side"
    fi
    git -C "$wt" submodule update -q --init rustafits solvemyastro
    sed 's/use SM_PATH as sm;/use athenaeum_core::star_measure as sm;/' "$HERE/port_digest.rs" \
        > "$wt/crates/athenaeum-core/examples/port_digest.rs"
    echo "== $side: $(git -C "$wt" log --oneline -1), rustafits $(git -C "$wt/rustafits" log --oneline -1 | cut -c1-60)"
    (cd "$wt" && cargo build -q -p athenaeum-core --example port_digest)
    "$CARGO_TARGET_DIR/debug/examples/port_digest" "$WORK/digest-$side.txt" "$FRAMES" 2>/dev/null
    rm "$wt/crates/athenaeum-core/examples/port_digest.rs"
done

echo
echo "== 1. digest ($(uname -m))"
shasum -a 256 "$WORK/digest-main.txt" "$WORK/digest-$BRANCH.txt"
if cmp -s "$WORK/digest-main.txt" "$WORK/digest-$BRANCH.txt"; then
    echo "PASS: byte-identical"
else
    echo "FAIL: the digests differ (first difference below)"
    cmp "$WORK/digest-main.txt" "$WORK/digest-$BRANCH.txt" || true
fi

echo
echo "== 2. suites on $BRANCH ($(uname -m), exact pins on aarch64)"
cd "$WORK/$BRANCH"
cargo test -q -p rustafits --features serde 2>&1 | grep -E "^test result|FAILED|panicked" || true
cargo test -q -p athenaeum-core --all-targets 2>&1 | grep -E "^test result|FAILED|panicked" || true

echo
echo "== 3. lab scorecard: run star_probe before (main) and after ($BRANCH) with the lab harness; the bar is 0 moved."
