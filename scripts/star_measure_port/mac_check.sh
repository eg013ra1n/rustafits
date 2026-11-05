#!/usr/bin/env bash
# Port plan phases 3 and 4: the gates that need the owner's Mac (aarch64), in one run.
#
#   1. Digest: `port_digest.rs` built against Athenaeum's `main` (its own star_measure, or the re-export of rustafits')
#      and against the branch under test (rustafits' star_measure through the re-export) must write byte-identical
#      files.
#   2. Suites with the exact aarch64 pins: `cargo test -p rustafits --features serde` and
#      `cargo test -p athenaeum-core --all-targets`, both on the branch.
#
# The lab scorecard (athenaeum-lab `star_probe`, before = `main`, after = the branch, 0 moved) is the third gate; run
# it with the lab's own harness, this script only reminds you.
#
# Usage:
#   [MAC_CHECK_BRANCH=<branch>] [MAC_CHECK_MAIN=<ref>] [MAC_CHECK_REMOTE=<remote>] \
#       scripts/star_measure_port/mac_check.sh <athenaeum checkout> [work dir]
#
#   MAC_CHECK_BRANCH  the Athenaeum branch under test (default star-measure-port-phase-3; phase 4 passes
#                     star-measure-port-phase-4). Its worktree is <work dir>/<branch>.
#   MAC_CHECK_MAIN    the "before" side, in <work dir>/main. `main` (default) is the fetched $MAC_CHECK_REMOTE/main; any
#                     other value is a ref or sha used as given, so a `main` that has moved does not change the
#                     "before" (phase 4's gate passes MAC_CHECK_MAIN=e6344ddd).
#   MAC_CHECK_REMOTE  where `main` and the branch are fetched from (default `origin`). In the owner's checkout `origin`
#                     is GitLab: push the branch there first, or set MAC_CHECK_REMOTE=github. The branch's rustafits
#                     submodule commit must be reachable from the submodule's own remote (`git submodule update`).
#
# The work dir (default ~/.athenaeum-acc/star-measure-port) holds two worktrees and one cargo target; give it a
# REGISTRY.tsv row and remove it when done (`git -C <athenaeum> worktree remove …`, then the directory).
set -euo pipefail

if [ $# -lt 1 ]; then
    sed -n '2,/^set -euo/p' "$0" | sed '$d'
    exit 2
fi
ATH=$(cd "$1" && pwd)
WORK=${2:-$HOME/.athenaeum-acc/star-measure-port}
HERE=$(cd "$(dirname "$0")" && pwd)
RUSTAFITS=$(cd "$HERE/../.." && pwd)
FRAMES=$RUSTAFITS/tests
BRANCH=${MAC_CHECK_BRANCH:-star-measure-port-phase-3}
MAIN=${MAC_CHECK_MAIN:-main}
REMOTE=${MAC_CHECK_REMOTE:-origin}

mkdir -p "$WORK"
export CARGO_TARGET_DIR=$WORK/target
git -C "$ATH" fetch -q "$REMOTE" main "$BRANCH"

for side in main "$BRANCH"; do
    wt=$WORK/$side
    if [ ! -d "$wt" ]; then
        if [ "$side" = main ]; then
            # `main` is the fetched remote main; any other MAC_CHECK_MAIN (a sha, a tag) is used as given.
            if [ "$MAIN" = main ]; then
                git -C "$ATH" worktree add -q --detach "$wt" "$REMOTE/main"
            else
                git -C "$ATH" worktree add -q --detach "$wt" "$MAIN"
            fi
        else
            git -C "$ATH" worktree add -q --detach "$wt" "$REMOTE/$side"
        fi
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
