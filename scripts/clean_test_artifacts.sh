#!/usr/bin/env bash
# Remove the output artifacts that test binaries and IB apps leave in the repo root.
#
# Run manually after `./scripts/run_tests.sh` (or any test/app invocation) when you
# want a clean tree. The list below is the exhaustive set of files produced by
# src/tests/*.f90 and src/app/* with their default arguments — keep it in sync
# with the `save_into_file` / `export_*` calls in those sources.
#
# Explicit list (no `rm *.stl`) so we never delete files added by hand.

set -u
cd "$(dirname "$0")/.."

if [[ -t 1 ]]; then
  CYAN=$'\033[0;36m'; DIM=$'\033[2m'; RESET=$'\033[0m'
else
  CYAN=''; DIM=''; RESET=''
fi

# Fixed-name artifacts.
artifacts=(
  fossil_test_clip.stl
  fossil_test_clip_remainder.stl
  fossil_test_distance-aabb.dat
  fossil_test_distance-brute.dat
  fossil_test_distance_aabb_tree.dat
  fossil_test_merge.stl
  fossil_test_mirror.stl
  fossil_test_resize-factor-centroid.stl
  fossil_test_resize-factor.stl
  fossil_test_resize-xyz.stl
  fossil_test_rotate.stl
  fossil_test_translate-delta.stl
  fossil_test_translate-xyz.stl
  ib-cart-block-aabb.dat
  ib.grd
  ib.ib
  ib.vts
  ib_aabb_tree.dat
)

# Tree-dump outputs whose count depends on AABB refinement levels, plus
# valgrind profiler dumps from fossil_test_distance_bench (one per pid).
# Scoped to repo root with explicit prefixes so the glob cannot escape.
artifact_globs=(
  'fossil_test_distance_aabb-l_*-b_*.stl'
  'ibaabb-l_*-b_*.stl'
  'cachegrind.out.*'
  'callgrind.out.*'
)

removed=0

for f in "${artifacts[@]}"; do
  if [[ -e "$f" ]]; then
    rm -f -- "$f"
    printf "${DIM}removed${RESET}  %s\n" "$f"
    removed=$((removed + 1))
  fi
done

shopt -s nullglob
for g in "${artifact_globs[@]}"; do
  for f in $g; do
    rm -f -- "$f"
    printf "${DIM}removed${RESET}  %s\n" "$f"
    removed=$((removed + 1))
  done
done
shopt -u nullglob

if [[ $removed -eq 0 ]]; then
  printf "${CYAN}nothing to clean${RESET}\n"
else
  printf "\n${CYAN}cleaned %d artifact(s)${RESET}\n" "$removed"
fi
