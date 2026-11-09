#!/usr/bin/env bash
# Performance baselines (issue #94): the bench src/tests/validation/fossil_val_perf.f90, built by
#   fobis build --mode validation-gnu        (-O2, no OpenMP, into exe-validation/)
# prints one line per kernel: group, kernel, problem size, items, t_min and t_median [s] of N repetitions after one
# warm-up, t_min per item, and a 64-bit hash of the bit patterns of the kernel's outputs.
#
# Usage:
#   scripts/run_perf.sh NAME [reps=N]   run the bench on one thread and write src/tests/validation/perf/NAME.txt
#                                       (header: commit, host, compiler, flags, date). An existing NAME is refused:
#                                       baselines are never overwritten.
#   scripts/run_perf.sh --compare A B   per kernel: t_min(A), t_min(B), speed-up A/B, hash status. Exit 1 if a hash
#                                       differs or a kernel of A is missing in B (or the reverse).
# FOSSIL_PERF_COMMIT overrides the commit recorded in the header (default: git rev-parse, +dirty if modified).
# A and B are names in src/tests/validation/perf/ or paths to files. Baselines are machine-specific: compare runs
# made on one machine (the header says where each was made).
set -uo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"
dir=src/tests/validation/perf

resolve() { if [ -f "$1" ]; then echo "$1"; else echo "$dir/$1.txt"; fi; }

if [ "${1:-}" = "--compare" ]; then
  [ $# -eq 3 ] || { echo "usage: $0 --compare A B"; exit 2; }
  a=$(resolve "$2"); b=$(resolve "$3")
  for f in "$a" "$b"; do [ -f "$f" ] || { echo "no baseline $f"; exit 2; }; done
  awk -v A="$2" -v B="$3" '
    # kernel lines: key = group/kernel; t_min is the 4th field from the end, the hash the last
    FNR == 1 { file++ }
    /^#/ || /^Are all tests passed/ || NF < 8 { next }
    { key = $1 "/" $2
      if (file == 1) { order[++n] = key; ta[key] = $(NF - 3); ha[key] = $NF }
      else { tb[key] = $(NF - 3); hb[key] = $NF; inb[key] = 1 } }
    END {
      printf "%-44s %12s %12s %9s  %s\n", "kernel", "t_min(" A ")", "t_min(" B ")", "A/B", "hash"
      bad = 0
      for (i = 1; i <= n; i++) {
        k = order[i]
        if (!(k in inb)) { printf "%-44s %12.4e %12s %9s  %s\n", k, ta[k], "-", "-", "MISSING in B"; bad = 1; continue }
        st = (ha[k] == hb[k]) ? "same" : "DIFFERS (" ha[k] " / " hb[k] ")"
        if (ha[k] != hb[k]) bad = 1
        printf "%-44s %12.4e %12.4e %9.3f  %s\n", k, ta[k], tb[k], (tb[k] > 0 ? ta[k] / tb[k] : 0), st
        seen[k] = 1
      }
      for (k in inb) if (!(k in seen)) { printf "%-44s %12s %12.4e %9s  %s\n", k, "-", tb[k], "-", "MISSING in A"; bad = 1 }
      exit bad
    }' "$a" "$b"
  rc=$?
  [ $rc -eq 0 ] && echo "every hash identical" || echo "hash differences or missing kernels"
  exit $rc
fi

[ $# -ge 1 ] && [ "${1#-}" = "$1" ] || { echo "usage: $0 NAME [reps=N] | --compare A B"; exit 2; }
name=$1; shift
out=$dir/$name.txt
if [ -e "$out" ]; then echo "refused: $out exists (baselines are never overwritten; pick another name)"; exit 1; fi
exe=exe-validation/fossil_val_perf
[ -x "$exe" ] || { echo "missing $exe (fobis build --mode validation-gnu)"; exit 1; }
if [ -n "${FOSSIL_PERF_COMMIT:-}" ]; then commit=$FOSSIL_PERF_COMMIT; else # override: a tree outside git
  commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)
  if [ -n "$(git status --porcelain --untracked-files=no 2>/dev/null)" ]; then commit="$commit+dirty"; fi
fi
mkdir -p "$dir"
tmp=$(mktemp)
OMP_NUM_THREADS=1 "$exe" "commit=$commit" "host=$(hostname)" "$@" | tee "$tmp"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then echo "the bench failed (exit $rc): no baseline written"; rm -f "$tmp"; exit 1; fi
mv "$tmp" "$out"
echo "baseline written: $out"
