#!/usr/bin/env bash
# Run the validation suite (issue #45): the programs of src/tests/validation/, built by
#   fobis build --mode validation-gnu        (-O2, into exe-validation/)
# The per-push regression suite is scripts/run_tests.sh (exe/); these tests are the heavy ones (benches, sweeps,
# convergence studies, large meshes) and run nightly and before a release (.github/workflows/validation.yml).
# A test passes on exit status 0 (fossil_xfail_* must fail). Prints one line per test with its time, in name order;
# a failing test also prints the last lines of its log (all of it in exe-validation/NAME.log).
#
# Usage: scripts/run_validation.sh [-j N] [test-name ...]   (default: every program of src/tests/validation/)
#   -j N   run N programs at a time (default 1). The programs are single-threaded and write no files, so they can
#          share the machine; their own timing lines then measure a loaded machine (fossil_val_perf baselines come
#          from scripts/run_perf.sh, one program at a time). The programs of SOLO run first, one at a time with
#          nothing else: their peak resident memory (GNU time, 2026-11-12) is fossil_test_distance_bench 11.9 GB,
#          fossil_val_cut_solver 5.0 GB, every other program under 1.6 GB; a GitHub runner has 16 GB.
set -uo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"
jobs=1
if [ "${1:-}" = "-j" ]; then
  jobs=${2:-}; shift 2 || { echo "usage: $0 [-j N] [test-name ...]"; exit 2; }
  [[ $jobs =~ ^[1-9][0-9]*$ ]] || { echo "usage: $0 [-j N] [test-name ...]"; exit 2; }
fi
if [ $# -gt 0 ]; then names=("$@"); else
  names=(); for src in src/tests/validation/*.f90; do
    grep -qiE '^ *program ' "$src" && names+=("$(sed -n 's/^ *[Pp][Rr][Oo][Gg][Rr][Aa][Mm] \+\([A-Za-z0-9_]*\).*/\1/p' "$src" | head -n 1)")
  done
fi
# fossil_val_perf (issue #94): one timed repetition here (baselines: scripts/run_perf.sh NAME, 5 repetitions).
export FOSSIL_PERF_REPS=${FOSSIL_PERF_REPS:-1}

run_one() {
  # Runs test $1 and writes "exit-status seconds" to exe-validation/$1.result.
  local name=$1 exe=exe-validation/$1 start end rc
  start=$(date +%s.%N); "$exe" > "exe-validation/$name.log" 2>&1; rc=$?; end=$(date +%s.%N)
  echo "$rc $(awk -v a="$start" -v b="$end" 'BEGIN {printf "%.1f", b - a}')" > "exe-validation/$name.result"
}

SOLO=(fossil_test_distance_bench fossil_val_cut_solver)
is_solo() { local s; for s in "${SOLO[@]}"; do [ "$1" = "$s" ] && return 0; done; return 1; }

for name in "${names[@]}"; do rm -f "exe-validation/$name.result"; done
for name in "${names[@]}"; do  # the memory-heavy programs first, alone
  [ -x "exe-validation/$name" ] && is_solo "$name" && run_one "$name"
done
running=0
for name in "${names[@]}"; do
  [ -x "exe-validation/$name" ] || continue
  is_solo "$name" && continue
  if [ "$jobs" -eq 1 ]; then run_one "$name"; continue; fi
  if [ "$running" -ge "$jobs" ]; then wait -n; running=$((running - 1)); fi
  run_one "$name" & running=$((running + 1))
done
wait

pass=0 ; fail=0
for name in "${names[@]}"; do
  if [ ! -x "exe-validation/$name" ]; then echo "  MISSING $name (fobis build --mode validation-gnu)"; fail=$((fail + 1)); continue; fi
  rc=255 ; t='?'
  [ -s "exe-validation/$name.result" ] && read -r rc t < "exe-validation/$name.result"
  if [[ $name == *_xfail_* ]]; then ok=$([ "$rc" -ne 0 ] && echo 1 || echo 0); else ok=$([ "$rc" -eq 0 ] && echo 1 || echo 0); fi
  if [ "$ok" = 1 ]; then echo "  PASS  $name (${t} s)"; pass=$((pass + 1));
  else
    echo "  FAIL  $name (${t} s, exit $rc; log exe-validation/$name.log), its last lines:"
    tail -n 40 "exe-validation/$name.log" | sed 's/^/    | /'
    fail=$((fail + 1))
  fi
done
echo ""
echo "validation: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
