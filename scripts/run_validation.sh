#!/usr/bin/env bash
# Run the validation suite (issue #45): the programs of src/tests/validation/, built by
#   fobis build --mode validation-gnu        (-O2, into exe-validation/)
# The per-push regression suite is scripts/run_tests.sh (exe/); these tests are the heavy ones (benches, sweeps,
# convergence studies, large meshes) and run nightly and before a release (.github/workflows/validation.yml).
# A test passes on exit status 0 (fossil_xfail_* must fail). Prints one line per test with its time.
#
# Usage: scripts/run_validation.sh [test-name ...]   (default: every program of src/tests/validation/)
set -uo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"
if [ $# -gt 0 ]; then names=("$@"); else
  names=(); for src in src/tests/validation/*.f90; do
    grep -qiE '^ *program ' "$src" && names+=("$(sed -n 's/^ *[Pp][Rr][Oo][Gg][Rr][Aa][Mm] \+\([A-Za-z0-9_]*\).*/\1/p' "$src" | head -n 1)")
  done
fi
# fossil_val_perf (issue #94): one timed repetition here (baselines: scripts/run_perf.sh NAME, 5 repetitions).
export FOSSIL_PERF_REPS=${FOSSIL_PERF_REPS:-1}
pass=0 ; fail=0
for name in "${names[@]}"; do
  exe=exe-validation/$name
  if [ ! -x "$exe" ]; then echo "  MISSING $name (fobis build --mode validation-gnu)"; fail=$((fail + 1)); continue; fi
  start=$(date +%s.%N); "$exe" > "exe-validation/$name.log" 2>&1; rc=$?; end=$(date +%s.%N)
  t=$(awk -v a="$start" -v b="$end" 'BEGIN {printf "%.1f", b - a}')
  if [[ $name == *_xfail_* ]]; then ok=$([ $rc -ne 0 ] && echo 1 || echo 0); else ok=$([ $rc -eq 0 ] && echo 1 || echo 0); fi
  if [ "$ok" = 1 ]; then echo "  PASS  $name (${t} s)"; pass=$((pass + 1));
  else echo "  FAIL  $name (${t} s, exit $rc; log exe-validation/$name.log)"; fail=$((fail + 1)); fi
done
echo ""
echo "validation: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
