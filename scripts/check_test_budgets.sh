#!/usr/bin/env bash
# check_test_budgets.sh — enforce the time budget of the regression tests (issue #45).
#
# Usage: scripts/check_test_budgets.sh [--factor F]
#
# Runs every executable of exe/ as scripts/run_tests.sh does (expected-failure tests included; their exit status is
# not checked here, run_tests.sh does that) and times it. A test (fossil_test_*, fossil_xfail_*) slower than its
# budget times the factor fails the check; the other executables (the apps, which run_tests.sh runs too) are timed
# and listed, with no budget. Run it on the build the budgets are set for: `tests-gnu-debug --coverage` (the `makecoverage` rule of fobos
# calls it right after run_tests.sh, which is what CI runs).
#
# Budgets: 2 s, or the seconds given for the test in src/tests/budgets.txt (lines "name seconds", # comments).
# Factor: --factor F, else $FOSSIL_BUDGET_FACTOR, else 2 on CI (CI=true: a GitHub runner is slower and noisier than a
# workstation; tighten once its times are known) and 1 elsewhere.
#
# A test over its limit is timed a second time and the faster run counts: a spike of the machine is not a regression.
# Always prints the five slowest tests, so a drift shows before a test breaks its budget.

set -u
cd "$(git rev-parse --show-toplevel 2>/dev/null || dirname "$0"/..)" || exit 2

DEFAULT_BUDGET=2
FACTOR="${FOSSIL_BUDGET_FACTOR:-}"
[[ -z "$FACTOR" && "${CI:-}" == "true" ]] && FACTOR=2
[[ -z "$FACTOR" ]] && FACTOR=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --factor | -f ) FACTOR="$2"; shift 2 ;;
    * ) printf "Unknown argument: %s\n" "$1" >&2; exit 2 ;;
  esac
done
[[ "$FACTOR" =~ ^[0-9]+([.][0-9]+)?$ ]] || { printf "Invalid factor: %s\n" "$FACTOR" >&2; exit 2; }

declare -A BUDGET=()
if [[ -f src/tests/budgets.txt ]]; then
  while read -r name secs _; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    BUDGET[$name]=$secs
  done < src/tests/budgets.txt
fi

run_timed() { # seconds taken by one run of executable $1 (name $2)
  local start end
  start=$(date +%s.%N)
  if [[ "$2" == *mpi* ]]; then mpirun -np 2 "$1" > /dev/null 2>&1 ; else "$1" > /dev/null 2>&1 ; fi
  end=$(date +%s.%N)
  awk -v a="$start" -v b="$end" 'BEGIN { printf "%.2f", b - a }'
}

times=$(mktemp) ; apps=$(mktemp) ; trap 'rm -f "$times" "$apps"' EXIT
over=0
shopt -s nullglob
for exe in exe/*; do
  [[ -f "$exe" && -x "$exe" ]] || continue
  name=$(basename "$exe")
  secs=$(run_timed "$exe" "$name")
  if [[ "$name" != fossil_test_* && "$name" != fossil_xfail_* ]]; then
    printf '%s %s\n' "$secs" "$name" >> "$apps"
    continue
  fi
  budget=${BUDGET[$name]:-$DEFAULT_BUDGET}
  limit=$(awk -v b="$budget" -v f="$FACTOR" 'BEGIN { printf "%.2f", b * f }')
  # one spike of the machine is not a regression: a test over its limit is timed again, the faster run counts
  if awk -v s="$secs" -v l="$limit" 'BEGIN { exit !(s > l) }'; then
    again=$(run_timed "$exe" "$name")
    secs=$(awk -v a="$secs" -v b="$again" 'BEGIN { printf "%.2f", (b < a ? b : a) }')
  fi
  printf '%s %s %s\n' "$secs" "$name" "$limit" >> "$times"
  if awk -v s="$secs" -v l="$limit" 'BEGIN { exit !(s > l) }'; then
    printf "  OVER  %-45s %6s s > %s s (budget %s s x %s)\n" "$name" "$secs" "$limit" "$budget" "$FACTOR"
    over=$((over + 1))
  fi
done

printf "\nslowest tests (factor %s):\n" "$FACTOR"
sort -rn "$times" | head -5 | awk '{ printf "  %6.2f s  %-45s (limit %s s)\n", $1, $2, $3 }'
total=$(awk '{ s += $1 } END { printf "%.1f", s }' "$times")
count=$(wc -l < "$times")
printf "%d tests, %s s in all; %d over budget\n" "$count" "$total" "$over"
if [[ -s "$apps" ]]; then
  printf "other executables (no budget):\n"
  sort -rn "$apps" | awk '{ printf "  %6.2f s  %s\n", $1, $2 }'
fi
exit $((over > 0))
