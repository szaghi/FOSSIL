#!/usr/bin/env bash
# Regenerate every documentation figure under docs/pictures/advanced/.
#
# This is the single entry point that ties together the three pieces of
# the doc-image pipeline:
#
#   1. exe/gen_doc_fixtures   — Fortran program (src/app/gen_doc_fixtures.f90)
#                               that emits the STL + VTU fixtures into
#                               docs/pictures/advanced/_fixtures/.
#   2. render_doc_images.sh   — STL  -> shaded-surface PNG.
#   3. render_scalar_field.sh — VTU  -> colour-mapped-scalar PNG.
#
# Most advanced-feature figures were rendered by hand and committed; this
# script covers the discrete-differential-geometry pages (smoothing,
# curvature, cotangent Laplacian) whose fixtures are reproducible from
# the library itself. Extend the FIGURES / SCALAR_FIGURES tables below as
# new pages need imagery.
#
# Usage:
#   scripts/make_doc_images.sh
#
# Pre-requisite: build the fixture generator first —
#   fobis build --mode tests-gnu
# (gen_doc_fixtures is an app target, built alongside the test suite).
#
# CI-tolerant: a missing renderer (pvbatch / gnuplot) is a warning, not a
# hard failure — the script reports what it could not produce and exits 0
# so a docs build can still proceed with missing-image markers.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

GEN="exe/gen_doc_fixtures"
FIXDIR="docs/pictures/advanced/_fixtures"
OUTDIR="docs/pictures/advanced"

if [[ ! -x "$GEN" ]]; then
  echo "error: $GEN not found — build it first with:" >&2
  echo "       fobis build --mode tests-gnu" >&2
  exit 2
fi

echo "==> generating fixtures"
"$GEN"

# ---- STL figures: <fixture.stl> <output.png> <view> --------------------------
# Smoothing before/after pair.
declare -a STL_FIGURES=(
  "$FIXDIR/smoothing-bunny-noisy.stl  $OUTDIR/smoothing-bunny-noisy.png  iso"
  "$FIXDIR/smoothing-bunny-taubin.stl $OUTDIR/smoothing-bunny-taubin.png iso"
)

# ---- Scalar-field figures: <fixture.vtu> <scalar> <output.png> <view> <pct> --
# The curvature fields carry extreme outliers at the handful of
# tessellation-degenerate vertices on the coarse marching-cubes sphere
# (discrete curvature spikes there — the documented coarse-mesh
# limitation). A 2nd/98th-percentile clamp keeps the bulk signal
# visible; the figure captions note the clamp. The Laplacian field has
# no such outliers, so it renders at the raw min/max range (pct 0).
declare -a SCALAR_FIGURES=(
  "$FIXDIR/curvature-sphere.vtu K  $OUTDIR/curvature-sphere-gaussian.png iso 2"
  "$FIXDIR/curvature-sphere.vtu H  $OUTDIR/curvature-sphere-mean.png     iso 2"
  "$FIXDIR/laplacian-sphere.vtu Lx $OUTDIR/laplacian-sphere-Lx.png       iso 0"
)

n_ok=0
n_fail=0

echo "==> rendering STL figures"
for spec in "${STL_FIGURES[@]}"; do
  # shellcheck disable=SC2086
  set -- $spec
  if scripts/render_doc_images.sh "$1" "$2" "$3"; then
    n_ok=$((n_ok + 1))
  else
    echo "  (skipped: $2)" >&2
    n_fail=$((n_fail + 1))
  fi
done

echo "==> rendering scalar-field figures"
for spec in "${SCALAR_FIGURES[@]}"; do
  # shellcheck disable=SC2086
  set -- $spec
  if scripts/render_scalar_field.sh "$1" "$2" "$3" "$4" "$5"; then
    n_ok=$((n_ok + 1))
  else
    echo "  (skipped: $3)" >&2
    n_fail=$((n_fail + 1))
  fi
done

echo "==> done: $n_ok rendered, $n_fail skipped"
# CI-tolerant: never hard-fail on a missing renderer.
exit 0
