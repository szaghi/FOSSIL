#!/usr/bin/env bash
# Render an STL file as a PNG preview for use in documentation.
#
# Usage: scripts/render_doc_images.sh <input.stl> <output.png> [view]
#   view = "iso" (default), "front", "top", "side"
#
# Backends, in priority order:
#   1. pvbatch (ParaView)  — full shaded render with antialiasing
#   2. gnuplot 5+          — wireframe-only fallback (lower quality)
#   3. neither             — skip with a warning (CI tolerant: exit 0 in
#                            the "doc images" CI step, exit 1 here so the
#                            caller knows nothing was produced)
#
# Output: 1024x768 PNG at the requested view, neutral lighting,
# slightly off-axis isometric for "iso" so all three faces of an axis-
# aligned object are visible.
#
# Design choice: this script is intentionally idempotent and
# CI-tolerant. Doc builds that depend on these images should treat
# missing PNGs as "image not yet rendered" and still produce useful HTML.

set -euo pipefail

INPUT="${1:-}"
OUTPUT="${2:-}"
VIEW="${3:-iso}"

if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
  echo "usage: $0 <input.stl> <output.png> [view]" >&2
  echo "       view ∈ {iso, front, top, side}, default iso" >&2
  exit 2
fi

if [[ ! -f "$INPUT" ]]; then
  echo "error: input file not found: $INPUT" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUTPUT")"

# Locate pvbatch — prefer the system PATH, then fall back to the common
# /opt install location used on many HPC workstations.
PVBATCH=""
if command -v pvbatch >/dev/null 2>&1; then
  PVBATCH="$(command -v pvbatch)"
else
  for cand in /opt/paraview/*/bin/pvbatch /usr/local/paraview/bin/pvbatch; do
    if [[ -x "$cand" ]]; then
      PVBATCH="$cand"
      break
    fi
  done
fi

GNUPLOT=""
if command -v gnuplot >/dev/null 2>&1; then
  GNUPLOT="$(command -v gnuplot)"
fi

# ---- Backend: pvbatch ---------------------------------------------------------
if [[ -n "$PVBATCH" ]]; then
  PYSCRIPT="$(mktemp --suffix=.py)"
  trap 'rm -f "$PYSCRIPT"' EXIT
  cat > "$PYSCRIPT" <<PYEOF
from paraview.simple import (STLReader, Show, Render, ResetCamera, GetActiveCamera,
                             SaveScreenshot, GetActiveViewOrCreate, ColorBy)

reader = STLReader(FileNames=['$INPUT'])
view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1024, 768]
view.Background = [1.0, 1.0, 1.0]   # white background for docs
view.OrientationAxesVisibility = 0  # hide the X/Y/Z gizmo
view.UseColorPaletteForBackground = 0
# Two-sided lighting plus a brighter fill so unlit faces still read in the
# model's color rather than near-black ParaView default.
view.UseLight = 1
view.KeyLightWarmth = 0.5
view.KeyLightIntensity = 0.85
view.FillLightWarmth = 0.5
view.FillLightKFRatio = 1.5
view.BackLightWarmth = 0.5
view.BackLightKBRatio = 1.5
view.HeadLightWarmth = 0.5
view.HeadLightKHRatio = 1.5

display = Show(reader, view)
display.Representation = 'Surface With Edges'
# Note on color: ParaView's Show() auto-enables ColorBy(Normals) on STL
# files, which paints the surface a deep navy. We accept this default —
# it has good contrast against the white doc background and reads well
# at thumbnail scale. Overriding it cleanly through the PV6 colour
# pipeline (ColorBy(None) + DiffuseColor) ran into multiple
# version-specific quirks; the auto-color is fit for purpose for
# reference imagery.
display.EdgeColor      = [0.05, 0.10, 0.20]
display.LineWidth      = 0.4
# Render back faces too, so thin shells don't disappear at grazing angles.
display.BackfaceRepresentation = 'Follow Frontface'

ResetCamera()
cam = GetActiveCamera()

view_name = '$VIEW'
if view_name == 'iso':
    cam.Azimuth(30)
    cam.Elevation(20)
elif view_name == 'front':
    pass  # default
elif view_name == 'top':
    cam.Elevation(89)
elif view_name == 'side':
    cam.Azimuth(90)

ResetCamera()
Render()
SaveScreenshot('$OUTPUT', view, ImageResolution=[1024, 768],
               TransparentBackground=False)
PYEOF
  if "$PVBATCH" --force-offscreen-rendering "$PYSCRIPT" >/dev/null 2>&1; then
    echo "rendered (pvbatch): $OUTPUT"
    exit 0
  fi
  echo "warning: pvbatch failed; falling back to gnuplot" >&2
fi

# ---- Backend: gnuplot ---------------------------------------------------------
if [[ -n "$GNUPLOT" ]]; then
  # gnuplot can read STL via 'splot' if we convert to (x, y, z) line data.
  # Cheap conversion: emit the 4 vertices per facet (triangle + closing
  # edge to first vertex) with blank lines as record separators.
  TMP_DAT="$(mktemp --suffix=.dat)"
  trap 'rm -f "$TMP_DAT"' EXIT
  awk '
    /^[[:space:]]*facet normal/ { vc = 0; next }
    /^[[:space:]]*vertex/ {
      vc++
      printf("%s %s %s\n", $2, $3, $4)
      if (vc == 1) { v1 = $2 " " $3 " " $4 }
      if (vc == 3) { print v1; print "" }
    }
  ' "$INPUT" > "$TMP_DAT"
  if [[ ! -s "$TMP_DAT" ]]; then
    echo "warning: gnuplot fallback could not parse $INPUT (binary STL?); skipping" >&2
    exit 1
  fi

  case "$VIEW" in
    iso)   GP_VIEW="set view 60, 30, 1, 1" ;;
    front) GP_VIEW="set view 90, 0, 1, 1" ;;
    top)   GP_VIEW="set view 0, 0, 1, 1" ;;
    side)  GP_VIEW="set view 90, 90, 1, 1" ;;
    *)     GP_VIEW="set view 60, 30, 1, 1" ;;
  esac

  "$GNUPLOT" <<GPEOF
set terminal pngcairo size 1024,768 background 'white'
set output '$OUTPUT'
set hidden3d
set xyplane at 0
unset border
unset xtics
unset ytics
unset ztics
unset key
$GP_VIEW
splot '$TMP_DAT' with lines lc rgb '#1c2c4c' lw 0.5
GPEOF
  echo "rendered (gnuplot wireframe): $OUTPUT"
  exit 0
fi

# ---- No renderer available ----------------------------------------------------
echo "warning: no renderer available (need pvbatch or gnuplot); skipping" >&2
echo "         pages referencing $OUTPUT will render with a missing-image marker" >&2
exit 1
