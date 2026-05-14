#!/usr/bin/env bash
# Render a VTU surface as a colour-mapped PNG for use in documentation.
#
# Unlike scripts/render_doc_images.sh (which renders an STL as a plain
# shaded surface), this script colours the surface by a named
# point-scalar field carried in the VTU — the natural visualisation for
# discrete-differential-geometry quantities (curvature, Laplacian
# action) that are scalar fields ON a mesh rather than shapes.
#
# Usage:
#   scripts/render_scalar_field.sh <input.vtu> <scalar_name> <output.png> [view] [pct]
#     view = "iso" (default), "front", "top", "side"
#     pct  = percentile clamp for the colour range, default 0
#            0       -> raw [min, max] of the field
#            P (>0)  -> clamp the LUT to [P-th percentile, (100-P)-th
#                       percentile]. Values outside saturate to the end
#                       colours. Use this when a handful of outlier
#                       vertices (e.g. tessellation-degenerate vertices
#                       on a coarse marching-cubes mesh, where discrete
#                       curvature spikes) would otherwise crush the
#                       colour scale and hide the bulk signal.
#
# Backend: pvbatch (ParaView) only. A scalar-field heatmap has no
# meaningful gnuplot fallback, so if pvbatch is unavailable the script
# warns and exits 1 (CI-tolerant callers should treat this as "image
# not yet rendered" and still produce useful HTML, exactly as
# render_doc_images.sh documents).
#
# Output: 1024x768 PNG, white background, the named scalar mapped
# through ParaView's default Cool-to-Warm diverging colour map with a
# horizontal colour-bar legend. Diverging map is deliberate: curvature
# and Laplacian fields are signed, and a diverging map puts the zero
# level at the neutral midpoint so the sign pattern reads at a glance.
# When a percentile clamp is active the colour-bar shows the clamped
# range, not the raw data range — the figure caption should say so.

set -euo pipefail

INPUT="${1:-}"
SCALAR="${2:-}"
OUTPUT="${3:-}"
VIEW="${4:-iso}"
PCT="${5:-0}"

if [[ -z "$INPUT" || -z "$SCALAR" || -z "$OUTPUT" ]]; then
  echo "usage: $0 <input.vtu> <scalar_name> <output.png> [view] [pct]" >&2
  echo "       view ∈ {iso, front, top, side}, default iso" >&2
  echo "       pct  = percentile clamp (0 = raw min/max), default 0" >&2
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

if [[ -z "$PVBATCH" ]]; then
  echo "warning: no renderer available (need pvbatch); skipping" >&2
  echo "         pages referencing $OUTPUT will render with a missing-image marker" >&2
  exit 1
fi

PYSCRIPT="$(mktemp --suffix=.py)"
trap 'rm -f "$PYSCRIPT"' EXIT
cat > "$PYSCRIPT" <<PYEOF
from paraview.simple import (XMLUnstructuredGridReader, Show, Render, ResetCamera,
                             GetActiveCamera, SaveScreenshot, GetActiveViewOrCreate,
                             ColorBy, GetColorTransferFunction, GetScalarBar,
                             servermanager)

reader = XMLUnstructuredGridReader(FileName=['$INPUT'])
reader.UpdatePipeline()

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1024, 768]
view.Background = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 0
view.UseColorPaletteForBackground = 0

display = Show(reader, view)
display.Representation = 'Surface'

# Colour by the named point-scalar field.
ColorBy(display, ('POINTS', '$SCALAR'))
display.SetScalarBarVisibility(view, True)

lut = GetColorTransferFunction('$SCALAR')

pct = float('$PCT')
if pct > 0.0:
    # Percentile clamp: pull the raw scalar values out of the dataset,
    # sort, and clamp the LUT to [P, 100-P]. Values outside saturate to
    # the end colours. Needed when a few tessellation-degenerate
    # vertices carry extreme outliers that would otherwise crush the
    # colour scale and hide the bulk signal.
    data = servermanager.Fetch(reader)
    arr = data.GetPointData().GetArray('$SCALAR')
    vals = sorted(arr.GetValue(i) for i in range(arr.GetNumberOfTuples()))
    n = len(vals)
    lo_i = max(0, min(n - 1, int(round((pct / 100.0) * (n - 1)))))
    hi_i = max(0, min(n - 1, int(round((1.0 - pct / 100.0) * (n - 1)))))
    lo, hi = vals[lo_i], vals[hi_i]
    if hi <= lo:                       # degenerate field — fall back to raw range
        lo, hi = vals[0], vals[-1]
    lut.RescaleTransferFunction(lo, hi)
    print('percentile clamp P%g: [%g, %g] (raw [%g, %g])'
          % (pct, lo, hi, vals[0], vals[-1]))
else:
    display.RescaleTransferFunctionToDataRange(True, False)

# Cool-to-Warm diverging map — the standard signed-field choice; puts
# the data midpoint at the neutral colour so the sign pattern reads.
lut.ApplyPreset('Cool to Warm', True)

# Horizontal colour bar, bottom-centre, compact.
bar = GetScalarBar(lut, view)
bar.Title = '$SCALAR'
bar.ComponentTitle = ''
bar.Orientation = 'Horizontal'
bar.WindowLocation = 'Lower Center'
bar.ScalarBarLength = 0.4

ResetCamera()
cam = GetActiveCamera()

view_name = '$VIEW'
if view_name == 'iso':
    cam.Azimuth(30)
    cam.Elevation(20)
elif view_name == 'front':
    pass
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
  echo "rendered (pvbatch, scalar='$SCALAR'): $OUTPUT"
  exit 0
fi

echo "warning: pvbatch failed to render $INPUT (scalar '$SCALAR'); skipping" >&2
exit 1
