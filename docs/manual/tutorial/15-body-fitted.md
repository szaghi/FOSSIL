# 15. Body-fitted grids

A body-fitted grid has its first grid line on the wall. FOSSIL does not generate such grids, but it supplies the
geometry their generators and solvers need: points put on the body without crossing its edges, a 2-D hyperbolic
O-grid marched out of a section, and a report of the wall spacing against a y₁⁺ target. The body of this chapter is
the NACA0012 wing `naca0012-ascii.stl` (chord x ∈ [-0.25, 0.75], span z ∈ [0, 1], a blunt trailing edge 0.0025
thick), whose sharp trailing edge is where an unconstrained projection goes wrong.

<<< @/examples/snippets/tut15_body_fitted-target.f90

- `projection_target(target)` builds, once, what every projection reads: the triangles grouped by named patch (one
  tree per patch), the **constraint curves** and the **corners**. The curves are the feature edges of the surface
  (sharp, boundary, non-manifold; `angle=` sets the dihedral threshold, 30 degrees by default) plus the boundaries
  between patches, chained between corners. A corner is a vertex where the constraint edges meet in a number other
  than two, or where a curve turns by more than `corner_angle` (30 degrees).
- The wing has one patch, 12 constraint curves (the root and tip outlines, the two trailing-edge lines, the edges of
  the blunt base) and 7 corners: the four trailing-edge corners and three points of the root and tip outlines near
  the leading edge.

<<< @/examples/snippets/tut15_body_fitted-points.f90

- `classify_points(target, x, y, z, tolerance, kinds, ids)` gives each point its constraint: `PROJ_CORNER` if a
  corner is within the tolerance, else `PROJ_CURVE` if a constraint curve is, else `PROJ_PATCH` with the patch of the
  closest facet. The three points sit 0.005 above the upper trailing-edge line, within the tolerance 0.02 of it: all
  three get `PROJ_CURVE`.
- `project_points(target, x, y, z, kinds, ids, fx, fy, fz)` moves each point to the closest point of its target: the
  whole surface (`PROJ_FREE`), one patch, one constraint curve (or any, id 0), one corner. Ties go to the lowest
  facet or segment index, so the result does not depend on the trees.
- Free, the feet land on the upper surface at x = 0.7494, ahead of the trailing edge; constrained, they land on the
  trailing-edge line itself, (0.75, 0.00126). A grid line meant to run along the edge would zig-zag across it with
  the free feet.

<<< @/examples/snippets/tut15_body_fitted-face.f90

- `project_block_face(block, face, target, tolerance, moved=, kinds=, status=)` classifies and projects every vertex
  of one face of a `curvilinear_block_t` (chapter 13), in place. The faces are numbered 1-6 = -i, +i, -j, +j, -k, +k;
  face 5 is k = 1, the layer 0.005 above the upper side. `moved` is the distance each vertex moved and `kinds` the
  constraint it took.
- Vertices within the tolerance of a wing corner go to it, those within it of a constraint curve to the curve, the
  others to the patch: 7 vertices on corners (the two trailing-edge corners of the face among them), 132 on curves,
  722 on the patch. The trailing-edge row lies on x = 0.75 and the root and tip rows on z = 0 and z = 1, exactly.
  In 2-D, `project_block_side` does the same for one side of a `curvilinear_block_2d_t`.

![the NACA0012 wing near its tip, translucent, seen from behind and above with the span drawn horizontal: the 41 × 21 face grid projected by project_block_face in yellow, its rows following the upper surface, its last chordwise row on the trailing edge and its last spanwise row on the tip outline](../../examples/images/tut15-face.svg){.fossil-figure}

## A hyperbolic O-grid

<<< @/examples/snippets/tut15_body_fitted-ogrid.f90

- `slice_surface` (chapter 7) gives the mid-span section, 60 vertices, counter-clockwise. Splitting each segment in
  two keeps the polygon and doubles the wall points to 120.
- The first height comes from the y₁⁺ target: y₁ = y₁⁺ ν / u_τ. `utau_flat_plate(u_inf, nu, x)` estimates u_τ from a
  turbulent flat plate, C_f = 0.0576 Re_x^(−1/5), u_τ = U sqrt(C_f / 2). At Re = 10⁶ and mid-chord it gives
  u_τ = 4.569e-2, so h₁ = 2.189e-5 for y₁⁺ = 1.
- `geometric_spacing(h1, ratio, layers)` returns the layer heights h₁ qᵏ: 40 layers with ratio 1.2 reach 0.161 from
  the wall.
- `hyperbolic_grid_2d(u, v, spacing, block, invalid=, status=)` marches the layers outward from a closed curve
  (counter-clockwise; `closed=.false.` for an open one). Each step is normal to the current layer, with the mitre
  length h / cos(β/2) at a point where the layer turns by β, so a layer stays at height h from both sides of a
  corner; the step directions are smoothed implicitly (`epsilon`, 0.5 by default). The result is a
  `curvilinear_block_2d_t` of 121 × 41 vertices (the first point repeated at i = 121); every cell is checked by the
  certified Jacobian of the [curvilinear blocks](/guide/advanced/curvilinear-blocks), and `status = 2` with
  `invalid > 0` would report a fold. Here there is none.

![the 121 × 41 O-grid marched from the mid-span section of the wing, every grid line drawn, the wall in yellow: the 40 layers packed against the wall, where h₁ is 2.2e-5, and opening out to 0.16, with fans of lines around the leading edge and the two corners of the blunt trailing edge](../../examples/images/tut15-ogrid.svg){.fossil-figure}

## The wall-spacing report

<<< @/examples/snippets/tut15_body_fitted-report.f90

- The report reads the wall from the block's topology records (chapter 13): every `TOPO_BOUNDARY` record named
  `wall`, here the side j = 1. `wall_cutter_2d(blocks, curve)` turns those sides into a curve, the true wall.
- `block_report_2d(block, cutter, nu, report, ...)` (`block_report` for a 3-D block and a surface cutter) measures
  per wall point:
  - `y1_line`, the first step along the wall-normal grid line, and `y1_wall`, the distance of the first off-wall point
    to the cutter: the true height, smaller when the grid line leaves the wall at an angle;
  - `y1p` = y1_wall u_τ / ν, with u_τ given (`u_tau=` for all points, `u_tau_wall=` per point) or estimated by the
    flat-plate correlation from `u_inf=` and the leading edge `x_le=` (flagged in `utau_estimated`);
  - `growth`, the largest ratio of consecutive steps along the grid line;
  - and over the block the cell quality of chapter 13: skewness, orthogonality, aspect, certified validity.
- `save_block_report(report, file_name)` writes the totals and one line per wall point as plain text.

## Running it

<<< @/examples/output/tut15.ansi{ansi}

- Free projection puts the trailing-edge points on the upper surface, at x = 0.7494; constrained, on the edge.
  The projected face has its trailing-edge, root and tip rows exactly on their lines; the largest move is 0.0085.
- With the u_τ the first height was designed with, y₁⁺ has minimum 0.9999 and mean 1.0065: the O-grid delivers its
  target. The maximum, 1.3303, is at a corner of the blunt trailing edge, (0.75, -0.00126): the mitre step puts the
  first point h₁ from both sides of the corner, so h₁ / cos(β/2) from the corner vertex itself. `above` counts every
  point with y₁⁺ > 1, round-off excess included: 0.752.
- With the flat-plate estimate at each point the picture changes: y₁⁺ is 0 at the leading edge, where x − x_le = 0
  and the estimate gives u_τ = 0, and 2.0163 just behind it, at x = -0.24955, where Re_x is small and the estimate
  large. Such estimates can be about 300 % off near a leading edge: use them to size a grid before a run, then pass
  the solver's own u_τ (`u_tau_wall=`) to check it after.
- The growth ratio stays at 1.2 (1.2031 at most), skewness reaches 0.9029 at the trailing-edge fans, and no cell
  fails the certified validity.
- The file holds the totals (`# ...` lines; the cell-quality line is skipped in the output above) and one line per
  wall point: indices, face, y1_line, y1_wall, y1+, growth.

::: tip What you learned
`projection_target` once, then `classify_points` / `project_points` or `project_block_face` to put grid points on a
body without crossing its edges; `geometric_spacing` from a y₁⁺ target and `hyperbolic_grid_2d` for an O-grid;
`block_report_2d` / `block_report` and `save_block_report` to check the wall spacing and the cells. Reference:
[Body-fitted grid support](/guide/advanced/body-fitted), [Curvilinear blocks](/guide/advanced/curvilinear-blocks).
:::

Next: [16. On the GPU](./16-gpu).
