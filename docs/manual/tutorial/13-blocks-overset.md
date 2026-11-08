# 13. Structured blocks and overset grids

The chapters so far put a body into a Cartesian grid. A body-fitted solver puts grids around the body instead:
curvilinear structured blocks, overlapping each other and a background. This chapter builds an O-grid around an
elliptic cylinder, checks it, locates points in it and interpolates from it, then assembles it with a Cartesian
background: holes, receptors, donors, fringes, orphans and the connectivity file a solver reads. Everything is 2-D;
the 3-D types and routines have the same names without `_2d`.

## An O-grid

<<< @/examples/snippets/tut13_blocks_overset-block.f90

- The O-grid has 49 × 17 vertices on the elliptic coordinates of an ellipse of semi-axes 1 and 0.4: i runs around
  the body, j out of it, with the layers clustered at the wall. The vertices are plain arrays `u(ni, nj)`,
  `v(ni, nj)` (3-D: `x, y, z(ni, nj, nk)`); cell (i, j) is spanned by vertices i..i+1, j..j+1.
- `initialize(u, v, name=)` copies them. The last column repeats the first exactly: the seam.
- `detect_topology_2d(blocks)` proposes the topology records of a set of blocks from coincident vertices (exact by
  default, `tolerance=` for grids written with round-off). Here it finds the seam, one `TOPO_ABUT` record on face 1
  (−i) and one on face 2 (+i), each pointing at the other, with the transform (1, 2): direction i of this side runs
  along direction i of the other. Walks and stencils cross the seam through these records. Periodic records come
  from `translations=`, polar axes as `TOPO_SINGULAR`.
- `add_record` appends records by hand. The physical boundaries are `TOPO_BOUNDARY` records with a name: the wall on
  face 3 (−j) and, on face 4 (+j), the boundary that receives from the background in the assembly; the assembly
  reads the names `wall` and `overlap`.

<<< @/examples/snippets/tut13_blocks_overset-validity.f90

- `geometry(area=, centroid=, face_i=, face_j=)` gives the exact cell geometry. The cells sum to −5.992: the sign is
  the block's handedness (i runs counter-clockwise, j outward, so the block is left-handed), not a defect. Its size
  is that of the polygon, 0.3 % below the area between the two ellipses, 6.010.
- `validity(status)` decides each cell exactly: `CELL_VALID`, `CELL_INVALID`, `CELL_SINGULAR` (a collapsed edge, a
  pole). In 2-D the edges are straight, so the four corner orientations (`orient2d`) decide; in 3-D the trilinear
  Jacobian is certified by its Bernstein coefficients.
- `quality(skewness=, orthogonality=, aspect=, min_jacobian=)` gives per-cell metrics: equiangle skewness,
  orthogonality of the faces to the line between neighbouring centroids, aspect ratio, and the certified minimum
  Jacobian over the mean.
- The folded copy drags vertex (13, 5) two layers out, past the vertices of the layer above it: the two cells above
  it, (12, 5) and (13, 5), turn inside out and `validity` flags them; the two below are stretched but stay valid.

## Locating points and interpolating

<<< @/examples/snippets/tut13_blocks_overset-locate.f90

- `build_search` stores the cell validity and the handedness and builds a bounding-box tree over the cells.
- `locate_many(pu, pv, cell, xi, status)` returns per point the cell (i, j), the local coordinates ξ in [0, 1]² and a
  status: `LOC_INSIDE`, `LOC_OUTSIDE`, or, when only a cell that is not valid contains the point, `LOC_WARPED` or
  `LOC_DEGENERATE`. Containment is exact in 2-D; a point on a face shared by two cells belongs to the low one.
- `walk=.true.` jumps from the previous point's cell towards the point (stencil jumping) and hands the point to the
  tree when it cannot settle it, for instance when the walk reaches the wall or the outer boundary: here 38 of the
  200 points fall back, `fallbacks=`. The result equals the robust one bit for bit.

<<< @/examples/snippets/tut13_blocks_overset-stencils.f90

- `stencils(cell, xi, q, nodes, weights, status)` builds a q × q Lagrange stencil on the vertices around each located
  point: `nodes(:, k, n)` are vertex indices, `weights(k, n)` their weights. Near the wall and the outer boundary it is
  shifted one-sided (`STENCIL_SHIFTED`); across the seam it takes the vertices of the other side through the records,
  as the figure shows.
- For q > 2, ξ is recomputed on the isoparametric map of degree q − 1 of the stencil's own vertices
  (`order_matched=`, on by default): the bilinear ξ alone would cap the interpolation at second order.

![the O-grid of 48 x 16 cells around the ellipse, the 200 points of a spiral, teal where locate_many finds them inside the block and ember where they fall inside the body or beyond the outer boundary, and in yellow the 4 x 4 stencil of a point just below the seam (ember line): columns 47, 48 and the seam on its side, column 2 across the seam](../../examples/images/tut13-block.svg){.fossil-figure}

## An overset assembly

The O-grid now overlaps a Cartesian background of 41 × 33 vertices at spacing 0.15. The background covers the body
too: its vertices inside the body must be blanked, those next to the blanked ones and the O-grid's outer boundary
must take their values from the other grid.

<<< @/examples/snippets/tut13_blocks_overset-assembly.f90

- `wall_cutter_2d(blocks, cutter)` builds the cutter, a `curve_object`, from the edges of the faces recorded as
  `wall`, oriented with the solid on the left.
- `hole_flags_2d(blocks, cutter, holes, cuts=)` flags the vertices inside the cutter by its fast winding number:
  exact on a closed wall. `cuts=` says which blocks the cutter cuts: the background, not the O-grid, whose wall lies
  on the cutter. `hole_map_2d` and `xray_flags_2d` are the Cartesian-map and X-ray alternatives.
- `donor_candidates_2d(blocks, donors)` locates every vertex of every block in every other block and keeps all
  candidates, each with its donor cell, ξ, the donor and receptor cell sizes and the certified minimum Jacobian of the
  donor.
- `classify_2d(blocks, donors, masks, SELECT_EXPLICIT, q, holes=)` gives the mask per vertex: `iota` = 1 for a field
  point, 0 for a hole, −k for a receptor from block k. With the explicit rule the ring around the holes and the
  vertices of the `overlap` faces must receive; a candidate is valid when its q × q donor stencil holds no hole and no
  receptor. Validity depends on the neighbours, so the classification is iterated to a fixed point; `circular=` counts
  receptor pairs in each other's stencils. `SELECT_IMPLICIT` chooses by donor cell size instead.
- `fringe_levels_2d` sets `masks(b)%fringe`: INNER 1 and 2 around the holes, OUTER 1 and 2 at the overlap boundary.
  `overlap_width_2d` is Chesshire's width, the least distance between the receptors of the two blocks.
- `build_connectivity(masks, donors, q, 2, conn)` makes the receptor list and `save_connectivity(conn, file)` writes
  it as text: per block the counts by class, then one line per receptor with its block and vertex, the donor block,
  the donor cell, ξ, the stencil width, a low-order flag and the donor-to-receptor size ratio.

The overlap has to be wide enough for the stencils. The same assembly with the O-grid cut to fewer layers:

<<< @/examples/snippets/tut13_blocks_overset-orphans.f90

- A receptor without a valid donor is an orphan; `masks(b)%orphan` gives the cause: `ORPHAN_STENCIL_HOLE` (every
  candidate stencil touches a hole), `ORPHAN_RECEPTOR_DONOR` (it holds a receptor, explicit interpolation),
  `ORPHAN_NO_CANDIDATE` (no donor cell at all), `ORPHAN_BOUNDARY`, `ORPHAN_WALL`.

![the O-grid cut to 11 cell layers over the Cartesian background: the background holes inside the body in grey, the background receptors around them in teal, the O-grid's outer vertices that receive in yellow and its orphans in ember, at the two ends of the long axis where the O-grid is thinnest](../../examples/images/tut13-overset.svg){.fossil-figure}

## Running it

<<< @/examples/output/tut13.ansi{ansi}

- The seam gives two records; all 768 cells are valid. The skewness is at most 0.103 and the orthogonality at
  least 0.986, since elliptic coordinates are orthogonal; the aspect ratio reaches 8.06 in the flat cells at the wall.
- 165 of the 200 spiral points are inside the block; the walk returns the same cells and ξ bit for bit.
- The interpolation error of sin(2x) cos(y) falls from 9.78e-3 at q = 2 to 1.75e-3 at q = 3 and 8.53e-4 at q = 4.
  Fewer stencils are centred as q grows (165, 149, 133): the others are shifted at the wall or the outer boundary.
- The cutter blanks 57 background vertices; the 40 around them receive from the O-grid and the 49 vertices of the
  O-grid's outer boundary from the background. No orphan, two iterations of the fixed point, no circular pair; the
  overlap width is 0.540. The connectivity file lists the 89 receptors; the first is the outer vertex 785 of the
  O-grid, (1, 17), in background cell (32, 16) at ξ = (0.097, 1.0): on that cell's top edge, which the low-side rule
  gives to the lower cell.
- With 13 layers 8 of the 49 outer vertices are orphans, with 11 layers 31, with 9 all of them. The O-grid is
  thinnest at the ends of the long axis, and the orphans start there (figure): their stencils reach the receptors
  around the holes, or the holes themselves. A solver cannot interpolate there; the grid has to change.

::: tip What you learned
A `curvilinear_block_2d_t` with its topology records (`detect_topology_2d`, `add_record`), exact cell geometry,
validity and quality, located points (`build_search`, `locate_many`, `walk=`) and Lagrange stencils (`stencils`). An
overset assembly from kernels: `wall_cutter_2d`, `hole_flags_2d`, `donor_candidates_2d`, `classify_2d`,
`fringe_levels_2d`, orphans by cause, `build_connectivity` and `save_connectivity`. Reference:
[Curvilinear blocks](/guide/advanced/curvilinear-blocks), [Overset assembly](/guide/advanced/overset).
:::

Next: [14. Contact: distance and collision](./14-contact).
