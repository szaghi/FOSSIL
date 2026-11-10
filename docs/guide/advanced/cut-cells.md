---
title: Cut-cell geometry
---

# Cut-cell geometry

## What it does

FOSSIL computes, on a rectilinear (possibly stretched) block of an AMR grid,
where the grid meets the wall:

- **Grid-line crossings.** Every crossing of an axis-aligned line with the
  surface (curve), sorted, with its facet (segment) and whether the line
  enters or leaves the solid. From them, for each cell centre, the fraction
  of the spacing toward each neighbour at which the wall is met, and the
  inside flag. These serve finite-difference closures: shortened spacing,
  ILW boundary data, fluxes at grid-line intersections.
- **Cut cells, 2-D:**
  - per cell: fluid area fraction, the four edge apertures, fluid centroid,
    and wetted length, centroid and normal;
  - the separate fluid pieces of a cell, each with its own measures.
- **Cut cells, 3-D:**
  - face fluid fractions and centroids from exact plane sections;
  - per cell: fluid volume fraction and centroid, and wetted area, centroid
    and normal;
  - the separate fluid pieces of a cell (`cut%pieces`), each with its own
    volume and centroid, wetted area, centroid and normal, and the fluid
    fraction of each of the six cell faces (`cut%piece_*`).

The fluid is outside the solid. These outputs serve conservative cut-cell
closures and sub-cell masks.

| Query | 3-D | 2-D |
|---|---|---|
| Crossings of one line | `surface%line_crossings` | `curve%line_crossings` |
| Crossings of many parallel lines, fixed buffers | `surface%line_crossings_many(_device)` | `curve%line_crossings_many(_device)` |
| Wall fractions along the grid lines of a block | `surface%cell_wall_fractions` | `curve%cell_wall_fractions` |
| Cut cells of a block | `cut_cells_3d(surface, x, y, z, cut)` | `curve%cut_cells(x, y, cut)` |
| Cut cells of a block on the device, copied back or resident | | `curve%cut_cells_device(x, y, cut, resident=)` |

![a solid ellipse, turned by 0.4 rad, on 8 by 6 cells of 0.125 seen from above: the cells shaded by fluid fraction, cream in the solid to sepia in the fluid; the 16 crossings of the grid lines with the ellipse as gold dots; on every edge the ellipse crosses, the fluid part of the edge in red, ending at the crossing; and from the wetted centroid of every cut cell a teal stroke along its wetted normal integral, pointing out of the ellipse](../../examples/images/ref-cut-cells-outputs.svg){.fossil-figure}

The 2-D outputs on one picture: 16 of the 48 cells are cut and the grid lines cross the ellipse 16 times. The solid
area of the cells is 0.22610, the polygon's area; the wetted lengths add up to 1.79520, and the aperture differences
equal the wetted normal integral in every cell, residual 0.00e+00 (`fig_ref_cut_cells_outputs.f90`).

### On the device

`curve%cut_cells_device(x, y, cut, status, intervals, resident)` (issue #100) returns the `cut` of `curve%cut_cells`
bit for bit:
- fractions, apertures, centroids;
- wetted lengths, centroids and normals;
- pieces and their measures;
- `node_in` and the fluid edge intervals.

Moments, boundaries and quadrature are host only, and the twin has no arguments for them.

How the work is split:
- **Host:** walks the segments through the grid and indexes the crossings.
- **Device:** computes the node flags (the exact tests of the host, one thread per node), classifies the cells, runs
  the integer prefix sums and the same per-cell kernels as `cut_cells` (count, then fill).

The two agree by construction: same kernels, and integer sums are exact in any order. The curve needs no
`enter_device`.

**Resident outputs.** With `resident=.true.` the outputs stay on the device:
- the host arrays of `cut` are allocated but hold garbage until `cut%update_host`;
- `cut%on_device` says so;
- `cut%exit_device` releases them, and must be called before the cut goes out of scope or is deallocated.

`cut_cells`, `cut_cells_device` and `cut_coarsen_2d` release a resident cut they are given before overwriting it. A
routine that reads a cut on the host (`cut_links_2d`, `cut_coarsen_2d`) needs `update_host` first. A copy `b = a` of
a resident cut is not resident. Device code reads the outputs through pointer aliases, the rule of the other device
routines:

```fortran
program resident_cut
use fossil, only : curve_object, cut_cells_2d_t
use penf,   only : I4P, R8P
implicit none
type(curve_object)           :: curve
type(cut_cells_2d_t), target :: cut
real(R8P),    allocatable    :: u(:), v(:), x(:), y(:)
real(R8P),    pointer        :: alpha(:,:)
real(R8P)                    :: fluid
integer(I4P)                 :: i, j, k, status

allocate(u(256), v(256), x(129), y(129))
do k = 1, 256
   u(k) = 0.5_R8P + 0.3_R8P * cos(0.0245436926_R8P * (k - 1)) ; v(k) = 0.5_R8P + 0.3_R8P * sin(0.0245436926_R8P * (k - 1))
enddo
call curve%adopt_loops(u=u, v=v, loop_sizes=[256_I4P])
do k = 1, 129
   x(k) = (k - 1) / 128._R8P ; y(k) = x(k)
enddo
call curve%cut_cells_device(x=x, y=y, cut=cut, status=status, resident=.true.)
alpha => cut%alpha                       ! a device consumer: the fluid area, summed on the device
fluid = 0._R8P
!$acc parallel loop collapse(2) reduction(+:fluid) present(alpha)
do j = 1, 128
   do i = 1, 128
      fluid = fluid + max(alpha(i, j), 0._R8P) / 128._R8P**2
   enddo
enddo
call cut%update_host                     ! only if the host needs the outputs
call cut%exit_device                     ! before cut goes away
print '(A,F8.5)', 'fluid area: ', fluid
endprogram resident_cut
```

Measured on an RTX 4070 (nvfortran 26.9, `-gpu=nofma`), best of 5 calls with the fluid intervals:

| curve | cells | CPU [ms] | device, copied back [ms] | device, resident [ms] |
|---|---|---|---|---|
| ellipse, 1024 segments | 1024² | 41 | 27 | 6.4 |
| NACA 0012, 1022 segments | 1024² | 42 | 30 | 6.6 |
| 50 circles, 3200 segments | 1024² | 68 | 32 | 9.3 |
| any of the three | 256² | 2.8-5.8 | 3.4-4.4 | 1.2-2.0 |

What bounds each path at 1024²:
- **Copied back:** the copy of about 250 MB of outputs to pageable host memory (12-14 GB/s into pages already written,
  3-4 GB/s into new ones: the twin reuses the host arrays of the cut it overwrites when the shape is unchanged).
- **Resident:** about 3 ms on the device (kernels 1.3 ms, 13 MB of indexes copied in) and about 3.6 ms on the host
  (walk, indexes, allocations).

**Per time step.** `fossil_test_curve_bench 256` also times a moving body: each step moves the curve
(`move`, 2 degrees about the block centre), tags the cells (`box_classify_many`, one box per cell) and computes the
cut cells, 20 steps, per-step best of 5:

| curve | cells | CPU [ms] | device, copied back [ms] | device, resident [ms] |
|---|---|---|---|---|
| ellipse | 1024² | 138 | 60 | 32 |
| NACA 0012 | 1024² | 111 | 63 | 31 |
| ellipse | 256² | 9.7 | 6.7 | 4.5 |
| NACA 0012 | 256² | 7.6 | 7.0 | 4.5 |

The resident step is bound by the tags, about 25 of its 32 ms at 1024², because `box_classify_many_device` uploads the
one million cell boxes and returns the codes at every step. The swept geometry and the piece links are not on the
device yet (issue #100, steps 4 and 5), so a whole solver step still goes back to the host.

## Example

The basic calls, on the dragon section and the dragon, are in
[tutorial chapter 10](/manual/tutorial/10-cut-cells), and the solver data built on them (links, SRD and merge
neighbourhoods) in the [cookbook](/manual/cookbook#neighbourhoods-for-a-small-cell-cure). This program runs the
[checks](#checks-you-can-make) of this page on small cases with known answers.

<<< @/examples/snippets/ref_cut_cells.f90

<<< @/examples/output/ref_cut_cells.ansi{ansi}

- **Thin strip, 2-D.** The solid strip, 0.05 thick, lies inside row 2 (v from 0.25 to 0.5), so each cell of the row
  holds two fluid pieces, of 0.4 of the cell each, and `alpha` is 0.8. `piece_aperture` splits the edges between
  them: each piece owns 0.4 of the -u and +u edges, the lower piece all of the -v edge and none of the +v edge, the
  upper piece the reverse.
- **Turned cube, 3-D.** `move` turns the unit cube 0.3 rad about (1, 1, 1) and shifts it, so no face is
  grid-aligned. On the 6³ block, 125 cells are cut and none is flagged. 7 cells hold two pieces or more: the body is
  convex, but where an edge of the cube runs through a cell, the solid wedge can split the cell's fluid in two.
- **Covering sum.** The solid volume summed over the block differs from `get_volume` by 5.55·10⁻¹⁶.
- **Divergence theorem.** The largest residual of the face-fraction differences against `wet_normal`, over every
  cell, is 3.22·10⁻¹⁶.
- **Small cells.** 5 cells have a fluid fraction below 10⁻², the smallest 1.576·10⁻⁵: the cells a
  solver's small-cell treatment has to handle.
- **Wall fractions.** `cell_wall_fractions` works on the cell centres, not the edges. Centre (2, 3, 3) is inside,
  and its -x line meets the wall at 0.270 of the spacing to the neighbouring centre; -1 marks the directions with no
  wall before the next centre.

The batched line query, on the dragon: 24 × 24 lines along x through its bounding box, into buffers of 16 crossings
per line, each line checked against the single-line call.

<<< @/examples/snippets/ref_cut_cells_lines.f90

<<< @/examples/output/ref_cut_cells_lines.ansi{ansi}

- **Buffers.** 337 of the 576 lines cross the dragon, 1460 crossings in all, the most crossed line 14 times: a buffer
  of 16 does not overflow (with 8, 10 lines overflowed). `overflow(i)` flags a line whose buffer was too small; its
  kept crossings are then not the first ones, so re-run it with a larger buffer.
- **Agreement.** Every line gives the same parameters and facets as `line_crossings`, and enters as often as it
  leaves: the dragon is not watertight (`is_watertight` is false), but none of these lines passes through an opening.

## The exactness model

- **Decisions are exact, constructed coordinates are rounded.**
  - Every topological decision is a comparison or the exact `orient2d`:
    - whether a line crosses a triangle or segment;
    - on which side of a grid plane or line a vertex lies;
    - in which cell a piece of a segment lies.
  - Crossing points are constructed in floating point.
  - So the topology is exact and the measures are exact to round-off.
- **One symbolic tie-break.** A point exactly on a grid line or plane counts
  as being on its low side (the line or plane is taken at x + ε). An
  axis-aligned line through an edge or a vertex of the surface is resolved
  by a symbolic perturbation of its position (Woop et al. 2013, made exact).
  - On a watertight surface every crossing is counted exactly once, also on
    grid-aligned bodies.
  - Every point lies in exactly one half-open cell, (x_i, x_i+1].

![a diamond on an 8 x 8 grid with its four vertices on grid nodes and its sides through grid nodes, the cells shaded by fluid fraction (cream solid, sepia fluid, the cut cells half and half), and the crossings of every horizontal grid line marked at the nodes, entering in teal on the left and leaving in ember on the right; the line through the bottom vertex has two crossings, the line through the top vertex none](../../examples/images/ref-cut-cells-crossings.svg){.fossil-figure}

Every crossing of the diamond lies on a grid node, the worst case for a rounded test. The line v = 1/8 touches
the bottom vertex: taken at v + ε it passes just above it and crosses both sides there, so it has 2 crossings, both at
u = 0.5. The line v = 7/8 touches the top vertex and, taken at v + ε, misses it: 0 crossings. The other lines have one
entering and one leaving crossing each, and the 12 cut cells are halved exactly: the solid area is 0.28125, the area
of the diamond.
- **2-D cut cells.**
  - Each segment is walked through the grid. Whether it leaves a cell
    through the vertical or the horizontal line first is decided by the
    exact side of the cell corner.
  - The measures come from Green's theorem over the wall pieces and the
    solid parts of the cell edges.
  - The pieces come from a union-find over the fluid boundary arcs, joined
    by the wall chains that connect them. A clockwise loop inside a cell (a
    hole) is a piece of its own; a counter-clockwise loop (an island) is
    subtracted from the piece that contains it.
- **3-D cut cells.**
  - Each grid plane is cut once (`slice_surface` with `on_plane_low`), and
    the section is cut by the 2-D machinery, so the two cells sharing a face
    see the same fraction.
  - Every triangle is clipped to every cell it overlaps, with the same rule.
  - Volume and first moments come from the divergence theorem over the
    clipped triangles and the solid face parts, in coordinates local to each
    cell.
  - **Pieces** come from a union-find over the parts of the boundary of the
    fluid:
    - the parts are the wall polygons (each triangle clipped to the cell) and
      the fluid intervals of the 12 cell edges;
    - the corner states are the exact inside flags of the perturbed nodes,
      and the piercings of an edge are its exact crossings under the global
      perturbation (between the perturbed corners), entering or leaving the
      solid by the sign of the normal along the edge, so the states
      alternate from corner to corner (issue #70; before, the piercings were
      the clipped polygon vertices on the edge, one per incident facet when a
      mesh edge or vertex lay on it);
    - the parts are joined where polygons share a vertex (the clipped points
      of a shared mesh edge are bit-identical), at a piercing, at a fluid
      corner, and where a solid island traced on a face is surrounded by
      fluid (a ray from the island to the nearest trace or face edge);
    - every joined set that touches the cell boundary is one piece; a closed
      wall set inside the cell adds a piece only if it encloses a void.
  - So a passage narrower than the cell, whose fluid meets the cell boundary
    in two separate disks, is one piece, and a thin plate across the cell
    leaves two.

![a plate 0.075 thick, 0.3 of the cell side, tilted across a 4 x 4 x 4 block and seen almost edge-on; the 22 cells it splits into two fluid pieces outlined in yellow, the centroid of the first piece of each in teal on one side of the plate and of the second in ember on the other](../../examples/images/ref-cut-cells-pieces.svg){.fossil-figure}

The plate cuts 32 cells and splits 22 of them in two; none is flagged. The two centroids of each split cell lie on
opposite sides of the plate.

  - **Piece measures.** The pieces of cell (i, j, k) are
    `first_piece(i, j, k)` to `first_piece(i, j, k) + pieces(i, j, k) - 1`;
    `piece_cell(:, p)` gives the cell back. For each piece:
    - the wall terms (area, centroid, normal, divergence-theorem volume and
      moments) of the wall polygons in its set;
    - the fluid area of each of the six faces by Green's theorem over the
      fluid edge intervals and the traces of its wall polygons on the face,
      oriented by the polygon winding: no sign test, so a wall lying in a
      face plane needs no special case;
    - volume and centroid from the high faces and the walls (coordinates
      local to the cell);
    - an inclusion (a closed solid inside the cell) belongs to the piece
      around it, found by a ray in +x from its rightmost vertex to the
      nearest wall, or to the +x face.
    A cell with one piece takes the cell totals.

## Checks you can make

- **Divergence theorem per cell.** The aperture-weighted face (edge) normals
  equal the wall normal integral. In 3-D:
  `[(face_x(i+1) - face_x(i)) dy dz, (face_y(j+1) - face_y(j)) dz dx, (face_z(k+1) - face_z(k)) dx dy] = wet_normal`.
  The residual is at round-off.
- **Covering sums.** Over a block covering the body, the solid volumes
  (areas) add up to the body volume (area).
- **Small cells.** Flag the cells whose fraction is below a threshold, for
  example `cut%alpha < 1.e-4`, the instability threshold the literature
  records for coarse cut cells. The small-cell treatment itself belongs to
  the solver.

![the turned unit cube of the example, opaque, zoomed on three of the cells whose fluid fraction is below 1e-2, outlined in red: each cell lies almost entirely inside the cube and only one corner of it pokes out through a face; the smallest is labelled alpha 1.58E-05](../../examples/images/ref-cut-cells-checks.svg){.fossil-figure}

On the turned cube of the example (6³ cells of 0.25) the covering sum is 5.55·10⁻¹⁶ and the largest divergence
residual 3.22·10⁻¹⁶. The 5 small cells are cells the cube almost fills: only a corner sticks out through a face, the
smallest, cell (5, 4, 3), with a fluid fraction of 1.576·10⁻⁵ (`fig_ref_cut_cells_checks.f90`).

![a plate 0.5 long and 0.03 thick, half a cell, translating and turning through 16 by 10 cells of 1/16 over 24 frames: the cut cells shaded by fluid fraction, two gold dots at the piece centroids of every cell the plate splits in two, and the cells holding a piece below 0.1 of the cell outlined in red; the caption counts them and gives the covering sum of the frame](../../examples/images/anim-ref-cut-cells-plate.svg){.fossil-figure}

The checks on a moving body, frame by frame (`anim_ref_cut_cells_plate.f90`): the plate cuts 9 to 18 cells, splits 3
to 8 of them into two pieces and leaves up to 7 cells with a piece below 0.1 of the cell; the covering sum stays
within 7.98·10⁻¹⁷ of the plate's area in every frame.

## Results

| Case | Result |
|---|---|
| 3-D: 9 sign-benchmark bodies, 13³ centres | parity inside = solid angle off the surface, every line consistent |
| 3-D: unit cube with lines along faces, through edges and vertices | no odd lines, parity right, fractions exact |
| 2-D: offset, rotated squares, 128-gon, uniform and stretched | fraction, moment, apertures, wetted length ≤ 5·10⁻¹⁵ against clipping |
| 2-D: diamond with vertices and sides through grid nodes | exact |
| 2-D: thin solid strip across a row | 2 fluid pieces per cell, areas right |
| 2-D: island, hole, island in a hole inside one cell | 1 piece, right fraction |
| 3-D faces: rotated cube, pyramid, 32-face cone | face fractions and moments ≤ 3.4·10⁻¹⁵ against half-space clipping |
| 3-D cells: axis-aligned box | fractions and moments ≤ 5·10⁻¹⁶ |
| 3-D cells: half-space | fractions ≤ 1.7·10⁻¹³ against the closed-form volume |
| 3-D cells: rotated cube, pyramid, cone | covering sums exact, per-cell divergence residual ≤ 2·10⁻¹⁶ |
| 3-D pieces: prisms (thin wedge rotated about z, square tube, L-prism), 410 cells | equal to the 2-D pieces of the section on every cell; 32 cells with 2 pieces |
| 3-D pieces: tube hole narrower than the cell; passage cutting a cell corner | 1 piece (the cell-boundary fluid alone would give 2) |
| 3-D pieces: thin plate; rods along x, y, z; void in the solid; hollow body; inclusion | 2; 1; 1; 2; 1 |
| 3-D pieces: cube on a grid through its faces, edges and corners | equal to 2-D, no cell flagged |
| 3-D pieces: bunny (open, 69451 facets), 24³ | holed cells flagged; no other cell with fluid and 0 pieces |
| 3-D piece measures: sums over the pieces vs cell totals, 16922 cells (115 with 2+ pieces; tilted wedge, cone, tube, L-prism generic and grid-aligned, plate, hollow box) | volume 5.5·10⁻¹⁵, moments 5.5·10⁻¹⁴, wall terms 3.8·10⁻¹⁵, face fractions 1.6·10⁻¹⁴ |
| 3-D piece measures: prisms vs the 2-D pieces of the section, 410 pieces | volume, first moments, ±z face fractions ≤ 1.8·10⁻¹⁴ |
| 3-D piece measures: thin plate; hollow box; inclusion in the void, above a plate, beside a wall | exact volumes (0.5 and 0.499 of the cell; 98 and 1; 0.992; 0.782; 0.792) |
| `orient3d`: near-coplanar grid (6144 points), coplanar dyadic points | exact signs and zeros (the naive determinant is wrong on 2022) |
| line crossings, device path vs CPU | 0 mismatches (device path run on the host) |
| 2-D cut cells, `cut_cells_device` vs `cut_cells`: the sections of 9 benchmark bodies and 4 curves up to 1024², with and without intervals (RTX 4070) | 0 mismatches, every output bitwise (`fossil_test_curve_bench 256`) |
| 2-D cut cells resident (`resident=.true.`, then `update_host`) vs `cut_cells`, same cases; the classification and prefix sums of the device path vs the host's, integer for integer | 0 mismatches (`fossil_test_curve_bench 256`, `fossil_test_cut_cells_device_twin`) |
| a moving ellipse and NACA 0012, 3 steps of move + tag + cut, device tags and resident cut vs the CPU at every step | 0 mismatches (`cut cells 2-D chain device vs serial mismatches`) |
| 3-D pieces: cube-01, an octahedron and 40 lattice tetrahedra on grids of step 1/8, 1/4, 1/16 (1.57 million cells) | 19 cells untrusted (947 before issue #70), none outside [0, 1]; per-piece closure 3.5·10⁻¹⁴, piece sums 1.4·10⁻¹⁴, block volume 1.6·10⁻¹⁵ |
| node parities along x, y, z on a grid-aligned cube, an octahedron, 20 lattice tetrahedra (2-D: square, diamond, 20 triangles) | 0 of 108 086 and 0 of 6358 disagree (108 and 45 before issue #69) |
| centre inside flags within 3 ulp of a tilted face: tetrahedron (3-D), triangle (2-D), quadruple-precision reference | 0 of 3168 and 0 of 3200 wrong (249 and 295 before issue #69) |

## Limitations

- **Zero-volume pieces.** A wall lying on the high plane of a cell belongs to
  that cell, so the open face (aperture 1) bounds a fluid piece of zero
  volume. It is counted, as in 2-D, with `piece_volume` 0 (centroid at the
  cell centre); a solver can drop pieces with no volume.
- **Flagged cells (`pieces = -1`).** Cells a boundary edge of the mesh
  overlaps, cells whose fraction is outside [0, 1] (an open surface corrupts
  the face parities beyond the cells it flags), cells whose edge states
  do not close, and (issue #70) multi-piece cells whose split does not
  close: the face fractions of the pieces must add up to the plane
  sections, their volumes to the cell's, and every piece must close per
  axis.
- **Points exactly on the surface.**
  - The inside flag of a cell centre is the parity of the crossings along
    its x line before it. Since issue #69 "before" is exact: with
    n = (B − A) × (C − A), orient3d(A, B, C, P) = −n_x (c − t), so one
    `orient3d` and the crossing's entering flag decide it (2-D: one
    `orient2d`). A centre on the wall counts as past the crossing (taken at
    x + ε). Before #69 the rounded crossing position decided it, and on a
    tilted tetrahedron 249 of 3168 centres within 3 ulp of a face landed on
    the wrong side (295 of 3200 in 2-D).
  - **One global perturbation** (issue #69). Every grid point is taken at
    P + (ε₁, ε₂, ε₃), ε₁ ≫ ε₂ ≫ ε₃: the line kernels break ties on a
    projected edge or vertex in that order (on y lines x before z; before
    #69 each axis ordered its own pair cyclically, with no global order),
    and a point on a facet's plane takes the sign of the perturbed
    `orient3d`, minus the first nonzero component of the normal. The
    parities of a grid node along x, y and z are then the point-in-solid
    test of one perturbed point and agree everywhere: 0 disagreements over
    108 086 nodes of a grid-aligned cube, an octahedron and 20 lattice
    tetrahedra (108 before), 0 over 6358 nodes in 2-D (45 before). The 2-D
    sections of the 3-D cut cells keep their own (u, v) order.
  - **Pieces on bodies through grid nodes** (issue #70). Before #70, on
    tetrahedra with vertices on the grid nodes, the piece classification
    flagged up to 38 cells per body (947 of 1.57 million cells over 42
    bodies and three grids): a mesh edge or vertex on a cell edge gave one
    piercing per incident facet. The edges now take their exact crossings,
    and the face traces are read on the perturbed plane: a polygon lying
    in the face plane bounds the fluid only if its solid lies beyond the
    face, and then not along the cell edges; a facet ending on the plane
    counts there only if the facet across its edge crosses the plane. 19
    cells remain flagged, where a facet touches the cell only along a
    segment of a cell edge (its clipped polygon is empty, so its trace on
    the face is missing); the closure test flags them, none is wrong. The
    "fractions outside [0, 1]" seen before were 1 + 4·10⁻¹⁶, round-off of
    the face sums, inside the 10⁻¹² tolerance.
- **Constructed points.** The 3-D sections and the clipping decide on the
  rounded crossing points, consistently across kernels. A body passing
  exactly through grid lines at non-representable points can then give a
  face that does not close; it is flagged, never silently wrong.
- **Watertight input.** Cells whose boundary does not close (open curves or
  surfaces) are flagged in `closed` and get no measures.
- **On the device: the 2-D cut cells only.** `curve%cut_cells_device` and the line crossings have device paths; the
  3-D cut cells, the swept geometry and the piece links are CPU only (issue #100 plans them).
