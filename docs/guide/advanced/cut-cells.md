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

## Usage

```fortran
program cut_geometry
use fossil, only : curve_object, cut_cells_2d_t, surface_stl_object, cut_cells_3d_t, cut_cells_3d
use penf,   only : I4P, R8P
implicit none
type(curve_object)       :: c
type(cut_cells_2d_t)     :: cut2
type(surface_stl_object) :: body
type(cut_cells_3d_t)     :: cut3
real(R8P)                :: e(9), xc(8), frac(6, 8, 8, 8)
logical                  :: inside(8, 8, 8)
integer(I4P)             :: k

do k = 1, 9
   e(k) = (k - 1) / 8._R8P ! cell edges
enddo
xc = 0.5_R8P * (e(1:8) + e(2:9)) ! cell centres

! 2-D: a square, cut cells of an 8 x 8 block
call c%adopt_loops(u=[0.3_R8P, 0.7_R8P, 0.7_R8P, 0.3_R8P], v=[0.3_R8P, 0.3_R8P, 0.7_R8P, 0.7_R8P], loop_sizes=[4_I4P])
call c%cut_cells(x=e, y=e, cut=cut2)
print '(A,F8.5)', 'fluid fraction of cell (3, 3): ', cut2%alpha(3, 3)

! 3-D: wall fractions at the cell centres, then the cut cells
call body%load_from_file(file_name='body.stl', guess_format=.true.)
call body%cell_wall_fractions(x=xc, y=xc, z=xc, frac=frac, inside=inside)
call cut_cells_3d(surface=body, x=e, y=e, z=e, cut=cut3)
print '(A,I0)', 'fluid pieces of cell (4, 4, 4): ', cut3%pieces(4, 4, 4) ! -1: cell not trusted
if (cut3%pieces(4, 4, 4) > 0) print '(A,ES12.4)', 'volume of its first piece: ', &
   cut3%piece_volume(cut3%first_piece(4, 4, 4))
endprogram cut_geometry
```

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
    - the clip flags every vertex it makes with its cell plane, so the
      vertices on a cell edge (two planes) are the edge's piercings, entering
      or leaving the solid by the sign of the normal along the edge;
    - the parts are joined where polygons share a vertex (the clipped points
      of a shared mesh edge are bit-identical), at a piercing, at a fluid
      corner, and where a solid island traced on a face is surrounded by
      fluid (a ray from the island to the nearest trace or face edge);
    - every joined set that touches the cell boundary is one piece; a closed
      wall set inside the cell adds a piece only if it encloses a void.
  - So a passage narrower than the cell, whose fluid meets the cell boundary
    in two separate disks, is one piece, and a thin plate across the cell
    leaves two.
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
| node parities along x, y, z on a grid-aligned cube, an octahedron, 20 lattice tetrahedra (2-D: square, diamond, 20 triangles) | 0 of 108 086 and 0 of 6358 disagree (108 and 45 before issue #69) |
| centre inside flags within 3 ulp of a tilted face: tetrahedron (3-D), triangle (2-D), quadruple-precision reference | 0 of 3168 and 0 of 3200 wrong (249 and 295 before issue #69) |

## Limitations

- **Zero-volume pieces.** A wall lying on the high plane of a cell belongs to
  that cell, so the open face (aperture 1) bounds a fluid piece of zero
  volume. It is counted, as in 2-D, with `piece_volume` 0 (centroid at the
  cell centre); a solver can drop pieces with no volume.
- **Flagged cells (`pieces = -1`).** Cells a boundary edge of the mesh
  overlaps, cells whose fraction is outside [0, 1] (an open surface corrupts
  the face parities beyond the cells it flags), and cells whose edge states
  do not close.
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
  - **Pieces on bodies through grid nodes.** On tetrahedra with vertices
    on the grid nodes (and so edges on grid lines and faces through them),
    every face and cell closes and the solid volume is right to 1e-15, but
    the piece classification flags 0 to 38 cells per body as untrusted
    (`pieces = -1`), and up to 2 cells per body get a volume fraction
    outside [0, 1]. Unchanged by issue #69 (identical before it): a
    follow-up.
- **Constructed points.** The 3-D sections and the clipping decide on the
  rounded crossing points, consistently across kernels. A body passing
  exactly through grid lines at non-representable points can then give a
  face that does not close; it is flagged, never silently wrong.
- **Watertight input.** Cells whose boundary does not close (open curves or
  surfaces) are flagged in `closed` and get no measures.
- **CPU only for cut cells.** The line crossings have a device path; the GPU
  gate has not run yet.
