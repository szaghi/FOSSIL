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
    and normal.

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
| line crossings, device path vs CPU | 0 mismatches (device path run on the host) |

## Limitations

- **No separate fluid pieces in 3-D yet.** Split and tunnel cells are
  measured as a whole; the 2-D cut cells report pieces.
- **Points exactly on the surface.**
  - The inside flag of a cell centre comes from the crossings along its x
    line. Whether a crossing lies before or after the centre is decided on
    the rounded crossing position, so a centre within round-off of the wall
    can land on either side.
  - The tie-break of the line crossings orders the two transverse
    perturbations cyclically per axis. That is exact for the crossing counts
    of each line, but not one global perturbation.
  - Exact agreement across axes for points on the surface would need an
    exact `orient3d` with symbolic perturbation, which FOSSIL does not have.
- **Constructed points.** The 3-D sections and the clipping decide on the
  rounded crossing points, consistently across kernels. A body passing
  exactly through grid lines at non-representable points can then give a
  face that does not close; it is flagged, never silently wrong.
- **Watertight input.** Cells whose boundary does not close (open curves or
  surfaces) are flagged in `closed` and get no measures.
- **CPU only for cut cells.** The line crossings have a device path; the GPU
  gate has not run yet.
