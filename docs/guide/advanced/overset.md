# Overset assembly kernels

Overset (Chimera) grids cover a domain with overlapping blocks. Assembling them is geometry between blocks and
bodies: which points lie inside a body (**holes**), which take their values from another block (**receptors**), from
which cells (**donors**) and with what quality, and which have no valid donor (**orphans**). Issue #72 (roadmap #68,
WP3) adds the kernels for it on top of the [curvilinear blocks](./curvilinear-blocks) of WP2, in the module
`fossil_overset`.

::: info Design and reference (issue #72)
FOSSIL gives kernels and measures. The **assembly policy** (explicit or implicit hole cutting, fringe width,
priorities) stays with the caller: a reference assembly built from the kernels lives in the tests, the cookbook and
the validation suite, not in the library. Sources are the llm-wiki bundle `overset-grid-methods` (pages cited by
path); statements from general knowledge are marked *(gk)*.
:::

## At a glance

![an O-grid of 32 by 12 cells around a unit cylinder crossing a Cartesian background of 60 by 32 cells: per frame, the background holes inside the cylinder in red, the background receptors around them in teal, the O-grid's outer receptors in yellow, and ember squares on the orphans where the O-grid's outer boundary leaves the background](../../examples/images/anim-overset.svg){.fossil-figure}

Each frame moves the O-grid and recomputes the assembly with the kernels of this page (`wall_cutter_2d`,
`hole_flags_2d`, `donor_candidates_2d`, `classify_2d` with the explicit rule and 3-point stencils): 137 to 141
background holes and 56 to 58 background receptors per frame. Orphans appear where the O-grid's outer boundary leaves
the background: 9 in the first frame, 3 while it grazes the top, 10 in the last (`anim_overset.f90`).

## Principles

1. **Points are vertices.** Every kernel works on the vertices of a set of `curvilinear_block_t` (2-D:
   `curvilinear_block_2d_t`), and the cells of those blocks are the donors. Cell-centred data (Xall, ADAM) use the
   blocks of their cell centres: the caller passes the centres with one ghost layer, or `cell_centre_block` builds
   the outer layer by linear extrapolation from the boundary face centres.
2. **Exact where it decides.** Hole flags come from the fast winding number on the real walls (exact on closed walls,
   robust on open or dirty ones); the hole map's cut elements are exact (`BOX_CUT`), its lines of sight and the
   X-rays use exact crossings. Everything approximate is measured against the exact flags.
3. **The standard output** (`definitions/domain-connectivity.md`, `definitions/point-classification.md`): a mask ι
   per point, 1 for a field point, −k′ for a receptor interpolated from block k′, 0 for a hole (Lee & Baeder); per
   receptor the donor block, the stencil corner, ξ, the width per direction, a low-order flag and the donor quality.
   Fringe levels and orphan causes are separate arrays, so ι keeps its meaning.
4. **CPU and device bit-exact** for the batched queries (hole flags, donor candidates), as every device path of
   FOSSIL.
5. **Deterministic.** Blocks, points and candidates are visited in index order; the result does not depend on the
   order in which blocks are passed beyond their numbering.

## Cutters

A cutter is the wall of a body as a `surface_stl_object` (2-D: a `curve_object`), built from the wall faces of the
blocks (`wall_cutter`): the faces covered by `TOPO_BOUNDARY` records with the given name (default `wall`).

- Each bilinear face is split into two triangles along its shorter diagonal; the distance between the triangles and
  the bilinear face is reported (`gap`), so the caller knows how far the cutter is from the grid's own wall.
- Every triangle (segment) is oriented with the block interior, the fluid, on its outer side: the solid is inside
  the cutter, as the winding number needs.
- **Overlapping wall patches** of different blocks (collar grids, the two caps of a Yin–Yang sphere) cover the same
  wall twice; a plain union counts the overlap twice in the winding number *(gk/inference)*. They are trimmed: an
  element of the coarser patch whose vertices all lie within the tolerance of the finer patch is dropped (the finer
  patch wins). Elements straddling the edge of the other patch stay: the cutter may keep a thin double-covered strip,
  never a gap.
- **Shared walls**: two blocks whose wall patches overlap share a wall. In Lee & Baeder's pair matrix they overlap
  *without mutual cutting*: they interpolate from each other but never cut each other's points
  (`methods/implicit-hole-cutting.md`).

## Hole cutting

Three families (`methods/hole-cutting.md`), all measured against the first:

| Kernel | What decides | Exact? |
|---|---|---|
| `hole_flags` | the fast winding number of the cutter at each vertex (CPU or device); `offset=` also blanks the points within a signed distance of the wall (PEGASUS 5 OFFSET, Noack) | yes on closed walls |
| `hole_map` | a Cartesian map over the cutter: elements cut by the wall (exact `BOX_CUT`), the rest painted from the corners; a point in a cut element looks along a segment to the centre of a neighbouring uncut element and takes its status, flipped by the parity of the exact wall crossings | yes, with the parity; PEGASUS 5's rule (status of a visible element) fails inside bodies thinner than the map, which the tests show |
| `xray_flags` | object X-rays (Meakin): the exact crossings of the wall along a raster of lines parallel to an axis; a point takes the parity of the crossings below it on the nearest ray | no: wrong only within the raster spacing of the wall |

A point exactly on the wall is a field point (signed distance 0 is not inside), so a body-fitted block never blanks
its own wall. Which cutter cuts which block is the caller's choice (`cuts=`, default every block). Pass the blocks
that contributed the cutter's patches as not cut when they move with a pose: their wall vertices, mapped to the world
frame and back, land off the wall by round-off and would be blanked.

![the rear half of a lens 0.008 thick at most and of chord 1 in a 241 by 41 background block, the vertical scale stretched 6 times: the background vertices inside the lens in teal (exact flags), six of them in ember near u = 0.78, the holes PEGASUS 5's line of sight misses on a 20 by 20 hole map, and in gold the vertices where X-rays through 41 rays differ from the exact flags, at u = 0.79 and at the trailing edge](../../examples/images/ref-overset-holes.svg){.fossil-figure}

On this lens the exact flags find 433 holes. The hole map with the crossing parity agrees with them everywhere; with
PEGASUS 5's line of sight it leaves 12 vertices unresolved and misses 12 holes, at u = 0.21 and u = 0.78: the rule
takes the status of a visible uncut element, and inside a body thinner than the map elements there is none. X-rays
through 41 rays differ from the exact flags at 7 vertices (`fig_ref_overset_holes.f90`).

## Donor candidates

`donor_candidates` locates every vertex of every block in every other block the pair matrix allows (WP2's robust
location, per block pair), and keeps **all** candidates (`methods/donor-search.md`), each with:

- the donor block, cell and ξ, and the location status;
- the donor cell volume and the receptor volume (the mean volume of the cells around the vertex), and their ratio;
- the certified min J of the donor cell over its mean (WP2's validity).

PEGASUS 5's cell-difference parameter is **not** used as printed: it is negative when the donor is larger
(`methods/overlap-optimization.md`). The wall distance of every vertex (`wall_distance_many`) comes with it, for the
wall-distance criterion. The pair matrix is three-valued: two blocks **overlap** if a boundary vertex of one finds a
donor in the other, or one sits inside the other; blocks with a shared wall **overlap without cutting**. On the
device the blocks are entered once (`enter_device`) and every block pair runs the same `_oac` kernels as the CPU.

## Selection, validity and the fixed point

`classify` turns holes and candidates into the mask (`definitions/point-classification.md`):

- **explicit**: holes from the cutters; the points next to holes (the fringe) and the points on overlap boundaries
  (faces named `overlap`) must receive;
- **implicit** (Lee & Baeder, `methods/implicit-hole-cutting.md`): a point receives when a candidate's donor cell is
  smaller than its own cells; overlap-boundary points always receive; points within a wall-normal immunity layer
  never receive;
- **wall distance** (Oh et al. 2024, `methods/wall-distance-hole-cutting.md`): `donor_candidates(..., wall=)`
  returns the exact wall distance of every vertex, so a caller can classify by it before the donor search; FOSSIL
  does not fix the threshold (a policy).

A receptor is **valid** when its donor stencil of width q (WP2's stencils) holds no hole, and, for explicit
interpolation, no receptor either (Chesshire & Henshaw); an invalid candidate gives way to the next one. Validity
depends on the neighbours, so the classification is iterated to a fixed point (CMPGRD), preferring the fewest
receptors. **Circular interpolation** (receptors in each other's stencils) is detected and reported (Sherer & Scott).

## Fringes, orphans, overlap width

- **Fringe layers** (`definitions/fringe-point.md`) n_f deep, by SUGGAR's rules: node-centred (a cell with a hole
  vertex makes its other vertices INNER 1, and so on), or cell-centred on the cell-centre blocks; overlap-boundary
  points are OUTER 1, their cell neighbours OUTER 2.
- **Orphans** (`definitions/orphan-point.md`) carry their cause: no candidate (insufficient overlap), the stencil hits
  a hole, the donor rejected near a boundary, or a low-quality donor whose own vertices receive (Roget & Sitaraman).
- **Overlap width** as Chesshire's d: the distance between the innermost receptors of two blocks; d = 0 makes the
  interpolation system singular, and a stencil one point short loses one order of accuracy.

## Wall projection for viscous walls

Two overlapping wall patches on a curved wall are each discretized linearly and do not coincide (Rogers et al. 2003,
`methods/inter-grid-interpolation.md`). On a concave wall the near-wall receptors fall outside the donor grid
(orphans); on a convex wall they find donors much farther from the wall. `wall_projection` moves each wall vertex of
the receptor block onto the donor block's wall (FOSSIL's closest point on it) and carries the shift up the
wall-normal grid line, decaying linearly to zero at a maximum distance. The shifted points are used for donor
search and weights only; the grid is never moved.

![an annulus unrolled, arc length across and 40 times the distance from its outer circle r = 2 upward, over 0.67 pi to 0.83 pi: three cells of the coarse block's two outer layers, its wall chords (gold) sagging below the circle; the fine block's wall on the circle, its vertices between the coarse vertices (ember) outside the coarse block, each moved by a short gold segment to its teal image on the coarse chord](../../examples/images/ref-overset-projection.svg){.fossil-figure}

Two blocks share a concave wall: a fine one, 41 vertices over 0.9 π, and a coarse one, 9 over 0.4 π. Without
projection, 12 fine wall vertices inside the coarse block's range have no donor; `wall_projection_2d` moves them by
6.104e-3 at most, the chord sag, and every one finds a donor (`fig_ref_overset_projection.f90`).

## Moving grids

`methods/moving-overset-grid.md`:

- blocks moving with a body keep their coordinates in the body frame and take a pose, so their search structures are
  built once (SUGGAR);
- the search is incremental: each point starts its walk from its previous donor (WP2's walk, equal to the robust
  path);
- **exposed points**, holes or receptors at step n that become field points at step n + 1, are listed (Henshaw &
  Schwendeman give them old values by interpolation);
- the check D^{n+1} ∩ H^n = ∅ (Meakin & Suhs: no hole becomes a field point without first being a receptor) is
  reported, which bounds how far a hole moves per step.

## Connectivity output

`connectivity_t` holds the mask, the fringe levels and orphan causes per block, and the receptor list (receptor block
and vertex, donor block, stencil corner, ξ, width per direction, low-order flag, quality); `save_connectivity` writes
it as plain text, one receptor per line, in the form the founding codes share (`definitions/domain-connectivity.md`).

## Example

A full assembly of the cylinder of the [cookbook recipe](/manual/cookbook#assemble-an-overset-grid) (an O-grid of
32 × 12 cells to r = 2 over a Cartesian background of spacing 0.15), with the kernels the recipe skips: the hole
offset and the distances, the pair matrix, the donor quality, the implicit rule, two fringe levels, the overlap width
and the connectivity file.

<<< @/examples/snippets/ref_overset.f90

<<< @/examples/output/ref_overset.ansi{ansi}

- **Holes.** The O-grid's wall cuts the background only (`cuts=[.false., .true.]`): 137 holes. `offset=0.3` also
  blanks the vertices within 0.3 of the wall, 241 in all. `distances=` returns the signed distance of every vertex:
  the nearest field vertex lies 0.0107 from the wall.
- **Pairs and donors.** `overlap_matrix_2d` gives `PAIR_OVERLAP` (1) for the two blocks, printed column by column.
  `donor_candidates_2d` keeps every candidate: all 429 O-grid vertices have a background cell, 416 background
  vertices an O-grid cell. The donor/receptor volume ratio of the O-grid vertices runs from 0.7067 to 1.3286.
- **Explicit rule.** With the offset holes and 3-point stencils, 33 O-grid vertices (its overlap boundary) and 72
  background vertices receive, no orphan. The fixed point settles in 2 iterations, with no circular pair.
- **Implicit rule.** Without holes, a vertex receives when a smaller donor cell exists: 191 and 132 receptors, a
  different assembly from the same candidates, after 6 iterations.
- **Fringes.** `levels=2` gives the two INNER levels around the holes (72 and 80 background vertices) and the two OUTER
  levels on the O-grid (33 and 33). `overlap_width_2d` gives Chesshire's d = 0.5088 between the innermost receptors.
- **Output.** `save_connectivity` writes the counts per block, then one receptor per line: receptor block and vertex,
  donor block, stencil corner, ξ, width per direction, low-order flag and quality (the line is cut at 100 characters
  here).

## API

```fortran
call wall_cutter(blocks, cutter, gap=gap, shared=shared)                  ! 2-D: wall_cutter_2d -> curve_object
call hole_flags(blocks, cutter, holes, offset=0, cuts=cuts, device=.true.) ! exact; or hole_map / xray_flags
call blocks(b)%build_search()                                             ! every block, then (device) enter_device
call overlap_matrix(blocks, pair, shared=shared)
call donor_candidates(blocks, donors, pair=pair, wall=cutter, device=, shifts=, poses=, previous=)
call classify(blocks, donors, masks, SELECT_EXPLICIT, q, holes=holes, layers=, immune=, explicit_interpolation=)
call fringe_levels(blocks, masks) ; call overlap_width(blocks, masks, d)
call wall_projection(blocks, a, b, shift, max_distance)                    ! shifts(a, b)%v = shift
call exposed_points(masks_old, masks_new, exposed, violations)
call build_connectivity(masks, donors, q, 3, conn) ; call save_connectivity(conn, 'connectivity.txt')
```

The mask `masks(b)%iota(v)` is 1 (field), 0 (hole) or −k (receptor from block k) at the linear vertex index v (i
fastest); `masks(b)%donor(v)` is the chosen candidate in `donors(b)`, `orphan(v)` an `ORPHAN_*` cause, `fringe(v)` a
`FRINGE_*` level. A cell-centred solver passes the blocks of its cell centres.

## Measured

- `fossil_test_overset_*`: holes equal to exact inside tests (2-D, 3-D, device = CPU); the hole map with the crossing
  parity equal to the exact flags, PEGASUS 5's line of sight missing 12 holes in a lens 0.008 thick on 0.06 elements;
  candidates equal to a brute force; the cylinder's receptor rings exact; circular pairs counted one per row; each
  orphan cause produced on purpose; overlap width of concentric annuli 0.5; the projection removes the concave-wall
  orphans and puts convex-wall donors on the wall (η below 1e-7); the incremental search equal to the search from
  scratch, bit for bit; a 0.6 h step exposes points without a hole becoming a field point, a 3.6 h step does.
- `fossil_val_overset`: interpolation through the assembly converges at the stencil width, q = 2, 3, 4 at orders 2.0,
  3.2, 4.4 (a smooth function) and 2.0, 3.2, 4.5 (an isentropic vortex across the overlap), no orphan; a spherical
  shell of 64 × 32 × 16 cells in a 49³ background assembles in about 1.2 s on one core (holes 1.1 s, candidates 0.09
  s, classification 0.05 s).
- GPU (RTX 4070): hole flags and donor candidates equal to the CPU, bit for bit.

## Proved and measured

| Claim | Kind |
|---|---|
| hole flags on closed walls | exact (winding number), CPU = device bit for bit |
| hole map with crossing parity, X-rays on their rays | exact; X-rays between rays measured (errors within the raster spacing) |
| donor candidates | WP2's location; equal to a brute force; CPU = device bit for bit |
| classification fixed point | Jacobi iteration, independent of the block order (tested); circular pairs detected |
| interpolation order through the overlap | measured: the rate equals q (fixed overlap); the order loss of a stencil one point short needs a discretization, left to the solver |

## Left out

- Conservative-interface geometry (`methods/conservative-overset-interface.md`): a later WP decides.
- Overlap optimization beyond the measures, and parallel (MPI) assembly: the caller's.
- Stencils crossing into a third block: reported low-order, not followed.
