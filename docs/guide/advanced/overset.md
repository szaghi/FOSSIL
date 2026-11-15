# Overset assembly kernels

Overset (Chimera) grids cover a domain with overlapping blocks. Assembling them is geometry between blocks and
bodies: which points lie inside a body (**holes**), which take their values from another block (**receptors**), from
which cells (**donors**) and with what quality, and which have no valid donor (**orphans**). Issue #72 (roadmap #68,
WP3) adds the kernels for it on top of the [curvilinear blocks](./curvilinear-blocks) of WP2, in the module
`fossil_overset`. Issue #106 (roadmap #92, WP4) adds the geometry of a
[conservative interface](#conservative-interface) on a rectilinear background, in `fossil_overset_patch`.

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

![a close-up of the cylinder assembly where the O-grid (tan curved lines) crosses the Cartesian background (grey lines): each background vertex carries its iota, 0 (red) inside the cylinder, -1 (teal) on the ring around the holes that receives from the O-grid, 1 (grey) elsewhere; the O-grid's outer vertices that receive from the background carry -2 (gold); one O-grid cell donating to a background vertex is filled teal, one background cell donating to an O-grid vertex gold](../../examples/images/ref-overset-mask.svg){.fossil-figure}

The mask of the cylinder assembly of the Example (explicit rule, holes without offset, 3-point stencils): 137
background holes, 56 background receptors from the O-grid and 33 O-grid receptors from the background. Background
vertex (27, 24) receives from O-grid cell (3, 1); O-grid vertex (4, 13) from background cell (32, 28)
(`fig_ref_overset_mask.f90`).
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
  element of the coarser patch (larger mean element size over the overlap) is dropped only if the finer patch covers
  it, that is its vertices and its centroid (in 2-D, its ends and midpoint) all lie within the tolerance of the finer
  patch; the tolerance defaults to 0.1 of the element size. The elements straddling the edge of the finer patch stay,
  so no gap opens: the cutter keeps a double-covered strip at most one coarse element wide, which raises the winding
  number inside by the doubly covered solid angle over 4π (issue #91; no gap is the safe side for hole cutting).
  `strict=.false.` restores the centroid-only rule, which also drops a straddling element whose centroid lies on the
  finer patch and leaves a gap there.
- **Shared walls**: two blocks whose wall patches overlap share a wall. In Lee & Baeder's pair matrix they overlap
  *without mutual cutting*: they interpolate from each other but never cut each other's points
  (`methods/implicit-hole-cutting.md`).

![two partial annuli around the unit cylinder: a fine one of 24 cells from 0 to 200 degrees on top, its wall in gold, and a coarse one of 11 cells from 153 to 373 degrees below, its kept wall segments in ember drawn just outside the circle; the two coarse segments dropped by the trimming are drawn in grey just inside it; on the left the kept coarse segment overlaps the fine wall from 193 to 200 degrees, on the right the kept coarse segment from 353 to 373 degrees overlaps it from 0 to 13, no gap; teal ticks on every kept segment point to the fluid](../../examples/images/ref-overset-cutter.svg){.fossil-figure}

Two blocks share the cylinder's wall (`shared(1, 2)` is true). `wall_cutter_2d` keeps the fine patch whole and 9 of
the 11 coarse segments: the two from 153 to 193 degrees lie on the fine patch and go; the segment from 193 to 213 has
its end at 213 off the fine patch and stays (the wall from 193 to 200 is covered twice), and so does the one from 353
to 373, whose end at 353 is off it (0 to 13 covered twice). The double cover is 20 degrees, the gap 0
(`fig_ref_overset_cutter.f90`; `fossil_test_overset_cutters` checks every angle of the coarse range exactly). The
centroid rule (`strict=.false.`) dropped the segment from 353 to 373, whose midpoint lies on the fine patch, and left
the wall from 353 to 360 uncovered. Both strips are within one coarse segment of the patch edge.

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

![the outer ring of the cylinder's O-grid (tan) at angles -11.25, 0 and 11.25 degrees over the background of spacing 0.15 (grey) and a refined patch of spacing 0.075 (teal); each of the three outer vertices (red dots) has two candidates, a background cell filled brown and a patch cell filled gold, labelled with the donor / receptor volume ratios 0.71 and 0.18](../../examples/images/ref-overset-donors.svg){.fossil-figure}

With a third block, a patch of spacing 0.075 over [1.65, 2.4] × [-0.6, 0.6], `donor_candidates_2d` returns 411 O-grid
vertices with one candidate and 18 with two. For the outer-ring vertices the background cell is 0.7067 times the
receptor's cells, the patch cell 0.1767 times; the vertex at angle 0 lies on a background cell edge (ξ = 1) and the
lower cell is taken (`fig_ref_overset_donors.f90`). Which candidate wins is the classification's choice, below.

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

![the cylinder assembly twice, from the same candidates: left the explicit rule with the holes cut at offset 0.3 (red), the ring of background receptors around them (teal) and the O-grid's outer ring receiving (gold); right the implicit rule without holes, a wider band of background receptors (teal) inside the O-grid and the O-grid's outer five or so rings receiving (gold) where the background cells are smaller](../../examples/images/ref-overset-selection.svg){.fossil-figure}

The explicit rule gives 33 O-grid and 72 background receptors and settles in 2 iterations; the implicit rule, from the
same candidates, 191 and 132 after 6: each vertex receives from whichever block has the smaller cell there, so the
O-grid hands over to the background well inside its outer boundary (`fig_ref_overset_mask.f90`).

## Fringes, orphans, overlap width

- **Fringe layers** (`definitions/fringe-point.md`) n_f deep, by SUGGAR's rules: node-centred (a cell with a hole
  vertex makes its other vertices INNER 1, and so on), or cell-centred on the cell-centre blocks; overlap-boundary
  points are OUTER 1, their cell neighbours OUTER 2.
- **Orphans** (`definitions/orphan-point.md`) carry their cause: no candidate (insufficient overlap), the stencil hits
  a hole, the donor rejected near a boundary, or a low-quality donor whose own vertices receive (Roget & Sitaraman).
- **Overlap width** as Chesshire's d: the distance between the innermost receptors of two blocks; d = 0 makes the
  interpolation system singular, and a stencil one point short loses one order of accuracy.

![the cylinder assembly over a background cut short at x = -1.8, so the O-grid's outer boundary leaves it on the left: holes (brown) inside the cylinder, the background levels INNER 1 (teal) and INNER 2 (gold) around them, the O-grid levels OUTER 1 (teal) on its outer boundary and OUTER 2 (gold) one ring in, five ember orphans on the left where the outer boundary has no background cell, and a red segment marking the overlap width d between the nearest receptors of the two blocks](../../examples/images/ref-overset-fringes.svg){.fossil-figure}

With the background starting at x = −1.8, `fringe_levels_2d(levels=2)` gives 56 INNER 1 and 64 INNER 2 background
vertices and 33 OUTER 1 and 33 OUTER 2 O-grid vertices; 5 OUTER 1 vertices on the left are orphans with
`ORPHAN_NO_CANDIDATE`. `overlap_width_2d` gives d = 0.8400, equal to the brute-force minimum over the receptor pairs
(`fig_ref_overset_mask.f90`).

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

![the cylinder moving over the background of spacing h = 0.15 in one step, its wall before (grey circle, n) and after (red circle, n + 1): left a step of 0.6 h, right of 3.6 h. Background vertices after the step: holes brown, receptors teal; exposed points, receptors before and field points after, gold on the trailing side; on the right, holes that became field points without being receptors, ember](../../examples/images/ref-overset-exposed.svg){.fossil-figure}

The O-grid moves with a pose (`poses=` and `cutter_pose=`, searches built once, `previous=` walks from the last
donors). A step of 0.6 h exposes 8 background vertices, all receptors before: no violation. A step of 3.6 h exposes 53,
26 of which were holes: the step is too long for the fringe (`fig_ref_overset_exposed.f90`).

## Connectivity output

`connectivity_t` holds the mask, the fringe levels and orphan causes per block, and the receptor list (receptor block
and vertex, donor block, stencil corner, ξ, width per direction, low-order flag, quality); `save_connectivity` writes
it as plain text, one receptor per line, in the form the founding codes share (`definitions/domain-connectivity.md`).

![three receptors on the O-grid's outer ring (gold dots) over the background grid, each joined by teal lines to the nine vertices of its 3 by 3 background stencil, its corner cell filled brown, and labelled with its record: receptor block 1, vertices 398, 399, 400, donor block 2, corner cells (34, 23), (33, 26), (32, 28), xi and width 3 3](../../examples/images/ref-overset-connectivity.svg){.fossil-figure}

The records of three O-grid receptors, as `build_connectivity` stores them: receptor block and vertex (linear index, i
fastest), donor block, corner cell, ξ in it, width 3 × 3. Vertex 400 reads `1 400 2 32 28 0.086 0.408 3 3`: its stencil
is shifted towards the low side of the corner cell, where ξ lies (`fig_ref_overset_mask.f90`).

## Conservative interface

Interpolation couples the blocks of an overset grid without conserving anything across them: each block conserves its
own cells, the composite does not (Wang 1995, `methods/conservative-overset-interface.md`). Two families of methods
are called "conservative overset", and only the first one is conservative:

- **A. The patched interface** (Wang 1995). The minor grid's outer boundary cuts the major grid exactly: the major
  cells inside it become holes, the cells it crosses become cut cells with retained face fractions W_f and clipped
  volumes, and the boundary, split at the major grid's faces, gives the **smallest faces**. The two grids exchange
  fluxes through the smallest faces instead of interpolated values. The overlap then holds one solution, and the
  composite conserves whatever each grid conserves. Wang measured a mass imbalance of 6·10⁻⁶ % against 0.38 % for the
  original Chimera on a transonic bump channel.
- **B. Intersection weights** (Picard et al. 2021). The overlap is kept. A receptor cell takes the donor states
  averaged over its intersections with the donor cells, q_c = Σ w (q_c′ + ∇q_c′ · (x_cc′ − x_c′)) with
  w = |V_c ∩ V_c′| / |V_c| and x_cc′ the centroid of V_c ∩ V_c′. The transfer is second order on any mesh and needs no
  index space. Its authors call it *not* globally conservative, because it transfers states, not fluxes. They measured
  a mass defect of 0.1 to 0.6 % for a shock on a non-coincident interface.

FOSSIL computes the geometry of both. The flux, the reconstruction and the update stay with the solver. The **major
grid** (A) and the **donor grid** (B) are **rectilinear**: cell edges x, y (z) along each axis, strictly increasing,
stretched or not (Cartesian or AMR blocks). The minor grid and the receptors are curvilinear blocks. Both families
then reduce to the [exact cut cells](./cut-cells), so the geometry is exact and needs no tolerance. Curvilinear on both
sides would need the clipping of two trilinear cells against each other, which is left out.

![a quarter of the cylinder O-grid's outer boundary, the circle of radius 2 drawn as its 32 straight faces, over the Cartesian background of spacing 0.15: the background cells inside it are holes (teal), the cells it crosses are cut cells shaded by their solid fraction, and the boundary is split at the background's grid lines into smallest faces, alternately gold and ember by O-grid face, each with a short teal tick along its area vector, which points into the O-grid](../../examples/images/ref-overset-patch.svg){.fossil-figure}

The O-grid of the Example (32 × 12 cells to r = 2) patched into a background of spacing 0.15: 500 holes, 108 cut cells
and 139 smallest faces for its 32 outer faces (`fig_ref_overset_patch.f90`). Each tick shows the direction of
`sf_area`: out of the background piece, into the O-grid cell.

### A: the patched interface

`overset_patch_2d(blocks, minor, x, y, patch)` and `overset_patch_3d(blocks, minor, x, y, z, patch)` take the minor
blocks' faces covered by the `TOPO_BOUNDARY` records named `overlap` (`boundary_name=` to change it, or an explicit
face list `faces=`), orient each face with the minor region on its solid side, chain the faces into closed loops (2-D)
or check that every directed edge has exactly one reverse (3-D), and cut the major grid with the result:
`boundary%cut_cells(boundaries=.true.)` in 2-D, `cut_cells_3d(..., links=, walls=)` in 3-D, bit for bit what those
routines give. The cut labels its own wall pieces with the boundary element they lie on (2-D `piece_seg_src`, 3-D
`walls=`: one term per cell and facet, with its piece), so each smallest face knows its minor face with no geometric
matching. The outputs are:

- `hole(i, j[, k])`: the major cells inside the minor region, with no fluid piece; `cut`, the full cut-cell result:
  fluid fractions, apertures W_f (`piece_aperture` in 2-D, `piece_face` in 3-D), piece areas or volumes and centroids;
  in 3-D also `links`, the piece links a solver uses for the interior fluxes and the merges.
- Per smallest face k: its major piece `sf_piece(k)`, its minor face `sf_face(k)` (an index into `face_id`, which holds
  (block, i, j[, k], face)), the minor cell face `sf_minor(:, k)`, `sf_area(:, k)` and `sf_centroid(:, k)`; in 3-D
  also `sf_facet` and `sf_tri`, the triangle of the face's split it is a part of, and `sf_wet`, its area. Two CSR
  views index them: per major piece (`piece_sf_start`) and per minor face (`face_sf_start`, `face_sf`; in 2-D in order
  along the face).
- Per minor face f: `face_vector(:, f)`, its length (area) times its unit normal out of the minor cell, and
  `face_covered(:, f)`, minus the sum of its smallest faces' `sf_area`.

**Signs and units.** Lengths, areas and volumes are in the grid's own units. `sf_area` is the length (area) times the
unit normal **out of the major piece and into the minor cell**: a flux F · sf_area leaves piece `sf_piece` and enters
the minor cell `sf_minor`. The update is Wang's: a cut piece sums its regular face fluxes, each weighted by its W_f,
and its smallest faces; a minor boundary face sums the fluxes of its smallest faces. Two identities hold to round-off,
and the tests check both:

- every major piece closes: Σ_f W_f S_f n_f + Σ sf_area = 0, where the sum of a piece's `sf_area` is
  −`piece_wet_normal`, bit for bit when the piece has one wall chain (2-D) or one fan per extra polygon (3-D);
- every minor face inside the major box is covered: `face_covered` = `face_vector`. Where a face leaves the box,
  `face_covered` is shorter, and the minor grid's own boundary condition applies to the rest.

**3-D faces are split along the shorter diagonal.** A bilinear face is not planar. Each minor boundary face is split
into two triangles along its shorter diagonal (c1 c3 on a tie), `wall_cutter`'s rule, so neighbouring faces agree
and the boundary is a closed triangulated surface. Areas and volumes are those of this polyhedral surface, not of the
trilinear cells, and differ from the trilinear ones at second order in the warp. `face_vector` is
½ (c3 − c1) × (c4 − c2), which equals the sum of the two triangles of either split, so `face_covered` matches it to
round-off.

**Checks and refusals.** The boundary must be valid before anything is cut. Only `OVERSET_PATCH_OK` (and, in 3-D,
`OVERSET_PATCH_UNTRUSTED`) computes the cut:

| Status | Meaning |
|---|---|
| `OVERSET_PATCH_BAD_INPUT` | a minor block or face out of range, or major grid edges that are not strictly increasing |
| `OVERSET_PATCH_NO_BOUNDARY` | no boundary face of non-zero length (area) |
| `OVERSET_PATCH_OPEN` | the faces do not close: a 2-D chain without its end, a 3-D edge without exactly one reverse |
| `OVERSET_PATCH_OVERLAP` | two boundary elements cross or touch (exact tests). Two overlapping minor regions need a minor-minor hierarchy: the cut counts regions by parity, so a lens covered twice would come out as fluid |
| `OVERSET_PATCH_NESTING` | a loop's orientation (a shell's signed volume) does not match its nesting depth, or is degenerate |
| `OVERSET_PATCH_UNTRUSTED` | 3-D only, with everything computed: some cut cells holding smallest faces have no trusted pieces (`pieces` = −1). Their smallest faces are listed per cell (`uc_cell`, `uc_sf_start`, `sf_piece` = 0); the solver takes such a cell as one piece, since its totals are right, or refines the grid. None occurred on the validation spheres: 0 of 25 688 cut cells at 96³ |

Several minor blocks are allowed when their loops (shells) are disjoint or nested. CPU only.

### Degenerate pieces

A minor boundary that lies along a grid line (plane) is not cut through any cell. By the cut cells' low-side rule,
the cells on its solid side hold it: each becomes a cut cell with a **degenerate** piece of zero area (volume), whose
wall and aperture lie on the same edge. FOSSIL flags such pieces and leaves them in place: `degenerate(p)` is true,
and `merge_into(p)` names the piece across the coincident edge (the linked piece of largest shared length or area that
is not degenerate itself). The geometry stays exact and every closure holds as it is. Folding is the solver's step,
because a Wang update divides by the piece area: it moves the smallest faces of p to `merge_into(p)` and drops the
shared edge. `merge_into` is 0 for a point-like piece left where a boundary passes through a cell corner. Such a piece
has no fluid edge interval, so no aperture, no smallest face and no flux. Smallest faces of zero length are dropped.
Those of ulp length, where a crossing is computed twice at a grid node, are kept, because the closures need their
area vectors.

![a square ring on a background of unit cells, its outer boundary along the grid lines x = 3, x = 9, y = 3, y = 9: the 25 cells inside it are holes (teal), the boundary is gold, and the 11 cells just inside its top and right sides hold degenerate pieces (ember dots), each joined by a gold segment to the cell across the coincident edge that takes its fluxes](../../examples/images/ref-overset-degenerate.svg){.fossil-figure}

The square [3, 9]² lies on the grid lines. Cells are (x_i, x_i+1], so the left and bottom sides belong to the fluid
cells on their outer side, which only gain a wall on an edge. The top and right sides belong to the 11 cells inside,
which hold degenerate pieces, each with its `merge_into` across the side (`fig_ref_overset_patch.f90`). The validation
folds 59 of them on a 60² grid and conserves mass to 2.8·10⁻¹⁶.

## Intersection weights

`overlap_weights_2d(block, cells, x, y, w)` and `overlap_weights_3d(block, cells, x, y, z, w)` return, for every
receptor cell `cells(:, r)` of a curvilinear block, its intersections with the cells of a rectilinear donor grid, in
CSR form (`start`): the donor cell `donor(:, k)`, `volume(k)` = |V_c ∩ V_c′|, `weight(k)` = volume / receptor volume,
and `centroid(:, k)`, the centroid of V_c ∩ V_c′, which Picard's gradient term needs.

- **The receptor is a polygon (polyhedron).** In 2-D it is the polygon of its 4 vertices, a collapsed edge giving a
  triangle. In 3-D it is the 12 triangles of its 6 faces, each split along its shorter diagonal, with ties broken by
  the least global vertex index so that two receptors split a shared face alike: the receptors of a block tile the
  region they cover, with no gap and no overlap. Exactly collapsed triangles (a pole) are dropped.
- **How it is cut.** The receptor is cut as a hole in a solid box against the donor sub-grid around its bounding box
  (`curve%cut_cells`, `cut_cells_3d`), in coordinates shifted to the sub-grid's low corner. Green's integrals about the
  origin lose relative accuracy as |x| / h grows (1·10⁻¹² at 100 unshifted); the shift keeps them at round-off. The
  fluid of donor cell c′ is then V_c ∩ V_c′, and its area (volume) and centroid are the cut's cell totals.
- **Sums.** `receptor_volume` is the compensated sum of the receptor's intersections, so its weights sum to 1 to
  2 eps. `covered` is exactly 1 when the receptor lies in the closed donor box, and below 1 otherwise. Every volume
  is ≥ 0 and every weight lies in [0, 1].
- **Zero-area donors.** The donor list is the exact cut's: every donor cell with a fluid piece, in (j, i) or (k, j, i)
  order. Under the cut's symbolic perturbation a receptor side on an interior grid line, or through a grid node, can
  leave a donor whose intersection has zero area. It is listed with a volume of zero (or round-off), and its centroid
  is the donor cell's centre.
- **Flags.** Each receptor gets a `flag`. `OVERLAP_WEIGHTS_DEGENERATE` means fewer than 3 distinct vertices (4
  triangles in 3-D), a non-positive area (volume) in the block's handedness, or sides that cross or fold (exact). With
  `OVERLAP_WEIGHTS_OUTSIDE` no intersection has positive volume. `OVERLAP_WEIGHTS_OPEN` means a donor cell of the cut
  did not close, which is not expected. Flagged receptors have no donors. `status` is `OVERLAP_WEIGHTS_BAD_INPUT`
  (nothing computed) for donor edges that are not finite and strictly increasing, or for a receptor out of range.
  Serial, CPU only.

![a receptor cell of a coarse O-grid (16 by 4 cells, r from 1 to 2), between 45 and 67.5 degrees on its outer ring, over a Cartesian donor grid of spacing 0.15: the receptor outlined in gold, its 18 intersections with the donor cells filled from cream to sepia by weight and labelled with the weights above 0.01, their centroids as teal dots, and one more donor cell, whose corner the receptor's side passes through, marked by an ember dot at its centre as the zero-area donor](../../examples/images/ref-overset-weights.svg){.fossil-figure}

Receptor 3 of the outer ring has 19 donors. One is a donor of zero area: the receptor's side along the 45-degree
diagonal passes exactly through that cell's corner. The figure program clips every receptor against its donor cells
on its own and finds the library's areas to 3.1·10⁻¹⁶ (`fig_ref_overset_weights.f90`).

B transfers states. For a linear field with its exact gradient the transfer is exact, and per transfer it conserves
the donor's reconstructed mass over the receptor to round-off. Each time step still overwrites the receptor values
instead of exchanging fluxes, so the composite does not conserve. The validation below separates the two effects.

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

## Example: the conservative interface

The O-grid of the [Example](#example) above, its outer boundary named `overlap`, over the same background: A by
`overset_patch_2d`, then B for its outer ring of cells by `overlap_weights_2d`.

<<< @/examples/snippets/ref_overset_conservative.f90

<<< @/examples/output/ref_overset_conservative.ansi{ansi}

- **A.** The 32 outer faces cut the background into 500 holes and 108 cut cells. Each of the 1100 pieces closes to
  3.3·10⁻¹⁵ h, every face is covered by its smallest faces to 1.4·10⁻¹⁶, and the fluid area equals the box less the
  32-gon to 1.7·10⁻¹⁴. Minor face 1 (cell (1, 12), angle 0 to 11.25 degrees) is split into 3 smallest faces in the
  background column i = 34. Each `sf_area` points to −x, into the O-grid.
- **B.** The 32 receptors of the outer ring have 214 intersections, all within the donor box. Their weights sum to 1
  within 2.2·10⁻¹⁶, and a linear field with Picard's gradient term is exact to 7.2·10⁻¹⁵ at the receptor centroid.
  Receptor 1 has 6 donors in two columns of the background.

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
! conservative interface: a rectilinear major / donor grid of cell edges x, y (z)
call overset_patch_2d(blocks, minor, x, y, patch, boundary_name='overlap', faces=, status=)    ! A; 3-D: overset_patch_3d
call overlap_weights_2d(block, cells, x, y, w, status=)                                         ! B; 3-D: overlap_weights_3d
```

The mask `masks(b)%iota(v)` is 1 (field), 0 (hole) or −k (receptor from block k) at the linear vertex index v (i
fastest); `masks(b)%donor(v)` is the chosen candidate in `donors(b)`, `orphan(v)` an `ORPHAN_*` cause, `fringe(v)` a
`FRINGE_*` level. A cell-centred solver passes the blocks of its cell centres.
The conservative interface returns `overset_patch_2d_t` / `overset_patch_3d_t` (holes, the cut, the smallest faces by
piece and by minor face, `degenerate`, `merge_into`; status `OVERSET_PATCH_*`) and `overlap_weights_t` (CSR receptor
→ donor cell, volume, weight, centroid; `receptor_volume`, `covered`, `flag`; status `OVERLAP_WEIGHTS_*`).

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
- `fossil_test_overset_patch_2d/3d`, `fossil_test_overlap_weights_2d/3d` (conservative interface): the cut equal to
  `cut_cells` / `cut_cells_3d` of the same boundary bit for bit; every major piece closes to 1e-15 h (2-D) and
  1e-12 h² (3-D); per piece Σ `sf_area` = −`piece_wet_normal`; covered faces equal their face vectors (1e-14, 1e-13);
  the pieces' total = box − minor region; every degenerate piece has a non-degenerate merge target; each refusal
  produced on purpose. Weights sum to 1 within 2 eps; a constant field exact (4 eps); a linear field exact with
  Picard's gradient term (1e-13 of its variation); the weighted centroids give the receptor's; receptors 100 h from
  the origin as accurate as at it.
- `fossil_val_overset_conservative`: a passive scalar, first-order upwind with forward Euler, advected through the
  cylinder O-grid of the step-0 probe (64 m × 8 m cells, r in [1, 2]) in the background [−3, 3]² of 60 m × 60 m cells,
  in the potential flow past the cylinder, coupled three ways: A the patched interface (small pieces merged by
  `merge_targets`, degenerate ones by `merge_into`), I bilinear interpolation (holes at centre radius < 1.4, receptors
  below 1.6), B the same assembly with Picard's weights. Each number below measures a different thing:

  | Measure, m = 1 / m = 2 | A | B | I |
  |---|---|---|---|
  | composite defect, max over the steps of \|M + out − M₀\| / M₀ | 1.4·10⁻¹⁶ / 2.8·10⁻¹⁶ | 2.1·10⁻² / 1.0·10⁻² | 2.2·10⁻² / 1.2·10⁻² |
  | the same on the exact partition | (the same) | 3.1·10⁻² / 1.0·10⁻² | 3.0·10⁻² / 1.2·10⁻² |
  | one transfer: receptor mass vs donor-side mass | no transfer | 1.1·10⁻¹⁵ / 3.7·10⁻¹⁵ | 1.9·10⁻⁴ / 5.0·10⁻⁵ |
  | one transfer: receptor mass vs exact integral | no transfer | 1.2·10⁻⁴ / 3.0·10⁻⁵ | 7.3·10⁻⁵ / 2.0·10⁻⁵ |

  - **Composite defect.** M is the composite mass and "out" the mass that left through the box. For A, M covers its
    pieces and the minor cells, an exact partition of the box. For I and B, M covers the major cells whose centre
    lies at r ≥ 2 plus every minor cell, a bookkeeping that counts a cell straddling r = 2 whole or not at all.
  - **Exact partition.** Each major cell is weighted by its fluid fraction outside the minor polygon (the cut of
    `overset_patch_2d`), plus every minor cell. This removes the bookkeeping and leaves the coupling alone. It barely
    changes the defect: the coupling is the cause.
  - **One transfer.** A smooth field is transferred once, with no time step, in both directions: the major receptors
    from the minor cells, then the minor outer row from the major cells. The first row compares the mass the transfer
    puts in the receptors with the donor-side mass over the same region, the donor's linear reconstruction integrated
    over the intersections. B conserves that mass per transfer to round-off; I, which interpolates point values,
    misses it at second order. The second row compares with the exact integral (a 7-point rule per triangle): both
    transfers are second order.

  B's composite defect is thus no better than interpolation's here. Its weights conserve each transfer, but every time
  step overwrites the receptors instead of exchanging fluxes, so the field cells next to them gain or lose mass that
  nothing balances. This is Picard et al.'s own statement that the method transfers states, not fluxes. Their
  0.1-0.6 % (a stationary shock, second-order scheme) and these 1-3 % (a blob advected through a fringe, first order)
  are different tests and do not compare. On a square ring whose boundary lies on grid lines (60², 59 degenerate
  pieces folded) A conserves to 2.8·10⁻¹⁶. 3-D, a sphere shell (r in [1, 2], 32 × 16 × 8 cells, poles collapsed) in
  40³: `overset_patch_3d` in 0.08 s, pieces closed to 2.5·10⁻¹⁴ h², faces covered to 3.6·10⁻¹⁵, the pieces' volume
  = box − polyhedron to 4.7·10⁻¹⁶, Σ `sf_area` = −`piece_wet_normal` bit for bit on 54 194 of 55 792 pieces (the
  rest to 4.7·10⁻¹⁶ h²); the scalar's defect 1.2·10⁻¹⁶ for A and 4.1·10⁻² for I. B in 3-D, with the whole shell
  (24 × 12 × 4) as receptors over 24³: 18 430 intersections in 0.44 s, weights within 2.2·10⁻¹⁶ of 1, a linear
  field exact to 3.7·10⁻¹⁶, and per donor cell the receptors' intersections sum to |V_c′ ∩ shell| (from
  `overset_patch_3d` of the shell's two faces) within 8.9·10⁻¹⁵ h³: the receptors split their shared faces alike and
  agree with the boundary's split.
- `fossil_val_overset_conservative_proto` (the step-0 design probe, kept): the variants the library did not take (B
  by the complement cut, and unshifted), timing per receptor, the patched interface at 48³ and 96³ (0 untrusted cut
  cells of 25 688 at 96³), and the overlapping-minor case that `OVERSET_PATCH_OVERLAP` now refuses.

## Proved and measured

| Claim | Kind |
|---|---|
| hole flags on closed walls | exact (winding number), CPU = device bit for bit |
| hole map with crossing parity, X-rays on their rays | exact; X-rays between rays measured (errors within the raster spacing) |
| donor candidates | WP2's location; equal to a brute force; CPU = device bit for bit |
| classification fixed point | Jacobi iteration, independent of the block order (tested); circular pairs detected |
| interpolation order through the overlap | measured: the rate equals q (fixed overlap); the order loss of a stencil one point short needs a discretization, left to the solver |
| patched interface (A): pieces close, faces covered, total volume | exact geometry (the cut cells); identities to round-off, tested in 2-D and 3-D |
| A conserves the composite | measured to round-off with a passive scalar (2-D, 3-D, degenerate pieces folded) |
| intersection weights (B) sum to 1, linear fields exact | exact geometry (the cut cells); tested to 2 eps and 1e-13 |
| B conserves per transfer, not over the composite | measured: round-off per transfer; 1-3 % composite drift, as interpolation (Picard et al.'s statement) |

## Left out

- Conservative interface: curvilinear major or donor grids (polyhedron clipping of trilinear cells), the device port,
  moving minor grids (the swept geometry of the minor boundary), a hierarchy of overlapping minor grids, and the
  solver itself (fluxes, reconstruction, the fold of degenerate pieces). Other routes the wiki page records, not built:
  Berger's conservative interpolation, Chandar & Sitaraman's global mass constraint, Devlin et al.'s virtual faces.
- Overlap optimization beyond the measures, and parallel (MPI) assembly: the caller's.
- Stencils crossing into a third block: reported low-order, not followed.
