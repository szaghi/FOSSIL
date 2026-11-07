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
its own wall. Which cutter cuts which block is the caller's choice; the default lets a cutter cut every block except
those that contributed its patches.

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
- **wall distance** (Oh et al. 2024, `methods/wall-distance-hole-cutting.md`): a point receives from the block in
  which it is closer to the wall, relative to the block's wall spacing.

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

## Proved and measured

| Claim | Kind |
|---|---|
| hole flags on closed walls | exact (winding number), CPU = device bit for bit |
| hole map with crossing parity, X-rays on their rays | exact; X-rays between rays measured (errors within the raster spacing) |
| donor candidates | WP2's location; equal to a brute force; CPU = device bit for bit |
| classification fixed point | independent of the iteration order (tested); circular pairs detected |
| interpolation order through the overlap | measured: the rate equals q when the overlap suffices, one order lower with a stencil one point short |

## Left out

- Conservative-interface geometry (`methods/conservative-overset-interface.md`): a later WP decides.
- Overlap optimization beyond the measures, and parallel (MPI) assembly: the caller's.
- Stencils crossing into a third block: reported low-order, not followed.
