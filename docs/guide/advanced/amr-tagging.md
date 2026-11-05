---
title: AMR tagging
---

# AMR tagging

## What it does

An immersed-boundary solver on an octree or a quadtree rebuilds its
geometry after every regrid. Tagging by sampling the level set φ needs the
field over whole blocks before it can decide where to refine. FOSSIL answers
the tagging questions from the BVH alone, before any field exists:

- does this block, or this cell, **meet** the geometry?
- does the geometry come **within r** of it?
- if neither, is it **inside or outside** the solid?
- how many cells is each cell from the nearest cut cell (a **tag band**)?

| Query | 3-D (`surface_stl_object`) | 2-D (`curve_object`) |
|---|---|---|
| Classify boxes: cut, near, inside, outside | `box_classify_many`, `box_classify_many_device` | `box_classify_many`, `box_classify_many_device` |
| Does one box meet the geometry? | `box_intersects_surface` | `box_intersects_curve` |
| Signed tag bands on a block of cells | `tag_block(x, y, z, band, tags)` | `tag_block(x, y, band, tags)` |

## Usage

```fortran
program tag_blocks
use fossil, only : surface_stl_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: body
real(R8P)                :: lo(3, 2), hi(3, 2), x(17)
integer(I4P)             :: codes(2), tags(16, 16, 16), k

call body%load_from_file(file_name='body.stl', guess_format=.true.)

! two blocks: cut / near (within 0.05) / inside / outside
lo(:, 1) = [0._R8P, 0._R8P, 0._R8P] ; hi(:, 1) = [0.5_R8P, 0.5_R8P, 0.5_R8P]
lo(:, 2) = [2._R8P, 2._R8P, 2._R8P] ; hi(:, 2) = [2.5_R8P, 2.5_R8P, 2.5_R8P]
call body%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.05_R8P)

! tag bands on a 16^3 block: 0 cut, -d inside, +d outside, |d| <= 5
do k = 1, 17
   x(k) = (k - 1) / 32._R8P
enddo
call body%tag_block(x=x, y=x, z=x, band=4, tags=tags)
endprogram tag_blocks
```

## Definitions

- **Cut.**
  - **3-D:** the separating-axis test of a box against each triangle
    (Akenine-Möller 2001: the 3 box normals, the triangle normal, the 9
    edge cross products). It runs in floating point on the box inflated
    outward by `box_tolerance`, 10⁻¹² of the larger of the geometry's
    diagonal and coordinate magnitude. So it is **conservative**: it never
    misses a box that touches the surface, and may call a box closer than
    that tolerance cut. For tagging this is the safe direction: a false
    positive costs one refined cell, a false negative loses a cut cell.
  - **2-D:** **exact.** A segment meets a closed box if an endpoint lies
    in it or it shares a point with a box edge, decided with the exact
    `orient2d`. The segments come from the curve's own vertices, so no
    rounding enters.
- **Near** (`BOX_NEAR`, with `radius` > 0): not cut, but the geometry comes
  within the radius.
  - **3-D:** a triangle meets the box inflated by the radius. It never
    misses a triangle within the radius, and may report one up to √3·radius
    away, beyond a box corner.
  - **2-D:** the segment–box distance is at most the radius (the distance in
    floating point, the cut test still exact).
- **Inside / outside.** A box that meets no facet lies on one side of a
  closed surface, so one point query at its centre decides it, with
  `sign_algorithm`.
  - The default is `SIGN_FAST_WINDING_NUMBER` (issue #39): the benchmark
    found it the most robust on dirty geometry.
  - With one query per box or component, its cost does not matter.
  - Pass `SIGN_PSEUDO_NORMAL` for the previous behaviour.
  - The optional `inside(:)` output also gives the centre sign of near
    boxes.
- **Tag band.** On a rectilinear block given its cell edges (stretched
  allowed):
  - 0 on a cut cell;
  - otherwise the Chebyshev cell distance d to the nearest cut cell, capped
    at band + 1, negative inside (−d) and positive outside (+d).

  `band = 4` gives the cut, neighbour (±1), extension (±3) and
  reinitialisation (±4) rings of Hoppe et al. (2022).

## How it works

- **Box queries.** A BVH traversal prunes every node whose box does not meet
  the query box inflated by the reach, and stops at the first cut. The
  kernels (`bvh_box_contact_oac`, `curve_box_contact_oac`) are shared by the
  CPU (OpenMP over the boxes) and the device (OpenACC, after
  `enter_device`). On an octree the 3-D query falls back to brute force.
- **Tag bands.**
  - A block that meets nothing costs one point query: ±(band + 1)
    everywhere.
  - Otherwise each cell gets the cut test, and a breadth-first search over
    the 26 neighbours (8 in 2-D) gives the Chebyshev distances.
  - Each connected component of uncut cells is signed by one point query.
    Two closed cells that share a point and are both uncut lie on the same
    side of a closed surface: a path between their centres through the
    shared point stays inside the two cells, which the surface does not
    meet.

## Refinement criteria from the literature

| Criterion | With FOSSIL |
|---|---|
| Refine every cut cell to the finest level (Natarajan) | `box_classify_many` on the cells: refine `BOX_CUT` |
| Refine where the smallest \|φ\| in a cell is below the diagonal of its children (Matsushita & Aoki 2021) | `radius` = child diagonal: refine `BOX_CUT` and `BOX_NEAR` |
| Keep coarse–fine boundaries away from cut cells (Giuliani et al. 2022) | `tag_block`: refine where \|tag\| ≤ the required buffer |
| Tag bands for cut cells, neighbours, extension, reinitialisation (Hoppe et al. 2022) | `tag_block(..., band=4)` |

## Results

| Case | Result |
|---|---|
| cube-02 (2000 boxes), bunny (400 boxes) | codes = brute force over all facets; all four codes occur |
| NACA 0012, annulus with an island (3000 boxes each) | codes = brute force over all segments |
| uncut boxes, cube-02 and curves | all corners on the side of the centre |
| cube-02 on the octree | codes = BVH codes |
| touching: box on a cube face plane; 2-D box on a side | cut; 10⁻⁶ off (3-D) or one ulp off (2-D) not cut |
| device path vs CPU, 3-D and 2-D, both signs | 0 mismatches (device path run on the host) |
| tag bands: cube-02 uniform 16³ and stretched 12³, bunny 12³, blocks inside and outside | = brute-force reference (sign checked at every cell centre) |
| tag bands: NACA 0012 stretched 48 × 20, 200-gon 32 × 32 | = brute-force reference |

## Refining along feature edges

Sharp edges and corners need more resolution than the smooth wall. `surface%build_features` followed by
`tag_block(..., feature_band=, feature_tags=)` adds a second tag field: the Chebyshev cell distance to the nearest cell
touched by a feature edge (issue #61). `box_feature_many(_device)` answers the same question for any list of boxes. See
[feature edges](/guide/advanced/feature-edges#refinement-along-the-features).

## Limitations

- **GPU not validated yet.** `box_classify_many_device` (3-D and 2-D) runs
  on the host in `tests-gnu`. The GPU benches print
  `box device vs serial mismatches` for the manual gate, which has not run
  yet.
- **Tag bands see one block.** A cut cell in a neighbouring block does not
  count, so bands that must cross block faces need the block padded with
  ghost cells.
- **The sign needs a closed geometry.** The default fast winding number
  tolerates holes and gaps. Even so, the one-query-per-component rule
  assumes the uncut cells of a component are on one side.
- **No device `tag_block`.** The breadth-first search runs on the CPU, once
  per block; the per-cell cut tests can come from `box_classify_many_device`.
