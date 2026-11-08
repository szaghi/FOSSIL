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

![a layer of 40 by 40 boxes through the middle of the dragon, the part of the dragon above the layer cut away: the 179 boxes the surface crosses in cream, the 222 boxes within half a box side of it in ember around them, the 46 boxes inside the body in teal; the 1153 outside boxes are left blank within the layer's outline](../../examples/images/ref-amr-tagging-codes.svg){.fossil-figure}

One `box_classify_many` batch over the layer gives every code at once, with no field computed. The near ring is
`radius` = half a box side wide; the inside boxes are the ones whose centre the fast winding number puts in the body.
`box_intersects_surface` asks the same cut question for a single box:

<<< @/examples/output/figref-amr-tagging-codes.ansi{ansi}

## Example

The NACA 0012 section (the wing cut at z = 0.5, chord from −0.25 to 0.75) in 2-D: a quadtree built level by level
from `box_classify_many`, then tag bands on one stretched block.
[Tutorial chapter 9](/manual/tutorial/09-amr-tagging) classifies one level of blocks and tags a uniform block around
the dragon; the cookbook [tags the cells near a cube](/manual/cookbook#tag-the-cells-near-the-wall). Here the
classification drives a whole hierarchy, and the block is stretched. The 3-D calls are the same, with `lo(3, n)`,
`hi(3, n)` and a `z` edge array.

<<< @/examples/snippets/ref_amr_tagging-tree.f90

- Each level is one `box_classify_many` batch over all its boxes. A box is refined when it is `BOX_CUT`, or
  `BOX_NEAR` within `radius` = the diagonal of its children: the criterion of Matsushita & Aoki (2021), decided from
  the tree before any field exists.
- The rest become leaves. Their code says what they are: outside or inside (one fast-winding-number query at the
  centre), or, at the last level, cut and near.

<<< @/examples/snippets/ref_amr_tagging-tags.f90

- `tag_block` takes the cell edges `x(49)`, `y(25)`; the y spacing starts at 0.006 at the chord line and grows by
  1.1. The tags count cells, not lengths, so a band of 3 is thin where the cells are thin.
- `feature_band=` and `feature_tags=` add the second field: the Chebyshev cell distance to the cells holding a corner
  of the curve, a vertex where it turns by more than `feature_angle` (default 30°). The only one here is the trailing
  edge.

<<< @/examples/output/ref_amr_tagging.ansi{ansi}

- Six levels end with 508 leaves; the uniform grid of the finest size (side 0.0234) has 4096 cells. No box is inside
  before level 6: none of side 0.047 or more fits in an airfoil 0.12 thick.
- The leaves of the last level are 96 cut, 126 near and 48 inside boxes: the near ring is wider than the cut ring
  because the radius is the children's diagonal, 0.0166, against a box side of 0.0234.
- On the block, 112 cells are cut; the counts per ring grow outward with the stretched cells, and 286 cells beyond
  the band get the cap, +4. Inside, 112 cells get −4: the airfoil is more than three of the thin cells thick there.
- 4 cells hold the trailing edge, and 36 are within 2 cells of it.

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

![a block of 48 by 24 cells around the NACA 0012 section, stretched away from the chord line, coloured by its tags with band 3: the cut cells in cream along the airfoil, the inside cells in teal, the outside rings in darker ember up to the cap +4; the cells within 2 cells of the trailing-edge corner outlined in teal](../../examples/images/ref-amr-tagging-tags.svg){.fossil-figure}

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

![a 24 by 18 block around a square ring with a square island in its hole, every cell labelled with its tag: 0 (cream, unlabelled) on the cells the loops cross, then 1, 2 and the cap 3 outward in ember and -1, -2 inward in teal; the uncut cells form 4 components (outside the ring, the ring's solid, the hole, the island's solid) and the first cell of each, in column-major order, is outlined in gold: the one cell whose centre gets a point query](../../examples/images/ref-amr-tagging-bands.svg){.fossil-figure}

On this block the breadth-first search gives every count, and four point queries give every sign: the 180 cells
outside the ring, the 112 of its solid, the 32 of the hole and the 4 of the island's solid are one component each.
The distances are capped at band + 1 = 3, which only the outer columns reach.

<<< @/examples/output/figref-amr-tagging-bands.ansi{ansi}

## Refinement criteria from the literature

| Criterion | With FOSSIL |
|---|---|
| Refine every cut cell to the finest level (Natarajan) | `box_classify_many` on the cells: refine `BOX_CUT` |
| Refine where the smallest \|φ\| in a cell is below the diagonal of its children (Matsushita & Aoki 2021) | `radius` = child diagonal: refine `BOX_CUT` and `BOX_NEAR` |
| Keep coarse–fine boundaries away from cut cells (Giuliani et al. 2022) | `tag_block`: refine where \|tag\| ≤ the required buffer |
| Tag bands for cut cells, neighbours, extension, reinitialisation (Hoppe et al. 2022) | `tag_block(..., band=4)` |

![the quadtree of the example around the NACA 0012 section, six levels from a root box of side 1.5: the outside leaves outlined, growing away from the airfoil; at the finest level the cut leaves in cream along the airfoil, the near leaves in ember around them and the inside leaves in teal](../../examples/images/ref-amr-tagging-quadtree.svg){.fossil-figure}

The quadtree of the [example](#example) uses the second criterion: every box is refined while it is cut or near
within its children's diagonal, so the finest level wraps the airfoil in a ring of near cells.

![the section of the dragon through the middle of its box (two loops, gold) under a quadtree built one level per frame: each frame classifies the boxes of one level in one batch, cut boxes in cream and near boxes in ember are split into the next level, inside boxes (teal) and outside boxes (outlined) stay as leaves; level 6 ends with 544 boxes, 144 cut, 158 near and 80 inside, and 709 leaves in all](../../examples/images/anim-ref-amr-tagging-levels.svg){.fossil-figure}

The same criterion level by level on the section of the dragon: no box is inside before level 5, when the boxes
become smaller than the body's cross-section; the 6 levels end with 709 leaves where the uniform grid of the finest
size has 4096 cells.

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

![a cylinder of 48 sides in a block of 16 cubed cells: the 88 cells touched by its two rims drawn as ember cubes, rings at the top and the bottom; the side, whose facets turn by 7.5 degrees, below the 30-degree threshold, carries no feature edge and no tag](../../examples/images/ref-amr-tagging-features.svg){.fossil-figure}

680 cells meet the cylinder, but only the 88 along the two rims meet a feature edge: the smooth side is resolved by
the surface tags alone, and the rims get the extra refinement. `box_feature_many` on the same 4096 cells marks the
same 88:

<<< @/examples/output/figref-amr-tagging-features.ansi{ansi}

## Refining thin walls and gaps

A wall or a gap thinner than the cell breaks an immersed boundary even where the surface is tagged. After
`surface%thickness` and `surface%build_thin`, `tag_block(..., thin_factor=, thin_band=, thin_tags=)` adds the Chebyshev
cell distance to the nearest cell that meets a facet whose local size (min of thickness and gap) is below
`thin_factor` times the cell's largest side (issue #62). Each block compares with its own cell size, so refinement
stops where the cells are fine enough. `box_thin_many(_device)` answers the same question for any list of boxes. See
[thickness and gaps](/guide/advanced/thickness#amr-criterion).

![a box 1 by 1 by 0.6 and, 0.06 off its +x face, a plate 0.04 thick (ember), seen from the front: the cells of side 0.1 that meet the plate or the box face across the gap outlined in gold, two columns of cells over the whole height and depth of the gap](../../examples/images/ref-amr-tagging-thin.svg){.fossil-figure}

The local size of a facet is the smaller of its thickness and its gap: 0.043 on the plate (the ray cone's bias over
0.04) and 0.064 on the box face across the gap. With `thin_factor=1` the cells of side 0.1 that meet those facets are
tagged, 154 of them; on a block of side 0.03 none is, because every local size is above the cell:

<<< @/examples/output/figref-amr-tagging-thin.ansi{ansi}

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
