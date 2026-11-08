---
title: Mesh decimation
---

# Mesh decimation (QEM)

## What it does

Reduces the number of triangles in a surface while preserving its
overall shape. Iteratively collapses the edge that introduces the
smallest geometric error, as measured by the **quadric error metric**
(Garland & Heckbert 1997).

Complementary to [isotropic remeshing](/guide/advanced/isotropic-remesh):
decimation reduces facet *count* and is shape-preserving (silhouette
priority); isotropic remeshing equalizes edge *length* without
changing the count. Reach for decimation when "this 2 M-facet STL
needs to fit in 100 k facets for cut-cell tagging." Reach for remesh
when "the tessellation is irregular and tripping AABB query balance."

A rounded box (half-sides 0.8, 0.45, 0.3, corner radius 0.15) by marching
cubes, 6460 facets, decimated to a sixteenth: 402 facets, the volume
1.971 to the printed digits. The flat faces are left with a few large
triangles; the rounded edges and corners keep the small ones:

![left, the rounded box by marching cubes, a regular mesh of 6460 small triangles; right, the same box after decimate to 402 facets: the flat faces covered by a few long triangles, the rounded edges by strips along them, the corners by fans of small triangles](../../examples/images/ref-decimate-box.svg){.fossil-figure}

## Pipeline

1. **Per-vertex 4×4 quadrics**. For each vertex, accumulate
   $Q_v = \sum_{p \in \text{incident planes}} K_p$ where $K_p = p\,p^T$
   for plane equation $p = (a, b, c, d)$ with normal $(a,b,c)$ and
   offset $d$. The quadric measures squared distance to all incident
   facet planes.
2. **Per-edge optimal collapse target**. For each edge $(u, v)$, the
   minimizer of $v^T (Q_u + Q_v) v$ is the optimal collapse position.
   3×3 Cramer's-rule solve; midpoint fallback if singular.
3. **Min-heap**, keyed on collapse cost. Decrease/increase-key
   supported (~80 LOC inline) so neighbour edges can have their costs
   re-evaluated after each collapse without rebuilding the heap.
4. **Iterative collapse** with **three safety checks** per candidate:
   - **Link condition** (Dey, Edelsbrunner, Guha, Nekhayev 1999): the
     edge $(a, b)$ collapses only if $\mathrm{Lk}(a) \cap \mathrm{Lk}(b) =
     \mathrm{Lk}(ab)$: the only common neighbours of $a$ and $b$ are the
     apexes of the edge's two facets, and no tetrahedron $(a, c, d)$,
     $(b, c, d)$ closes over it. On an open surface the links include a
     dummy vertex joined to the boundary, so two boundary vertices merge
     only along a boundary edge. A manifold stays a manifold of the same
     topology: no non-manifold edge, no pinched shell.
   - **Normal-flip rejection**: if any incident triangle's normal would
     flip, or the triangle would collapse to zero area, the collapse is
     skipped. A zero-area triangle puts two vertices at one place, which
     the rebuilt vertex pool merges into a pinch: on flat faces, where
     every collapse costs 0, this gave `cube-01.stl` decimated to 36
     facets 2 non-manifold edges and 2 shells with `DEC_STATUS_OK`.
   - **Degenerate-facet rejection**: if a surviving facet would get two
     equal vertices, skipped.
   When all remaining edges fail one of these checks, the algorithm
   returns `DEC_STATUS_NO_PROGRESS` with the smallest valid surface it
   reached.

Stages 1 and 2 computed on the rounded box above (by the figure program,
from the facet planes and the vertex ids of the pool): each facet
coloured by its cheapest edge, the cost $\bar m^T (Q_u + Q_v) \bar m$ at
the edge midpoint. On 5604 of the 6460 facets it is below 1e-10: on the
flat faces, and along the cylindrical edges, where a collapse along the
axis moves no vertex off its planes. Only the corners, curved both ways,
cost up to 7.8e-4. The heap takes the cream edges first:

![the rounded box with every facet coloured by its cheapest edge collapse on a log scale: cream on the flat faces and along the cylindrical edges, sepia only on the eight spherical corners](../../examples/images/ref-decimate-quadric.svg){.fossil-figure}

Decimated to targets halved one after the other (each a fresh call on the
input), the box loses its flat-face triangles first and the corners last;
the volume moves by 2e-14 down to 806 facets, by -1.1 % at 100:

![an animation of the rounded box decimated to 6460, 3230, 1614, 806, 402, 200 and 100 facets, each frame with the volume change: the flat faces empty first into a few large triangles, the cylindrical edges thin into long strips, the corners coarsen last](../../examples/images/anim-ref-decimate-target.svg){.fossil-figure}

## API

```fortran
call surface%decimate(target_facets, status)
   class(surface_stl_object), intent(inout)         :: self
   integer(I4P),              intent(in)            :: target_facets
   integer(I4P),              intent(out), optional :: status   ! DEC_STATUS_*
```

- **`target_facets`** is the **upper bound** on the output facet
  count. The algorithm collapses edges until it reaches this count or
  cannot safely collapse any more.
- **`status`**:
  - `DEC_STATUS_OK` — reached `target_facets`: the count is at or just
    below it (a collapse removes two facets, so an odd target can end
    one below: the dragon's 6588 / 4 = 1647 ends at 1646).
  - `DEC_STATUS_NO_PROGRESS` — couldn't reach the target; surface is
    at the smallest count the safety checks allowed.
  - `DEC_STATUS_BAD_INPUT` — the input surface has no facets, or
    `target_facets <= 0`; the surface is left untouched.

After return, `self`'s AABB tree, vertex pool, connectivity, and
pseudo-normals are all rebuilt — no further user action needed before
calling `distance`, `winding_number`, `boolean`, or any other query.

## Example

The cookbook's [Decimate a mesh](/manual/cookbook#decimate-a-mesh)
reduces the dragon once. Here the question is how far each level moves
the surface, on a body whose shape is known: a unit sphere by
[marching cubes](/guide/advanced/marching-cubes) (41³ samples of
[-1.5, 1.5]³, cleaned by `sanitize`):

<<< @/examples/snippets/ref_decimate-sphere.f90

Three levels, each a quarter of the previous, measured against the
input by the certified [Hausdorff distance](/guide/advanced/surface-deviation):

<<< @/examples/snippets/ref_decimate-levels.f90

<<< @/examples/output/ref_decimate.ansi{ansi}

- **Count.** Each level reaches its target with `DEC_STATUS_OK`. An odd
  target ends one below (103 asked, 102 given): a collapse removes two
  facets.
- **Distance.** The Hausdorff interval brackets the true distance to
  the input: 0.0033 to 0.0068 for a quarter of the facets, 0.0585 to
  0.0619 for 1/64. It grows about 4× per level, as the sagitta of a
  chord, $e^2/8$ for an edge $e$ on a unit sphere, does when the facet
  area quadruples.
- **Volume.** It drifts down level after level: 0.2 %, 1.1 % and
  4.5 % below the input's (see
  [Known limitations](#known-limitations)).
- **Topology.** Every level stays watertight and manifold: the link
  condition keeps the input's topology on any manifold input.

## Visual reference

The example's bunny (69 451 facets) is too heavy for a figure; the
same call on the dragon (`dragon.stl`), 6588 facets as loaded:

![the dragon as loaded, 6588 facets, flat-shaded with faint facet edges](../../examples/images/ref-decimate-dragon-input.svg){.fossil-figure}

After `decimate(target_facets=n0 / 4)` (1646 facets, from the same
camera):

![the dragon decimated to a quarter of its facets, 1646: the same silhouette with coarser facets on the body, horns, claws and tail kept](../../examples/images/ref-decimate-dragon-output.svg){.fossil-figure}

The silhouette is preserved; per-facet area roughly quadruples (each
remaining facet covers roughly 4× the area of an original). The
thin features (horns, claws, the tail's crest) are preserved because
those regions have high quadric error per edge, so the heap keeps them
out of the collapse queue until the budget forces it; the smooth body
takes most of the collapses.

## Known limitations

These are the **documented gaps** from §1.3's implementation issue,
deliberately surfaced here so the user knows what they're getting:

- **No tie-breaking heuristic.** On perfectly flat inputs (a flat
  cube face, an axis-aligned rectangle) every edge has cost ≈ 0 and
  the algorithm picks arbitrarily among the valid collapses.
  `cube-01.stl` (48 facets) decimated to 36 facets stays watertight,
  manifold and one shell with the unit volume; asked for 4 facets, it
  ends at a tetrahedron of volume 0.158.
- **No cost ceiling.** The loop collapses the cheapest edge until the
  target is met or no edge passes the checks, whatever the cost: a
  target too low for the shape is met by distorting it. On the bracket
  of [tutorial 12](/manual/tutorial/12-geometry) built with 128-sided
  prisms (1868 facets), the patch boundaries (three rims of 128
  vertices, never merged) leave few cheap collapses: decimating to 1000
  or 600 facets stops with `DEC_STATUS_NO_PROGRESS` at 1582 facets, a
  valid closed surface with 37 % of the volume lost: once the cheap
  collapses are spent, the base box is sheared along its edges.
  Compare the volume before and after.
- **No volume-preservation term.** Collapses pull vertices inward,
  so volume drifts down monotonically with iteration. The
  marching-cubes sphere of the [Example](#example) (6632 facets)
  decimated to 1/4, 1/16 and 1/64 of its facets loses 0.2 %, 1.1 % and
  4.5 % of its volume. Lindstrom-Turk 1998 adds a volume-
  preservation Lagrange multiplier; deferred.
- **No boundary preservation.** Edges on a mesh boundary (open
  shells) are treated like interior edges; the boundary moves
  during collapse (the link condition keeps it one loop per hole).
  Irrelevant for closed solids (the typical CFD case).

## See also

- [isotropic remeshing](/guide/advanced/isotropic-remesh) — the
  complementary primitive for tessellation regularity.
- [`surface%get_facets_number`](/guide/api/surface-stl-object) — read
  the count before/after decimation.
- [Constants](/guide/api/constants) — `DEC_STATUS_*`.

## References

- Garland & Heckbert, *Surface Simplification Using Quadric Error
  Metrics*, SIGGRAPH 1997. The reference algorithm.
- Dey, Edelsbrunner, Guha & Nekhayev, *Topology preserving edge
  contraction*, Publ. Inst. Math. (Beograd) 66, 1999. The link
  condition.
- Lindstrom & Turk, *Fast and memory efficient polygonal
  simplification*, IEEE Visualization 1998. Adds the volume-
  preservation term that addresses the documented sphere drift.
- Hoppe, *New quadric metric for simplifying meshes with appearance
  attributes*, IEEE Visualization 1999. Generalization to colored
  meshes; not relevant to FOSSIL's per-vertex-coordinate-only API.
