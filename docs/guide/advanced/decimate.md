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
5. **Patch rules** (see [Patches](#patches)): an edge that touches a
   patch rim has a constrained target. The heap holds the free
   optimum, a lower bound of the constrained cost, and an edge is
   re-costed when it reaches the top; on a surface with one patch
   nothing changes, bit for bit.
6. **Optional certified tolerance** `max_error` (see
   [A certified tolerance](#a-certified-tolerance)).

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
call surface%decimate(target_facets, max_error, rim_collapses, status)
   class(surface_stl_object), intent(inout)         :: self
   integer(I4P),              intent(in)            :: target_facets
   real(R8P),                 intent(in),  optional :: max_error     ! certified tolerance, a distance
   logical,                   intent(in),  optional :: rim_collapses ! default .true.
   integer(I4P),              intent(out), optional :: status        ! DEC_STATUS_*
```

- **`target_facets`** is the **upper bound** on the output facet
  count. The algorithm collapses edges until it reaches this count or
  cannot safely collapse any more.
- **`status`**:
  - `DEC_STATUS_OK` — reached `target_facets`: the count is at or just
    below it (a collapse removes two facets, so an odd target can end
    one below: the dragon's 6588 / 4 = 1647 ends at 1646).
  - `DEC_STATUS_NO_PROGRESS` — couldn't reach the target; surface is
    at the smallest count the safety checks allowed, or the cheapest
    collapse costs more than `max_error**2`.
- **`max_error`** (optional, a distance): every point of the result
  lies within `max_error` of the input, certified; a collapse that
  would break it is skipped, and the loop stops when the cheapest
  admissible collapse costs more than `max_error**2` (see
  [A certified tolerance](#a-certified-tolerance)). Absent: neither.
- **`rim_collapses`** (optional, default `.true.`): patch rims coarsen
  along themselves; `.false.` locks them (see [Patches](#patches)).
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

## Patches

Every facet keeps its [patch](/guide/advanced/patches). A vertex shared
by facets of two patches is on a *rim*, and the target of an edge that
touches a rim is constrained (issues #66, #81):

- **One end on a rim**: the collapse goes to the rim vertex, which
  stays where it is.
- **Both ends on a rim**: only along a rim edge (its two facets in
  different patches), to the cheaper of its two endpoints, so a rim
  coarsens as a curve whose vertices are a subset of the input rim's,
  bit for bit. A *corner*, a rim vertex with other than two rim edges
  (three or more patches meet there), never moves: an edge from a
  corner goes to the corner, and two corners never merge. Any other
  edge between two rim vertices is not collapsed, nor is the last
  facet of a patch.
- **Constraint planes**: each rim edge adds to the quadrics of its two
  ends the plane through it perpendicular to the mean normal of its
  facets (the boundary constraint of Garland & Heckbert), so a rim
  lying in a plane is not straightened for free.

`rim_collapses=.false.` locks the rims instead (the rule before #81):
two rim vertices never merge, and the rim vertices are the input ones
exactly.

The bracket of [tutorial 12](/manual/tutorial/12-geometry) built with
128-sided prisms (1868 facets, three rims of 128 vertices), each frame a
fresh call on the input. To 1000, 600 and 300 facets it reaches the
target with `DEC_STATUS_OK` and moves the volume by -6.0e-6, 2.2e-5
and 2.4e-4; the rims thin out along the prisms. With the rims locked,
asked for 600, it stops at 1580 facets with `DEC_STATUS_NO_PROGRESS`
and 42 % of the volume gone: once the cheap collapses are spent, the
loop takes costly ones (the box sheared along its edges) until none
passes the checks. With `max_error = 1e-2` it stops at 294 facets:

![an animation of the bracket coloured by patch, base sepia, boss teal, bore ember: the input with 1868 facets, then decimated to 1000, 600 and 300 facets with the rims of the boss and the bore coarsened into fewer sides, then with max_error 1e-2 at 294 facets, then with the rims locked: 1580 facets and the box sheared](../../examples/images/anim-ref-decimate-bracket.svg){.fossil-figure}

The validation `fossil_val_decimate_bracket` asserts it: OK, valid,
the patches kept and every output rim vertex an input one, |dV/V| ≤
1e-4; the regression `fossil_test_decimate_rims` takes `cube-01` with
a patch per face to the 12-facet cube (8 corners kept, before:
`NO_PROGRESS` at 36) and a 64-gon prism with three patches, every
vertex on a rim, to 64 facets (before: no collapse at all).

## A certified tolerance

With `max_error`, two rules act on the loop (roadmap
[#92](https://github.com/szaghi/FOSSIL/issues/92) WP5 made the first one
a certificate):

- **Every collapse is certified.** A collapse is applied only if its
  target and every new facet of its one-ring lie within `max_error` of
  the input surface. Each new facet is checked by a branch and bound
  against the input's AABB tree: a triangle is settled when the least of
  three upper bounds of the distance over it is at most `max_error`, else
  it is split in four at its edge midpoints. The bounds: the largest
  corner distance plus the triangle's radius (the distance is
  1-Lipschitz); the largest corner distance to one input facet F, among
  the corners' closest ones (d(x, input) ≤ d(x, F), and d(x, F) is convex
  in x, so its largest value over the triangle is at a corner); and the
  same per part over a split of the triangle along the edge between two
  of those facets, or into the wedges of the fan of an input vertex they
  share, each part bounded by its own facet. A midpoint farther than `max_error` refuses the collapse; so
  does a facet still undecided after 40 levels or 4096 triangles (never
  reached in any measured run). A refused collapse is **skipped and the
  loop goes on**: its edge stays alive and is queued again when a
  neighbouring collapse moves one of its ends. The guarantee holds to
  round-off: h(result, input) ≤ `max_error`.
- **The plane cost still stops the loop.** The cost of a collapse is the
  sum of the squared distances of the new vertex to the planes it
  carries (unit normals, weight 1): the facet planes merged into its two
  ends and the constraint planes. The loop stops with
  `DEC_STATUS_NO_PROGRESS` when the cheapest admissible collapse has
  $\sqrt{\text{cost}} > $ `max_error`. This bounds no distance of the
  surface, but it is the only guard of the other side, the input close
  to the result: cutting off a feature merges its planes into a vertex,
  and that costs. Without the stop, with the certificate alone, the
  dragon at 1e-2 of its diagonal goes to 1020 facets instead of 2960 and
  the input lies up to 3 `max_error` from the result.

Before the certificate the stop was the only rule, and a collapse could
put a new edge across a concave patch: `dragon.stl` at `max_error` = 3e-4
of its diagonal ended 1.70 to 1.74 `max_error` from the input (6534
facets; in the collapse 6546 → 6544 the new vertex is 0.04 `max_error`
from the input, the midpoint of a new edge 1.15), and is now 0.83 to 0.87
`max_error` from it (6538 facets, 2 collapses refused). `bunny.stl` at
3e-3 of its diagonal: 1.109 to 1.135 `max_error` before, 0.525 to 0.575
now (10 905 facets, 1 refusal). On the bracket, the wing and the
marching-cubes spheres the plane cost was already conservative and no
collapse is refused: their counts and distances are the ones before. The cost of the certificate is +15 to
+60 % of a decimation of the dragon (0.09 to 0.22 s) and +6 % on the
bunny; without `max_error` none of it runs.

With `max_error`, the open edges of the input (one facet) get
constraint planes too (through the edge, perpendicular to its facet):
without them a vertex of an open edge slides along its plane for free
(the NACA 0012 wing, asked for 1e-3 of its diagonal, had a tip vertex
moved 1.65 along the span; now its Hausdorff distance is 1.03e-3 for a
max_error of 1.42e-3). Without `max_error`, open edges stay
unconstrained, as before.

It certifies one side, the result close to the input: see
[Known limitations](#known-limitations). The certified
[Hausdorff distance](/guide/advanced/surface-deviation), both sides,
measured below `max_error` on closed bodies:

| Body | `max_error` | Facets | Certified Hausdorff ≤ |
| --- | --- | --- | --- |
| bracket (1868 facets) | 3e-3, 1e-2 | 458, 294 | 1.8e-3, 6.3e-3 |
| 64-gon prism, 3 patches (252) | 2e-2 | 102 | 1.26e-2 |
| marching-cubes sphere, 32³ (2252) | 1e-3, 3e-3, 1e-2, 3e-2 | 1804, 1604, 1260, 664 | 6.6e-4, 2.2e-3, 4.4e-3, 1.7e-2 |

The bracket and its decimations with `max_error = 3e-3` and `1e-2` (target 4: as few facets as the
tolerance allows), the Hausdorff distances certified to a quarter of `max_error`:

![three brackets side by side, coloured by patch: the input with 1868 facets, then the decimations with max_error 3e-3 and 1e-2, the flat faces of the box covered by a few large triangles, the prisms keeping their sides](../../examples/images/ref-decimate-ceiling.svg){.fossil-figure}

The cookbook's [Decimate to a tolerance](/manual/cookbook#decimate-to-a-tolerance)
decimates the wing to 0.1 % of its diagonal and checks the result.

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
- **No tolerance by default.** Without `max_error` the loop
  collapses the cheapest edge until the target is met or no edge
  passes the checks, whatever the cost: a target too low for the shape
  is met by distorting it (the bracket with locked rims above, 42 % of
  its volume gone). Pass `max_error`, or compare the volume before and
  after.
- **`max_error` certifies one side.** Every point of the result is
  within `max_error` of the input; the input points a collapse leaves
  behind are not certified, only held by the plane-cost stop (a small
  spike flattened costs little against its own planes). On the closed
  bodies of [A certified tolerance](#a-certified-tolerance) both sides
  stay below it; on `bunny.stl` at 1e-3 of its diagonal (2.5e-4) the
  input lies up to 3.1e-4 to 3.7e-4 from the result, on `dragon.stl` at
  2.67e-2 up to 0.11 (the result within 2.67e-2, certified). A fold, a
  fin, a thin part of the input can go: the first collapse of
  `dragon.stl` at 3e-4 of its diagonal removes a fold whose facets 399
  and 2777 lie back to back, and leaves the input 5.0 `max_error` from
  the result. Check with `hausdorff` when that side matters.
- **No volume-preservation term.** Collapses pull vertices inward,
  so volume drifts down monotonically with iteration. The
  marching-cubes sphere of the [Example](#example) (6632 facets)
  decimated to 1/4, 1/16 and 1/64 of its facets loses 0.2 %, 1.1 % and
  4.5 % of its volume. Lindstrom-Turk 1998 adds a volume-
  preservation Lagrange multiplier; deferred.
- **No boundary preservation by default.** Edges on a mesh boundary
  (open shells) are treated like interior edges; the boundary moves
  during collapse (the link condition keeps it one loop per hole).
  With `max_error` they get constraint planes (see
  [A certified tolerance](#a-certified-tolerance)). Irrelevant for closed solids
  (the typical CFD case).

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
- Klein, Liebich & Straßer, *Mesh reduction with error control*, IEEE
  Visualization 1996; Cohen et al., *Simplification envelopes*,
  SIGGRAPH 1996. Two-sided error control by tracking which input facets
  each output facet replaces: not implemented (FOSSIL certifies the
  result side only). Cited from general knowledge, not checked against
  the papers.
- Lindstrom & Turk, *Fast and memory efficient polygonal
  simplification*, IEEE Visualization 1998. Adds the volume-
  preservation term that addresses the documented sphere drift.
- Hoppe, *New quadric metric for simplifying meshes with appearance
  attributes*, IEEE Visualization 1999. Generalization to colored
  meshes; not relevant to FOSSIL's per-vertex-coordinate-only API.
