---
title: Mesh decimation (§1.3)
---

# §1.3 — Mesh decimation (QEM)

## What it does

Reduces the number of triangles in a surface while preserving its
overall shape. Iteratively collapses the edge that introduces the
smallest geometric error, as measured by the **quadric error metric**
(Garland & Heckbert 1997).

Complementary to [isotropic remeshing (§1.7)](/guide/advanced/isotropic-remesh):
decimation reduces facet *count* and is shape-preserving (silhouette
priority); isotropic remeshing equalizes edge *length* without
changing the count. Reach for decimation when "this 2 M-facet STL
needs to fit in 100 k facets for cut-cell tagging." Reach for remesh
when "the tessellation is irregular and tripping AABB query balance."

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
   - **Normal-flip rejection**: if any incident triangle's normal would
     flip, the collapse is skipped.
   - **Duplicate-vertex rejection**: if the collapse would merge two
     vertices that already exist on the same facet (creating a
     degenerate triangle), skipped.
   - **Non-manifold-edge rejection**: if the collapse would create an
     edge with > 2 incident facets, skipped.
   When all remaining edges fail one of these checks, the algorithm
   returns `DEC_STATUS_NO_PROGRESS` — the smallest count it could
   safely reach.

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
  - `DEC_STATUS_OK` — reached `target_facets` exactly.
  - `DEC_STATUS_NO_PROGRESS` — couldn't reach the target; surface is
    at the smallest count the safety checks allowed.
  - `DEC_STATUS_BAD_INPUT` — `target_facets <= 0`, or the input
    surface has no facets.

After return, `self`'s AABB tree, vertex pool, connectivity, and
pseudo-normals are all rebuilt — no further user action needed before
calling `distance`, `winding_number`, `boolean`, or any other query.

## Example

```fortran
program ex_decimate
use fossil
use penf, only : I4P, R8P
implicit none

type(surface_stl_object) :: bunny
integer(I4P)             :: n0, status

call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
n0 = bunny%get_facets_number()
print '(A,I0)', 'input bunny: n_facets = ', n0

! Decimate to 25% of original.
call bunny%decimate(target_facets=n0 / 4_I4P, status=status)
print '(A,I0,A,I0)', 'after 25%: n_facets = ', bunny%get_facets_number(), &
                     '  status=', status
endprogram ex_decimate
```

## Visual reference

Stanford bunny — input (69 451 facets):

![Stanford bunny — input](/pictures/advanced/decimate-bunny-input.png)

After decimation to 25 % of input (17 362 facets):

![Stanford bunny — decimated to 25%](/pictures/advanced/decimate-bunny-output.png)

The silhouette is preserved; per-facet area roughly quadruples (each
remaining facet covers roughly 4× the area of an original). The
ear-and-chin detail is preserved because those regions have high
quadric error per edge, so the heap keeps them out of the collapse
queue until the budget forces it.

## Known limitations

These are the **documented gaps** from §1.3's implementation issue,
deliberately surfaced here so the user knows what they're getting.
Production-grade quality requires fixing all four:

- **No tie-breaking heuristic.** On perfectly flat inputs (a flat
  cube face, an axis-aligned rectangle) every edge has cost ≈ 0 and
  the algorithm picks arbitrarily. Result: degenerates on
  cube-shaped inputs, where any collapse produces equally-valid (or
  equally-bad) output.
- **No volume-preservation term.** Collapses pull vertices inward,
  so volume drifts down monotonically with iteration. Sphere-via-MC
  decimated to 50 % / 25 % / 10 % loses ~9 % / ~13 % / ~24 % of
  volume respectively. Lindstrom-Turk 1998 adds a volume-
  preservation Lagrange multiplier; deferred.
- **No boundary preservation.** Edges on a mesh boundary (open
  shells) are treated like interior edges; the boundary moves
  during collapse. Irrelevant for closed solids (the typical CFD
  case); breaks open shell decimation.
- **No link condition.** Some edge collapses produce non-manifold
  vertices even when the per-edge safety checks pass. In practice,
  the three existing checks catch ~99 % of pathological cases on
  CFD-relevant inputs.

## See also

- [§1.7 isotropic remeshing](/guide/advanced/isotropic-remesh) — the
  complementary primitive for tessellation regularity.
- [`surface%get_facets_number`](/guide/api/surface-stl-object) — read
  the count before/after decimation.
- [Constants](/guide/api/constants) — `DEC_STATUS_*`.

## References

- Garland & Heckbert, *Surface Simplification Using Quadric Error
  Metrics*, SIGGRAPH 1997. The reference algorithm.
- Lindstrom & Turk, *Fast and memory efficient polygonal
  simplification*, IEEE Visualization 1998. Adds the volume-
  preservation term that addresses the documented sphere drift.
- Hoppe, *New quadric metric for simplifying meshes with appearance
  attributes*, IEEE Visualization 1999. Generalization to colored
  meshes; not relevant to FOSSIL's per-vertex-coordinate-only API.
