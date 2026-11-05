---
title: Feature edges
---

# Feature edges

## What it does

Finds the edges an immersed-boundary grid must resolve (issue #61, roadmap #58):

- **sharp edges**, where the two facets meet at more than an angle (30° by default);
- **boundary edges**, with one facet: the rims of holes and open surfaces;
- **non-manifold edges**, with three or more facets;
- the **corners**, where feature edges meet or turn sharply;
- the **polylines** that chain the edges from corner to corner.

IB methods smear sharp corners and edges unless the grid is refined there. This is the output that drives
feature refinement in CFD meshers (OpenFOAM's `surfaceFeatureExtract` feeds snappyHexMesh with it).

![a machined part, a block with a raised boss and a round hole drilled through it, its 155 feature edges in ember and its 16 corners in gold](../../examples/images/ref-feature-edges.svg){.fossil-figure}

The part is built with [`boolean`](/guide/advanced/booleans): a block, united with a boss, minus a 24-sided hole. Its
12 block edges, 8 boss edges and 4 edges where the boss meets the block are sharp, which gives 16 corners and 24
polylines. The hole's two rims are closed polylines without corners: the hole's facets meet at 15°, below the
threshold, so the hole wall itself has no features.

## Method

- **Edges.** The half-edges of all facets are sorted by (lower vertex, higher vertex, facet) with a stable merge sort,
  so the result is deterministic. Each unique edge records its number of facets and its first two.
- **Sharp.** A two-facet edge is sharp when the unit normals of its facets make an angle above `angle`:
  $n_1 \cdot n_2 < \cos(\text{angle})$. The normals come from the vertex order, so the facets must be consistently
  oriented, as [`sanitize`](/guide/usage) leaves them. A facet of zero area has no normal. Its edges are not
  classified as sharp, and are counted in `degenerate`.
- **Corners.** A feature vertex is a corner when it has a number of feature edges other than 2, or exactly 2 that
  turn by more than `corner_angle` (30° by default).
- **Polylines.** From every corner, each unused feature edge is followed through the other feature vertices until the
  next corner. The edges left form closed loops, and each loop starts and ends at the same vertex.
- **Shared with the remesher.** `isotropic_remesh(preserve_features=.true.)` locks its features with the same sharp
  predicate and the same 30° threshold (`is_sharp_dihedral`, `feature_cos` in `fossil_feature_edges`), so the two
  cannot drift apart. Its results are unchanged, bit for bit.

## API

```fortran
call surface%feature_edges(edges, kinds, corners, polyline_start, polyline_vertices, angle, corner_angle, &
                           include_boundary, include_nonmanifold, degenerate, status)
   integer(I4P), allocatable, intent(out), optional :: edges(:,:)           ! [2, ne] vertex-pool ids
   integer(I4P), allocatable, intent(out), optional :: kinds(:)             ! FEATURE_EDGE_SHARP / _BOUNDARY / _NONMANIFOLD
   integer(I4P), allocatable, intent(out), optional :: corners(:)           ! vertex-pool ids
   integer(I4P), allocatable, intent(out), optional :: polyline_start(:)    ! [np + 1], CSR offsets
   integer(I4P), allocatable, intent(out), optional :: polyline_vertices(:) ! polyline p: start(p) .. start(p+1)-1
   real(R8P),                 intent(in),  optional :: angle, corner_angle  ! degrees, default 30
   logical,                   intent(in),  optional :: include_boundary     ! default .true.
   logical,                   intent(in),  optional :: include_nonmanifold  ! default .true.
   integer(I4P),              intent(out), optional :: degenerate           ! edges skipped for a zero-area facet
   integer(I4P),              intent(out), optional :: status               ! STATUS_OK, STATUS_INVALID_INPUT
```

- Vertex ids index the vertex pool: `surface%get_vertex_pool()`, `pool%coord(id)`.
- Edges are listed in increasing (low id, high id) order.
- `status = STATUS_INVALID_INPUT` on an empty surface, or an angle outside [0°, 180°].

## Example

<<< @/examples/snippets/ck_feature_edges.f90

<<< @/examples/output/ck_feature_edges.ansi{ansi}

The dragon as loaded has 5104 sharp edges at 30°: its scan is coarse, and neighbouring facets often meet at steep
angles. It has no boundary edge: its 62 duplicated facets (removed by `sanitize`, tutorial chapter 2) put three or more
facets on 61 edges, which count as non-manifold. Raising the angle to 60° and excluding the other kinds keeps 2508
creases.

## Validation

`src/tests/fossil_test_feature_edges.f90`:

| Case | Expected |
|---|---|
| `cube-01` (48 facets) | 24 sharp edges, all on the cube's 12 edges, of total length 12; the 8 cube corners; 12 polylines |
| A 32-sided cylinder, caps fanned | the two cap rims (64 sharp edges), no corners, two closed polylines of 32 edges |
| `cube-01` without one facet | 3 boundary edges; none with `include_boundary=.false.` |
| Three triangles sharing one edge | 1 non-manifold edge, 6 boundary edges |
| `dragon.stl`, sanitized | the same count after a translation; non-increasing as the angle goes 10°, 20°, 30°, 45°, 60° |

`fossil_test_isotropic_remesh` gives the same output as before the shared predicate (compared line by line against a
build of the previous commit).

## See also

- [Isotropic remeshing](/guide/advanced/isotropic-remesh): `preserve_features` locks these edges.
- [AMR tagging](/guide/advanced/amr-tagging): the grid the features are meant for.
