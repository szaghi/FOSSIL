---
title: Isotropic remeshing
---

# Isotropic remeshing

## What it does

Equalizes **edge length** across a triangulated surface while preserving
its geometry. CAD-exported STL routinely mixes 0.01 mm and 10 mm
triangles on the same body — perfectly fine for visualization, **fatal**
for AABB-tree query balance and for cut-cell quality near refinement
bands. Isotropic remeshing rewrites the tessellation driving every edge
towards the band `(4L/5, 4L/3)` about a target length `L`, producing a
mesh whose facets have roughly the same shape and size (the band is a
goal, not a bound: on the dragon below, after 3 iterations at `L = 0.46`,
90 % of the edges lie in 0.30–0.57).

Complementary to [decimation (§1.3)](https://github.com/szaghi/FOSSIL/issues/18):
decimation reduces facet *count* while preserving the silhouette;
isotropic remeshing equalizes edge *length* without preferring either
direction. Use decimation to make a million-facet model tractable; use
isotropic remeshing to make any model's tessellation regular.

A torus by marching cubes on cells 2.75 times flatter than wide, so its
triangles are stretched (edges 0.001 to 0.137, 5280 facets), remeshed
with `L = 0.08`, 3 iterations (edges 0.037 to 0.118, 3124 facets):

![left, the torus by marching cubes, its facets coloured by their longest edge over L, a patchwork of dark long triangles and light small ones; right, the torus after isotropic_remesh, an even light colour: every facet's longest edge close to L](../../examples/images/ref-isotropic-remesh-torus.svg){.fossil-figure}

## Pipeline

Each outer iteration runs four passes (Botsch & Kobbelt 2004):

1. **Split** edges longer than `4L/3`. New vertex at the midpoint;
   incident triangles re-triangulated.
2. **Collapse** edges shorter than `4L/5`, into the midpoint (or into
   the endpoint that is locked, see feature preservation). A collapse
   is rejected if it would break the topology: a non-manifold edge, an
   interior edge between two boundary vertices, a tetrahedron, or two
   endpoints that share a neighbour other than the apexes of the edge's
   two facets (the link condition); or if it would create an edge
   longer than `4L/3` (which the next split would undo), flip a
   triangle's normal or duplicate a vertex. A vertex takes part in one
   collapse per pass.
3. **Flip** interior edges to drive vertex valence toward 6 (the
   regular-triangulation ideal). Rejects flips that would create a
   duplicate edge or a degenerate triangle.
4. **Tangential relaxation**: move each vertex toward the area-weighted
   centroid of its one-ring, keeping only the part of the move in the
   vertex's tangent plane (`(I - n nᵀ)` applied to it, n the vertex
   normal), then **project** back onto the original surface via the AABB
   tree's closest-point query. Before
   [issue #91](https://github.com/szaghi/FOSSIL/issues/91) the whole move
   was applied, and on a coarse mesh it pulled the vertices inward.

**Error bound.** Collapse, flip and relaxation refuse a change when a
facet it creates (sampled at its centroid and edge midpoints) lies
farther than `REMESH_ERROR_TOL · L` from the input surface, with
`REMESH_ERROR_TOL = 0.07` (issue #91, after Dunyach et al. 2013). A flip
or a move whose facets are already that far is not held back. Without
the bound a target near the size of the body ate it: the collapses merged
the corners of the cube away, the flips cut across its edges, the
relaxation slid the vertices off them, and `cube-01.stl` at `L = 1` ended
as a flat tetrahedron of volume 4e-4. The constant comes from a sweep of
0.03 to 0.3 on the cube, the sphere, the torus and the dragon below:

| `REMESH_ERROR_TOL` | cube `L = 1`, volume | cube `L = 0.2`, volume | sphere `L = 0.4`, facets | dragon `L = 1.5`, facets |
|---|---|---|---|---|
| none | 0.137 | 0.931 | 298 | 1382 |
| 0.2 | 0.480 | 0.937 | 298 | 1428 |
| 0.1 | 0.833 | 0.977 | 296 | 1806 |
| **0.07** | **0.949** | **0.986** | **332** | **2280** |
| 0.05 | 0.955 | 0.980 | 412 | 2888 |
| 0.03 | 0.990 (no collapse) | 0.983 | 642 | 4006 |

(features off, 3 iterations, 5 on the cube at `L = 0.2`; the sphere of the
[Example](#example), whose chords at `L = 0.4` already sag 0.02 = 0.05 L.)
Below 0.07 the cube keeps little more and the curved bodies stop
coarsening; above it the cube loses its corners.

Three to five outer iterations is typical. Connectivity is rebuilt from
scratch between passes — chosen after the first attempt's packed-array
incremental design produced uncatchable heap corruption; the
"rebuild every pass" approach trades speed for debuggability.

What the first iteration acts on, on the sphere of the [Example](#example)
(marching cubes on 25³ samples, `L = 0.1`): 1424 edges longer than
`4L/3` to split, 696 shorter than `4L/5` to collapse, 508 of the 1158
vertices with a valence other than 6 for the flips:

![a zoom on the marching-cubes sphere with its edges coloured by the pass that acts on them: long edges to split in ember, short edges to collapse in teal, the others grey; the vertices of valence other than 6 as gold points, about half of them](../../examples/images/ref-isotropic-remesh-passes.svg){.fossil-figure}

The same classes after each of five calls of one iteration. The long
edges go: 1424, 130, 26, 10, 4, 3. The short ones rise first, 696 to
1458, because the splits halve edges just above `4L/3` into edges just
below `4L/5` and a vertex collapses once per pass; then the collapses
take them: 422, 175, 55, 24. The facet count rises from 2312 to 3654
and falls back to 2760; 4113 of the 4140 edges end in the band:

![an animation of the sphere after 0 to 5 iterations of isotropic_remesh at L = 0.1, edges coloured as above: the ember edges disappear in the first two iterations, the teal ones multiply after the first and then fade, the facets settle at 2776](../../examples/images/anim-ref-isotropic-remesh-iterations.svg){.fossil-figure}

**Feature preservation** is optional and off by default. When enabled,
edges with dihedral angle > 30° are flagged as **sharp**, and their
endpoint vertices are **locked**: collapse skips them, flip skips any
edge with a locked vertex, relaxation skips them entirely. The features
are detected again after every split pass, so the midpoints inserted on a
sharp edge are locked like its endpoints (before issue #42 they were not,
and relaxation chamfered every edge: a unit cube remeshed to edge 0.17
lost 10 % of its volume). This is what keeps the 12 edges of a cube sharp
through remeshing, with volume and area unchanged to round-off; without
it, the cube rounds out into a ball over a few iterations.

## API

```fortran
call surface%isotropic_remesh(target_length, iterations, preserve_features, status)
   real(R8P),    intent(in), optional :: target_length     ! ≤ 0 → median input edge
   integer(I4P), intent(in), optional :: iterations        ! default 5
   logical,      intent(in), optional :: preserve_features ! default .false.
   integer(I4P), intent(out), optional :: status
```

- **`target_length`**: desired uniform edge length. Pass `≤ 0` (or
  omit) to use the input's *median* edge length — a remesh at the
  input's scale that evens the distribution (on an input with a long
  tail of long edges, such as the dragon below, the splits outnumber
  the collapses and the facet count grows). Pass a smaller
  value to densify, larger to coarsen (on the sphere of the
  [Example](#example), `L = 0.4` takes 2312 facets to 332 in 3
  iterations; the error bound stops a coarsening the curvature cannot
  carry: `L = 0.8` and `L = 1.2` stop at 296 facets, edges up to 0.60
  and 0.69, volume 3.94 and 3.93 of 4.15).
- **`iterations`**: outer-loop count. Default `5`. Convergence is
  typically visible after 3 passes; further iterations smooth the edge-
  length distribution further and cost a little volume each on curved
  bodies (see Limitations).
- **`preserve_features`**: lock vertices on edges with dihedral > 30°.
  **Essential for mechanical CAD inputs**; should usually be `.true.`
  on anything that isn't an organic shape.

Mutates `self` in place.

## Example

The cookbook's
[Remesh to a uniform edge length](/manual/cookbook#remesh-to-a-uniform-edge-length)
remeshes the cube once. Here, what `preserve_features` changes on it,
then a curved body watched one iteration at a time.

The unit cube (`cube-01.stl`, edges 0.5 and 0.71) remeshed to
`L = 0.2`, without and with feature preservation:

<<< @/examples/snippets/ref_isotropic_remesh-cube.f90

A unit sphere by [marching cubes](/guide/advanced/marching-cubes) on 25³
samples, whose edges range over an order of magnitude, remeshed to
`L = 0.1` by three calls of one iteration each (the state carries over
from call to call), then measured against the input by the certified
[Hausdorff distance](/guide/advanced/surface-deviation):

<<< @/examples/snippets/ref_isotropic_remesh-sphere.f90

<<< @/examples/output/ref_isotropic_remesh.ansi{ansi}

- **Features off.** Nothing locks the vertices of the cube's edges and
  corners, but the error bound refuses the changes that move a facet more
  than 0.07 L = 0.014 off the cube: the edges are chamfered only that
  much. Area 5.809 and volume 0.986 instead of 6 and 1 (4.833 and 0.872
  before [issue #91](https://github.com/szaghi/FOSSIL/issues/91)).
- **Features on.** The 12 edges (dihedral 90°, above 30°) are locked,
  the midpoints inserted on them too: area and volume stay 6 and 1 to
  the printed digits.
- **Sphere.** The splits raise the facet count from 2312 to 3654, then
  the collapses bring it down to 2880; the shortest edge rises from
  0.0173 to 0.058, the longest falls from 0.20 to 0.15. Each iteration
  costs volume: 0.30 %, 0.24 %, 0.15 % (see
  [Known limitations](#known-limitations)). The remeshed surface is
  within 0.0114 of the input.

## Dragon remesh — visual reference

Input dragon (`dragon.stl`, 6588 facets; the 69k-facet bunny is too heavy
for an inline figure), its tessellation mixing long slivers and small
triangles:

![the dragon as loaded, 6588 facets, a mix of long slivers and small triangles: per-facet edge median 0.462, standard deviation 0.389](../../examples/images/ref-remesh-dragon-input.svg){.fossil-figure}

After 3 iterations of `isotropic_remesh` at the default target (the
median input edge):

![the same dragon after 3 iterations of isotropic_remesh at the default target: 11584 facets, edge median 0.399, standard deviation 0.135](../../examples/images/ref-remesh-dragon-output.svg){.fossil-figure}

The silhouette is preserved; the per-facet edge-length distribution
tightens (standard deviation 0.389 → 0.135, 2.9×) and its median stays
near `L` (0.462 → 0.399). Remeshing isn't decimation, but neither does
it keep the facet count: the dragon's edges have a long tail above the
median (5 % are longer than 1.32, 2.9× the median), every edge above
`4L/3` is split, and the count grows (6588 → 11584), while the edges
below `4L/5` collapse where the error bound allows it. The thin features
(horns, claws) keep their short edges: a collapse there would move the
surface more than 0.07 L. The area falls by 4.7 % (742 → 708), the
volume by 0.3 % (before
[issue #91](https://github.com/szaghi/FOSSIL/issues/91): 7730 facets,
standard deviation 0.091, area −14 %, volume −2.1 %).

## Known limitations

- **Volume loss on curved bodies.** The vertices lie on the input, the
  facets between them are chords: on a convex region they lie inside, by
  about `L² / (8 R)` along an edge of a region of radius `R`. The
  marching-cubes sphere of the [Example](#example) loses 0.30 %, 0.24 %
  and 0.15 % of its volume in its first three iterations at `L = 0.1`.
  The relaxation is tangential, so this is the chords' share, not an
  inward pull (the full Laplacian of the releases before
  [issue #91](https://github.com/szaghi/FOSSIL/issues/91) lost 0.30 %,
  0.25 % and 0.15 % on this sphere, so it was the chords there too).
- **Feature preservation uses a hard 30° dihedral cutoff.** Edges below
  the threshold get smoothed even if they were intended as features
  (subtle creases on a fillet). The cutoff isn't user-tunable in the
  current API; if your input has sub-30° features that must survive,
  pre-sharpen them by inserting a tiny fillet before remeshing.
- **`target_length` is a guideline, not a contract.** The split (`>
  4L/3`) and collapse (`< 4L/5`) thresholds let edges in `[0.8L,
  1.33L]` survive untouched. Ratio fixed in the literature; not a
  parameter.
- **Connectivity is rebuilt every pass**, not incrementally maintained.
  Per-pass cost is `O(N log N)` instead of `O(δ_changed)`. Acceptable
  for typical inputs; a 70k-facet remesh runs in ~1 s. For million-
  facet meshes consider running fewer iterations.
- **Volume drift compounds across iterations and calls.** It shrinks per
  iteration (see the first item) but does not stop: verify it with
  `surface%get_volume()` before and after. Do not use the remesh in
  volume-critical pipelines without an explicit post-remesh
  volume-restoration pass.
- **A target near the size of the body.** The error bound keeps the
  body, not the target: `cube-01.stl` at `L = 1` without
  `preserve_features` keeps its 48 facets in one call of 3 iterations,
  area 5.503 and volume 0.949 (the collapses that would merge its corners
  are refused); the sphere above stops at 296 facets for any `L` from 0.8
  up. Each call bounds the error against its own input, so repeated calls
  add up: three calls of one iteration leave the cube at 40 facets, volume
  0.823. Before [issue #91](https://github.com/szaghi/FOSSIL/issues/91)
  the cube collapsed to a flat tetrahedron (volume 4e-4). With
  `preserve_features` the edges and corners are locked: 36 facets, area 6,
  volume 1.
- **The error bound trades coarsening for fidelity.** Where the
  curvature radius is below about `1.8 L` (the sagitta of an edge `L` is
  then above `0.07 L`) the collapses stop short of the target: the dragon
  at its median edge keeps 11584 facets, at `L = 1.5` 2280 (1382 without
  the bound). `REMESH_ERROR_TOL` is a named constant of `fossil_remesh`,
  not an argument.
- **Locked vertices do not collapse into each other.** Edges between two
  locked vertices (on a crease or a patch boundary) keep their length:
  the bracket of [tutorial 12](/manual/tutorial/12-geometry) keeps its
  0.0038 edges at the rims. A dirty input keeps short edges where its
  non-manifold edges block the collapses: the dragon at `L = 1.5`
  goes from 6588 to 2280 facets, with edges of 2e-4 left.
- Before [issue #79](https://github.com/szaghi/FOSSIL/issues/79) the
  collapse pass did not coarsen: its 1-ring check counted the apexes of
  the edge's own facets, so it rejected almost every collapse on a
  closed mesh. `fossil_test_remesh_coarsen` guards the fix.

## See also

- [`decimate`](/guide/api/surface-stl-object) (facet count
  reduction; complementary primitive).
- [`statistics`](/guide/api/surface-stl-object#statistics) — print
  edge-length distribution before and after to verify the remesh
  worked.
- The [§1.7 implementation issue](https://github.com/szaghi/FOSSIL/issues/18)
  has the full multi-commit history and design notes (rectangular vs.
  packed connectivity, the heap-corruption pivot, the relax-bias
  diagnosis).

## References

- Botsch & Kobbelt, *A Remeshing Approach to Multiresolution Modeling*,
  Symposium on Geometry Processing 2004. The four-pass pipeline and the
  tangential relaxation.
- Dunyach, Vanderhaeghe, Barthe & Botsch, *Adaptive Remeshing for
  Real-Time Mesh Deformation*, Eurographics 2013 short papers. Remeshing
  operations bounded by the approximation error.
- Botsch et al., *Polygon Mesh Processing*, AK Peters 2010, Chapter 6.
  Textbook treatment with implementation guidance.
- Alliez, Cohen-Steiner, Devillers, Lévy, & Desbrun, *Anisotropic
  Polygonal Remeshing*, SIGGRAPH 2003. The anisotropic
  generalization (out of FOSSIL scope; for surfaces with directional
  feature preservation requirements).
