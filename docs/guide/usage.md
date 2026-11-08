---
title: Usage
---

# Usage

All examples use `use fossil` which re-exports `surface_stl_object` and `facet_object`, plus the constants for sign algorithms (`SIGN_PSEUDO_NORMAL`, `SIGN_RAY_INTERSECTIONS`, `SIGN_SOLID_ANGLE`, `SIGN_FAST_WINDING_NUMBER`), tree kinds (`AABB_TREE_SAH_BVH`, `AABB_TREE_OCTREE`), and status codes (`STATUS_OK`, `STATUS_INVALID_INPUT`, …).

Numeric kinds (`I4P`, `R8P`) come from [PENF](https://github.com/szaghi/PENF); 3D vectors and unit versors (`ex_R8P`, `ey_R8P`, `ez_R8P`) come from [VecFor](https://github.com/szaghi/VecFor).

Every block on this page is a region of a program in `docs/examples/src/` (`ref_usage_*.f90`), compiled and run when
the site is built; the outputs shown are those runs. The bodies are the test fixtures `dragon.stl` (6588 facets, not
watertight), `cube-01.stl` (the unit cube) and `naca0012-ascii.stl`.

## Loading an STL file

### Auto-detect format

Pass `guess_format=.true.` to let FOSSIL determine whether the file is ASCII or binary automatically.

<<< @/examples/snippets/ref_usage_io-load.f90

`load_from_file` runs `analyze` internally — bounding box, connectivity, volume, centroid, AABB tree, and pseudo-normals are all populated before it returns. You typically call `sanitize` next to repair common STL defects.

![the dragon as loaded, inside its bounding box with the corners bmin and bmax labelled, its centroid marked, and in red the edges without a single neighbouring facet (fcon_edge = 0) at the jaws, the horns and the tail: on the dragon as loaded these are its 61 non-manifold edges, not open boundaries (see Validity predicates)](../examples/images/ref-usage-analyze.svg){.fossil-figure}

### Explicit format

<<< @/examples/snippets/ref_usage_io-explicit.f90

`is_ascii` defaults to `.true.` when `guess_format` is absent or `.false.`.

### Load with on-the-fly clipping

Only the facets that overlap the box are loaded, without a separate clip pass; a facet crossing the box is kept whole,
so the loaded surface can reach out of it:

<<< @/examples/snippets/ref_usage_io-clipped.f90

### Rejecting bad input

If any vertex coordinate is `NaN` or `Inf`, `load_from_file` refuses the file via the optional `status` argument rather than letting the bad value propagate:

<<< @/examples/snippets/ref_usage_io-reject.f90

The run of the loading blocks (and of the save block below, which writes `dragon-out-binary.stl`): the clipped load
keeps 3811 of the 6588 facets, and its x range ends at 1.2261, past the box face x = 0, because the facets crossing
that face are kept whole.

<<< @/examples/output/ref-usage-io.ansi{ansi}

---

## Printing statistics

<<< @/examples/snippets/ref_usage_analysis-statistics.f90

The output for the cube fixture, post-sanitize:

<<< @/examples/output/ref-usage-statistics.ansi{ansi}

The last line is the depth of the octree; the default SAH BVH has no levels and reports `AABB_AUTO_REFINEMENT` (-1).

The same cube drawn with what the statistics describe: the box of the extents, the centroid, and the 48 facet normals,
all pointing out, which is why the divergence-theorem volume is +1.0000:

![the cube fixture after sanitize drawn translucent inside its bounding box, the extents [0, 1] labelled, its centroid (0.5, 0.5, 0.5) marked, and the normal of each of its 48 facets drawn from the facet centroid, every one pointing outwards](../examples/images/ref-usage-statistics.svg){.fossil-figure}

::: tip
Volume is positive for a closed body with outward-pointing normals (the standard divergence-theorem convention; `sanitize` produces this state). If you see a negative volume after sanitize, the mesh is probably non-watertight — check the non-manifold-edge and disconnected-edge counts.
:::

---

## Saving an STL file

<<< @/examples/snippets/ref_usage_io-save.f90

`is_ascii` defaults to `.true.`: pass `is_ascii=.false.` for a binary file. A surface with several patches is written as
one `solid` block per patch in ASCII.

An STL file has no shared vertices: each facet is written on its own, its normal and its three vertices. The cube
fixture saved this way holds 144 vertex records for its 26 distinct vertices; `load_from_file` welds them back
(`analyze` builds the vertex pool), and the file loaded back has the same 48 facets and volume 1.0000:

![the cube fixture exploded, every facet facing the camera shrunk towards its centroid and pushed out along its normal, each with its own three vertices as gold points; one facet's normal is drawn and labelled, one of its vertices labelled](../examples/images/ref-usage-saving.svg){.fossil-figure}

The [Named patches](./advanced/patches.md#reading-and-writing) page shows a multi-patch file read back.

---

## Sanitize (the full repair pipeline)

`sanitize` is the orchestrator. It runs every repair pass in the right order and prints one stderr warning per defect class detected.

<<< @/examples/snippets/ref_usage_repair-sanitize.f90

On the dragon it removes 62 duplicate facets and reports what it cannot repair: facets with disconnected edges (the
holes) and 54 non-manifold edges remain.

<<< @/examples/output/ref-usage-sanitize.ansi{ansi}

The passes that run inside `sanitize`:

1. `remove_degenerate_facets` — drop zero-area / sliver triangles before their NaN normals contaminate downstream queries.
2. `connect_nearby_vertices` (if disconnected edges were detected) — union-find vertex deduplication on a tolerance.
3. `analyze` — rebuild metrix, connectivity, vertex occurrences.
4. `remove_duplicate_facets` — drop literal duplicates, including reversed-orientation ones.
5. `analyze` — rebuild after duplicate removal (only if any duplicates were dropped).
6. `sanitize_normals` — BFS-propagate winding consistency, then global flip so volume > 0.
7. `fill_holes` — only with `sanitize(fill_holes=.true.)`: close the boundary loops with minimum-weight patches
   (`refine_holes` and `fair_holes` refine and fair the patches; see [Hole filling](./advanced/hole-filling.md)).

Each individual pass is also callable directly if you want fine-grained control. Called one by one, in this order, on
the unit cube with one defect of each kind (a degenerate facet, a facet shifted by 1e-4, a duplicate facet, two facets
wound inwards and a missing facet), each pass removes its own defect:

![an animation of the unit cube with five defects labelled (a degenerate facet as three gold points on the top front edge, a crack along a shifted facet in red, a duplicate facet in brown, two facets wound inwards in ember, a hole), then one frame per pass: remove_degenerate_facets removes 1 facet, connect_nearby_vertices closes the crack, remove_duplicate_facets removes 1 and leaves the 3 boundary edges of the hole, sanitize_normals flips 2 facets, fill_holes closes the hole: watertight, volume 1.0000](../examples/images/anim-ref-usage-sanitize.svg){.fossil-figure}

On the dragon the passes show a less obvious effect: as loaded it has no open boundary at all, only 61 non-manifold
edges, several of them made by duplicate facets; removing the 62 duplicates turns those edges into 128 boundary edges
(see the next section).

---

## Validity predicates

After `sanitize`, three named predicates summarise the mesh state:

<<< @/examples/snippets/ref_usage_analysis-validity.f90

| Predicate | True iff |
|---|---|
| `is_watertight` | Every edge has exactly 2 incident facets |
| `is_manifold`   | No non-manifold edges (boundary edges allowed) |
| `is_volume`     | watertight AND volume > 0 AND finite centroid |

The block run on the dragon as loaded, then on the cube after `sanitize`:

<<< @/examples/output/ref-usage-validity.ansi{ansi}

What the predicates see on the dragon, from `feature_edges` (boundary edges have one facet, non-manifold edges three or
more): as loaded, 0 boundary and 61 non-manifold edges; after `sanitize`, 128 boundary and 54 non-manifold edges. All
three predicates stay false in both states:

![an animation of the dragon drawn translucent in two frames: as loaded, with its 61 non-manifold edges in gold and no boundary edge; after sanitize, with 128 boundary edges in red, where duplicate facets were removed, and 54 non-manifold edges in gold](../examples/images/anim-ref-usage-validity.svg){.fossil-figure}

---

::: warning Call analyze after transforming
`translate`, `rotate`, `mirror` and `resize` move the vertices only. The
centroid, volume, bounding box, AABB tree and pseudo-normals stay as they were
until `call surface%analyze`, so getters and distance queries return stale
values until then. One `analyze` after the last transform of a batch is enough.
See [After a transform](/guide/api/surface-stl-object#after-a-transform).
:::

The blocks below start each operation from a copy of the dragon as loaded (`surface = dragon`), declared as:

<<< @/examples/snippets/ref_usage_manipulation-declarations.f90

## Translate

Move the surface by a vectorial delta or by individual axis components (`delta` and `x`/`y`/`z` together are
ambiguous: `status = STATUS_AMBIGUOUS_ARGS`, or an error stop without `status`):

<<< @/examples/snippets/ref_usage_manipulation-translate.f90

![the dragon as loaded as a translucent ghost and, solid, the dragon translated by (4, 3.5, 0.5), the segment between the two centroids marked](../examples/images/ref-usage-translate.svg){.fossil-figure}

---

## Rotate

Rotate around an arbitrary axis by an angle in radians; the pivot is the origin:

<<< @/examples/snippets/ref_usage_manipulation-rotate.f90

Alternatively, provide a pre-built rotation matrix:

<<< @/examples/snippets/ref_usage_manipulation-rotate_matrix.f90

Both forms take an optional `center`, the pivot of the rotation:

<<< @/examples/snippets/ref_usage_manipulation-rotate_center.f90

![the dragon as loaded as a ghost above the x axis, and the dragon after a quarter turn about the x axis through the origin, lying below the axis](../examples/images/ref-usage-rotate.svg){.fossil-figure}

---

## Mirror

Mirror with respect to a plane through the origin, defined by its normal (normalised inside):

<<< @/examples/snippets/ref_usage_manipulation-mirror.f90

<<< @/examples/snippets/ref_usage_manipulation-mirror_oblique.f90

Or by a pre-built mirror matrix:

<<< @/examples/snippets/ref_usage_manipulation-mirror_matrix.f90

![the dragon seen from above: the original as a ghost lying along x, and its reflection through the plane x + y = 0, drawn as a line, lying along y](../examples/images/ref-usage-mirror.svg){.fossil-figure}

The facet winding is reversed together with the reflection, so an outward
surface stays outward (the same holds for `resize` with an odd number of
negative factors): the volumes of the mirrored dragons in the run below stay
+465.929. As for every transform, call `surface%analyze` before queries.

---

## Resize (scale)

Scale all facets by a uniform or per-axis factor (`factor` and `x`/`y`/`z` together are ambiguous, as for
`translate`). The default pivot is the origin:

<<< @/examples/snippets/ref_usage_manipulation-resize.f90

<<< @/examples/snippets/ref_usage_manipulation-resize_axes.f90

Pass `respect_centroid=.true.` to scale about the surface centroid instead (the centroid of the last `analyze`):

<<< @/examples/snippets/ref_usage_manipulation-resize_centroid.f90

![the dragon as loaded as a ghost and, solid, the dragon scaled by 2 about the origin, marked below: twice as large and moved away from the origin](../examples/images/ref-usage-resize.svg){.fossil-figure}

---

## Clip

Keep only the facets that overlap a given axis-aligned bounding box (a facet crossing the box is kept whole). The removed facets can be captured as a separate surface:

<<< @/examples/snippets/ref_usage_manipulation-clip.f90

![the clip box around the head half of the dragon: the 3811 kept facets solid, reaching out of the box where they cross its face, and the 2777 facets of the remainder as a ghost](../examples/images/ref-usage-clip.svg){.fossil-figure}

---

## Merge

Combine two surfaces into one. The result is stored in the first surface; here the two parts are the clip and its
remainder, saved and loaded back:

<<< @/examples/snippets/ref_usage_manipulation-merge.f90

![the merged surface: the facets of dragon_part_1 in stone and those of dragon_part_2 in ember, together the whole dragon again](../examples/images/ref-usage-merge.svg){.fossil-figure}

The run of every block of the manipulation sections: the box and the volume that `analyze` computes after each
operation. Rotations, mirrors and translations keep the volume 465.929; the uniform factor 2 multiplies it by 8
(3727.429), the factors (0.5, 1, 1.2) by 0.6 (279.557). The merge gives back the 6588 facets and the box of the
dragon.

<<< @/examples/output/ref-usage-manipulation.ansi{ansi}

---

## Distance queries

The first blocks query the unit cube `cube-01.stl`:

<<< @/examples/snippets/ref_usage_distance-load.f90

### Unsigned distance (squared, default)

<<< @/examples/snippets/ref_usage_distance-unsigned.f90

### Euclidean (non-squared) distance

<<< @/examples/snippets/ref_usage_distance-euclidean.f90

### Signed distance

Negative inside, positive outside. Default sign algorithm is `SIGN_PSEUDO_NORMAL` (Bærentzen–Aanæs); pass `sign_algorithm=` to override.

<<< @/examples/snippets/ref_usage_distance-signed.f90

### Point-in-polyhedron test

<<< @/examples/snippets/ref_usage_distance-inside.f90

### Closest point and normal

`compute_distance` returns the distance together with the closest facet, the local edge (0 if the closest point is not on
an edge) and the local vertex (0 if it is not on a vertex) of that facet, the closest point and the pseudo-normal there:

<<< @/examples/snippets/ref_usage_distance-closest.f90

![probe points around the dragon, each joined to its closest point on the surface, with the pseudo-normal drawn as a short segment from the closest point; where the closest point lies on an edge or a vertex the normal is not the facet normal](../examples/images/ref-usage-closest.svg){.fossil-figure}

### Batched distances

`distance_many` is a loop of `distance` over a point array, OpenMP-parallel, with the same options for every point;
`distances` is allocated by the caller, and a size mismatch returns `STATUS_INVALID_INPUT`:

<<< @/examples/snippets/ref_usage_distance-batched.f90

### Choosing the AABB tree kind

The default `AABB_TREE_SAH_BVH` is faster on every measured mesh (see [AABB tree](./api/aabb-tree-object.md) and
[tutorial 5](/manual/tutorial/05-batched)). The legacy octree is still selectable:

<<< @/examples/snippets/ref_usage_distance-tree_kind.f90

Both kinds produce bit-exact distance results — the difference is purely query speed. The run of the distance blocks:

<<< @/examples/output/ref-usage-distance.ansi{ansi}

---

## Connect nearby vertices (repair)

`sanitize` calls this internally when needed; you only invoke it directly if you want to repair a surface without
running the full pipeline. The surface here is the cube with its first facet shifted by 1e-4 along x, which
disconnects its three edges from the neighbours:

<<< @/examples/snippets/ref_usage_repair-crack.f90

<<< @/examples/snippets/ref_usage_repair-connect.f90

Each group of nearby vertices is replaced by its centroid, so the crack closes (no disconnected edges after) but the
merged vertices sit between the old positions: the volume moves from 0.9999958 to 0.9999948.

<<< @/examples/output/ref-usage-connect.ansi{ansi}

Zoomed on one corner of the shifted facet, at the point (0, 0.5, 0), 4e-4 across:

![an animation of two frames zoomed on a corner of the cube, 4e-4 across: before, the shifted facet's vertex and the neighbours' vertex as two gold points 1e-4 apart, with the disconnected edges in red; after connect_nearby_vertices, one vertex at the centroid of the group and no red edge](../examples/images/anim-ref-usage-connect.svg){.fossil-figure}

---

## Full workflow example

<<< @/examples/snippets/ref_usage_workflow.f90

<<< @/examples/output/ref-usage-workflow.ansi{ansi}

![the dragon after the workflow: its centroid at the origin, marked, with the x, y and z axes drawn from it, inside its bounding box, whose longest side (along x) is 1 and whose corner bmin is (-0.4695, -0.2440, -0.2455)](../examples/images/ref-usage-workflow.svg){.fossil-figure}
