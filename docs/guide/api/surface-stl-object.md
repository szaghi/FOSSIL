---
title: surface_stl_object
---

# `surface_stl_object`

The triangulated surface — the library's centrepiece. A `surface_stl_object`
owns:

- a private array of [`facet_object`](/guide/api/facet-object) holding the
  geometry,
- a private [`aabb_tree_object`](/guide/api/aabb-tree-object) accelerating
  distance and inside queries,
- a private unique-vertex pool deduplicating coincident vertices,
- cached scalar properties (volume, centroid, bounding box, manifold counts)
  populated by [`analyze`](#analyze).

Every public component is `private`. All state is read through getters
([`get_facets_number`](#getters-and-predicates),
[`get_bmin`](#getters-and-predicates),
[`get_volume`](#getters-and-predicates), …) and mutated through type-bound
procedures ([`load_from_file`](#load_from_file),
[`sanitize`](#sanitize),
[`translate`](#translate), …).

The TBPs below are grouped by user-flow: I/O, repair, analysis, queries,
manipulation, and a final catch-all section for the read-only accessors and
validity predicates.

[[toc]]

---

## I/O

### `load_from_file`

```fortran
call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
```

**Purpose.** Read an STL file from disk into the surface. On return the surface
holds the parsed facets and a full set of cached scalars: bounding box,
connectivity, volume, centroid, AABB tree, and per-vertex / per-edge
pseudo-normals. In other words, `load_from_file` is `read + analyze` fused
into one call — you can immediately query [`distance`](#distance) without an
intermediate [`analyze`](#analyze).

The format is detected automatically when `guess_format=.true.`: the file size
is compared against the size implied by the binary header's facet count, and
the format that matches wins. Pass `is_ascii=.true./.false.` explicitly if you
already know the format and want to skip detection.

**Arguments.**

| Argument                  | Intent       | Type               | Notes                                                                                              |
| ------------------------- | ------------ | ------------------ | -------------------------------------------------------------------------------------------------- |
| `file_name`               | `in`         | `character(*)`     | (required) Path to the STL file.                                                                   |
| `is_ascii`                | `in`, opt.   | `logical`          | Force ASCII or binary parsing. Ignored when `guess_format=.true.`. Defaults to ASCII otherwise.     |
| `guess_format`            | `in`, opt.   | `logical`          | Auto-detect ASCII vs binary from the file size. Recommended for end-user-supplied paths.            |
| `clip_min`, `clip_max`    | `in`, opt.   | `type(vector_R8P)` | Both must be present together. Facets entirely outside the AABB `[clip_min, clip_max]` are dropped during load. Useful for opening a single subregion of a huge file without paying for full I/O. |
| `aabb_refinement_levels`  | `in`, opt.   | `integer(I4P)`     | Octree depth, or `AABB_AUTO_REFINEMENT`. Ignored when the kind is `AABB_TREE_SAH_BVH`.              |
| `aabb_tree_kind`          | `in`, opt.   | `integer(I4P)`     | `AABB_TREE_SAH_BVH` (default) or `AABB_TREE_OCTREE`.                                                |
| `status`                  | `out`, opt.  | `integer(I4P)`     | `STATUS_OK`, `STATUS_FILE_NOT_FOUND`, `STATUS_INVALID_INPUT`, or `STATUS_ALLOC_FAIL`.                |

**Example.** Load the bundled cube, print summary, and check it's a valid solid.

```fortran
program ex_load_from_file
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface
integer(I4P)             :: stat

call surface%load_from_file(file_name='cube.stl', guess_format=.true., status=stat)
if (stat /= STATUS_OK) error stop 'load failed'
print '(A, I0)',    'facets : ', surface%get_facets_number()
print '(A, ES12.5)','volume : ', surface%get_volume()
print '(A, L1)',    'is_volume? ', surface%is_volume()
end program ex_load_from_file
```

**Example.** Load with an octree (overriding the default BVH).

```fortran
program ex_load_octree
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.,        &
                            aabb_tree_kind=AABB_TREE_OCTREE,                  &
                            aabb_refinement_levels=AABB_AUTO_REFINEMENT)
print '(A, I0)', 'octree levels = ', surface%aabb%get_refinement_levels()
end program ex_load_octree
```

**Notes.**

- NaN / Inf coordinate values cause an immediate `STATUS_INVALID_INPUT`; the
  surface is left untouched. This is policy: load refuses corrupt input rather
  than passing it downstream to be mis-repaired.
- A successful load **does not** imply a clean mesh. The file may still have
  degenerate, duplicate, or inconsistently-oriented facets. Call
  [`sanitize`](#sanitize) immediately after loading user-supplied input.
- Repeated calls reuse the surface — the previous content is destroyed first.

**See also.** [`save_into_file`](#save_into_file),
[`sanitize`](#sanitize),
[`analyze`](#analyze),
[`AABB_TREE_SAH_BVH`](/guide/api/constants#aabb_tree_sah_bvh),
[`STATUS_*`](/guide/api/constants#status-codes).

---

### `save_into_file`

```fortran
call surface%save_into_file(file_name='out.stl', is_ascii=.false.)
```

**Purpose.** Write the surface to disk in STL format. The current header (set
via [`set_header`](#getters-and-predicates) or carried over from
[`load_from_file`](#load_from_file)) is written first.

**Arguments.**

| Argument    | Intent       | Type             | Notes                                                              |
| ----------- | ------------ | ---------------- | ------------------------------------------------------------------ |
| `file_name` | `in`         | `character(*)`   | (required) Output path. The file is truncated.                     |
| `is_ascii`  | `in`, opt.   | `logical`        | Default `.true.` (ASCII output). Set `.false.` for binary STL.     |
| `status`    | `out`, opt.  | `integer(I4P)`   | `STATUS_OK` or `STATUS_FILE_OPEN_FAIL`.                            |

**Example.** Round-trip through binary STL.

```fortran
program ex_save_into_file
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%save_into_file(file_name='cube-binary.stl', is_ascii=.false.)
call surface%save_into_file(file_name='cube-ascii.stl',  is_ascii=.true.)
end program ex_save_into_file
```

**Notes.**

- Binary is ~6× smaller and ~20× faster to read back than ASCII. Use it for
  intermediate data; reserve ASCII for human inspection.
- The STL format is lossy w.r.t. connectivity — the file stores only triangles.
  After a save / load round trip you must call [`analyze`](#analyze) (or rely
  on `load_from_file`'s internal analyse) to rebuild connectivity.

**See also.** [`load_from_file`](#load_from_file),
[`set_header`](#getters-and-predicates).

---

## Repair pipeline

Real-world STL files are routinely broken — zero-area slivers, duplicate
triangles, unwelded vertices, mixed winding. FOSSIL's repair API is layered:
a single entry point ([`sanitize`](#sanitize)) that orchestrates every pass in
the right order, plus the individual passes exposed so you can run them
manually when the orchestrator's policy is wrong for your data.

### `sanitize`

```fortran
call surface%sanitize
```

**Purpose.** Full repair pipeline. Run this immediately after every
user-supplied load. The pipeline is the result of a long-running iteration on
the order in which repair passes must compose; running them in a different
order can leave subtle corruptions.

The fixed pipeline is:

1. `remove_degenerate_facets` — drop zero-area / sliver triangles first so
   their NaN normals do not poison later steps.
2. `connect_nearby_vertices` — only if disconnected edges remain. Snaps
   coincident vertices via union-find on a tolerance.
3. `analyze` — rebuild bounding box, connectivity, volume, centroid, AABB tree.
4. `remove_duplicate_facets` — drop literal duplicates (orientation-agnostic).
5. `analyze` again, if step 4 removed anything, so winding fixup starts from
   a clean state.
6. `sanitize_normals` — BFS-propagate winding consistency from facet 1, then
   flip globally if the signed volume came out negative.
7. **Warnings** — one line to `stderr` for every defect class with non-zero
   count. Suppress them by redirecting unit `error_unit`; preserve them in
   production logs so you notice when your inputs degrade.

**Arguments.**

| Argument       | Intent       | Type           | Notes                                                                                |
| -------------- | ------------ | -------------- | ------------------------------------------------------------------------------------ |
| `do_analysis`  | `in`, opt.   | `logical`      | If `.true.`, run an initial `analyze` before step 1. Default `.false.` — the caller normally arrives here straight from `load_from_file`, which has already analysed. |
| `status`       | `out`, opt.  | `integer(I4P)` | `STATUS_OK`; reserved for future status codes.                                       |

**Example.** Load a deliberately broken mesh and verify repair.

```fortran
program ex_sanitize
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube-inconsistent.stl', guess_format=.true.)
print '(A, L1)', 'is_volume before sanitize? ', surface%is_volume()
call surface%sanitize
print '(A, L1)', 'is_volume after  sanitize? ', surface%is_volume()
print '(A, I0)', 'degenerate facets removed = ', surface%get_degenerate_facets_removed()
print '(A, I0)', 'duplicate  facets removed = ', surface%get_duplicate_facets_removed()
end program ex_sanitize
```

**Notes.**

- `sanitize` is **idempotent on a clean mesh**: a second call is a no-op
  (every pass detects zero defects). Use this property to assert cleanliness
  in tests.
- The connect / duplicate-removal tolerances are derived from the bounding-box
  diagonal. If your mesh has wildly different feature scales, run the
  individual passes manually with explicit tolerances.

**See also.** [`remove_degenerate_facets`](#remove_degenerate_facets),
[`remove_duplicate_facets`](#remove_duplicate_facets),
[`connect_nearby_vertices`](#connect_nearby_vertices),
[`sanitize_normals`](#sanitize_normals),
[`is_volume`](#getters-and-predicates).

---

### `remove_degenerate_facets`

```fortran
call surface%remove_degenerate_facets
```

**Purpose.** Drop facets whose area is below a tolerance relative to the
mesh's bounding-box diagonal. Degenerate facets have undefined normals (often
NaN) and contaminate every downstream computation; removing them is the first
step of any repair pipeline.

**Arguments.** None.

**Example.**

```fortran
program ex_remove_degenerate
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%remove_degenerate_facets
print '(A, I0)', 'degenerate facets removed = ', surface%get_degenerate_facets_removed()
end program ex_remove_degenerate
```

**Notes.**

- The count is queryable through
  [`get_degenerate_facets_removed`](#getters-and-predicates) for post-run
  inspection.
- Removal does not re-run [`analyze`](#analyze) — call it yourself if you need
  fresh connectivity afterwards. The orchestrator [`sanitize`](#sanitize) does
  this automatically.

**See also.** [`remove_duplicate_facets`](#remove_duplicate_facets),
[`sanitize`](#sanitize).

---

### `remove_duplicate_facets`

```fortran
call surface%remove_duplicate_facets
```

**Purpose.** Drop facets that duplicate another facet up to vertex permutation
(including reversed winding). Duplicates skew volume, double-count surface
area, and bias the winding-consistency BFS.

**Arguments.** None.

**Example.**

```fortran
program ex_remove_duplicate
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%remove_duplicate_facets
print '(A, I0)', 'duplicate facets removed = ', surface%get_duplicate_facets_removed()
end program ex_remove_duplicate
```

**Notes.**

- Should run *after* [`remove_degenerate_facets`](#remove_degenerate_facets),
  because degenerate slivers tend to come in duplicate pairs.
- The detection key is the *sorted triple of unique vertex ids*; this catches
  reversed-winding duplicates that a naïve vertex-coordinate comparison would
  miss.

**See also.** [`sanitize`](#sanitize),
[`remove_degenerate_facets`](#remove_degenerate_facets).

---

### `connect_nearby_vertices`

```fortran
call surface%connect_nearby_vertices
```

**Purpose.** Snap coincident-but-not-shared vertices via union-find. STL is
soup-of-triangles by spec: nothing in the format enforces that two triangles
sharing an edge use the same vertex objects. After this pass, the
unique-vertex pool reflects geometric coincidence, and connectivity (computed
by [`analyze`](#analyze)) can resolve all the edges that visually look shared
but were previously disconnected.

**Arguments.** None.

**Example.**

```fortran
program ex_connect_nearby
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube_disconnect.stl', guess_format=.true.)
call surface%connect_nearby_vertices
call surface%analyze
print '(A, L1)', 'manifold after connect? ', surface%is_manifold()
end program ex_connect_nearby
```

**Notes.**

- Uses a bounding-box-diagonal-relative tolerance. The tolerance is fixed
  internally; for adversarial inputs with feature-scale separation you should
  instead pre-snap your vertices.
- Always re-run [`analyze`](#analyze) afterwards (the orchestrator does this
  automatically).

**See also.** [`sanitize`](#sanitize),
[`analyze`](#analyze).

---

### `sanitize_normals`

```fortran
call surface%sanitize_normals
```

**Purpose.** Make every facet's winding consistent with its neighbours, then
flip the whole surface if the signed volume is negative — yielding outward
normals on a closed body. This is the precondition for the default
[`SIGN_PSEUDO_NORMAL`](/guide/api/constants#sign_pseudo_normal) sign algorithm.

**Arguments.** None.

**Example.**

```fortran
program ex_sanitize_normals
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube-inconsistent.stl', guess_format=.true.)
call surface%sanitize_normals
print '(A, ES12.5)', 'volume after normals fix = ', surface%get_volume()
end program ex_sanitize_normals
```

**Notes.**

- The BFS starts at facet 1; if your mesh has disconnected components, each
  component is walked independently and the global sign is decided per
  component, by component-signed-volume.
- For open shells (non-closed surfaces), "outward" is ill-defined and the
  global flip is skipped; the BFS still produces locally-consistent winding.

**See also.** [`reverse_normals`](#reverse_normals),
[`SIGN_PSEUDO_NORMAL`](/guide/api/constants#sign_pseudo_normal),
[`sanitize`](#sanitize).

---

### `reverse_normals`

```fortran
call surface%reverse_normals
```

**Purpose.** Flip every facet's normal and reverse its vertex winding. Useful
when you load a mesh that is correctly *consistent* but oriented inward
(volume negative) and you want to fix it without running the full BFS in
[`sanitize_normals`](#sanitize_normals).

**Arguments.** None.

**Example.**

```fortran
program ex_reverse_normals
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface
real(R8P)                :: v_before, v_after

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
v_before = surface%get_volume()
call surface%reverse_normals
call surface%analyze
v_after  = surface%get_volume()
print '(A, ES12.5, A, ES12.5)', 'volume: before=', v_before, ' after=', v_after
end program ex_reverse_normals
```

**Notes.**

- Re-run [`analyze`](#analyze) afterwards to refresh the signed volume.

**See also.** [`sanitize_normals`](#sanitize_normals).

---

## Analysis

### `analyze`

```fortran
call surface%analyze
```

**Purpose.** Compute, in one pass, every derived quantity that the surface
caches: bounding box, edge connectivity, signed volume, centroid, AABB tree,
and per-edge / per-vertex pseudo-normals. Called automatically by
[`load_from_file`](#load_from_file), [`sanitize`](#sanitize), and
[`clip`](#clip); you only call it manually after a transformation that you
chose **not** to follow with `recompute_metrix=.true.`.

**Arguments.**

| Argument                 | Intent       | Type           | Notes                                                                              |
| ------------------------ | ------------ | -------------- | ---------------------------------------------------------------------------------- |
| `aabb_refinement_levels` | `in`, opt.   | `integer(I4P)` | Octree depth or `AABB_AUTO_REFINEMENT`. Ignored by `AABB_TREE_SAH_BVH`.             |
| `aabb_tree_kind`         | `in`, opt.   | `integer(I4P)` | `AABB_TREE_SAH_BVH` or `AABB_TREE_OCTREE`. Defaults to the surface's current kind.  |
| `status`                 | `out`, opt.  | `integer(I4P)` | `STATUS_OK` or `STATUS_ALLOC_FAIL`.                                                |

**Example.** Translate, then re-analyse to refresh the AABB tree and centroid.

```fortran
program ex_analyze
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface

type(vector_R8P)         :: c

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%translate(delta=10._R8P * ex_R8P)
call surface%analyze
c = surface%get_centroid()
print '(A, 3ES12.5)', 'new centroid = ', c%x, c%y, c%z
end program ex_analyze
```

**Notes.**

- Switching `aabb_tree_kind` between calls *does* rebuild the tree.
- Pure rotations / translations / mirrors leave connectivity intact; if
  performance matters you can call `recompute_metrix=.true.` on the
  transformation directly and avoid the full analyse.

**See also.** [`build_connectivity`](#build_connectivity),
[`statistics`](#statistics),
[`AABB_TREE_SAH_BVH`](/guide/api/constants#aabb_tree_sah_bvh).

---

### `build_connectivity`

```fortran
call surface%build_connectivity
```

**Purpose.** Rebuild only the facet edge-adjacency (which facet shares which
edge with which). Cheaper than [`analyze`](#analyze) when you don't need the
full bounding box / volume / AABB refresh.

**Arguments.** None.

**Example.**

```fortran
program ex_build_connectivity
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%build_connectivity
print '(A, I0)', 'non-manifold edges = ', surface%get_non_manifold_edges_number()
end program ex_build_connectivity
```

**See also.** [`analyze`](#analyze).

---

### `statistics`

```fortran
print '(A)', surface%statistics()
```

**Purpose.** Return a multi-line formatted string describing the surface:
bounding box, volume, centroid, facet count, manifold / degenerate /
duplicate counts, and AABB-tree configuration. The format is human-readable;
do not parse it.

**Arguments.** None.

**Example.**

```fortran
program ex_statistics
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
print '(A)', surface%statistics()
end program ex_statistics
```

**See also.** [`analyze`](#analyze).

---

## Distance and inside queries

### `distance`

```fortran
d = surface%distance(point=p)
```

**Purpose.** Return the (minimum) distance from a point to the triangulated
surface. By default the result is **unsigned and squared** — useful for
nearest-facet ranking without a square-root cost. Pass `is_square_root=.true.`
for Euclidean distance and `is_signed=.true.` to apply one of the
[`SIGN_*`](/guide/api/constants#sign-algorithms) sign algorithms.

**Arguments.**

| Argument          | Intent       | Type               | Notes                                                                                          |
| ----------------- | ------------ | ------------------ | ---------------------------------------------------------------------------------------------- |
| `point`           | `in`         | `type(vector_R8P)` | (required) Query point.                                                                        |
| `is_signed`       | `in`, opt.   | `logical`          | Default `.false.` (unsigned).                                                                  |
| `sign_algorithm`  | `in`, opt.   | `integer(I4P)`     | One of `SIGN_PSEUDO_NORMAL` (default), `SIGN_RAY_INTERSECTIONS`, `SIGN_SOLID_ANGLE`.            |
| `is_square_root`  | `in`, opt.   | `logical`          | Default `.false.` — returns d². Set `.true.` for Euclidean distance.                          |

**Returns.** `real(R8P)` — the distance.

**Example.** All four common combinations from a single query point.

```fortran
program ex_distance
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: p

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
p = 2._R8P * ex_R8P
print '(A, ES12.5)', 'unsigned d^2     = ', surface%distance(point=p)
print '(A, ES12.5)', 'unsigned d       = ', surface%distance(point=p, is_square_root=.true.)
print '(A, ES12.5)', 'signed   d (pn)  = ', surface%distance(point=p, is_signed=.true.,         &
                                                              is_square_root=.true.)
print '(A, ES12.5)', 'signed   d (ray) = ', surface%distance(point=p, is_signed=.true.,         &
                                                              sign_algorithm=SIGN_RAY_INTERSECTIONS,&
                                                              is_square_root=.true.)
end program ex_distance
```

**Notes.**

- This is a thin wrapper over [`compute_distance`](#compute_distance) — use
  that one when you want the closest facet / edge / vertex indices in
  addition to the distance.
- Tree-accelerated: O(log N) expected on uniform meshes. The traversal is
  best-first with d² pruning, so it never visits a node whose AABB is farther
  than the current best.
- Bit-exact against brute-force on the regression suite, by construction —
  the pruning predicate is conservative.

**See also.** [`compute_distance`](#compute_distance),
[`is_point_inside`](#is_point_inside),
[`SIGN_*`](/guide/api/constants#sign-algorithms).

---

### `compute_distance`

```fortran
call surface%compute_distance(point=p, distance=d)
```

**Purpose.** Same as [`distance`](#distance) but with the distance returned
through an `intent(out)` argument plus four optional outputs that locate the
closest point on the mesh: facet index, local edge index, and local vertex
index. Use this when you need to know *which* feature of the mesh produced the
closest distance — for visualisation, classification, or as input to a
projection.

**Arguments.**

| Argument          | Intent       | Type               | Notes                                                                                                 |
| ----------------- | ------------ | ------------------ | ----------------------------------------------------------------------------------------------------- |
| `point`           | `in`         | `type(vector_R8P)` | (required) Query point.                                                                               |
| `distance`        | `out`        | `real(R8P)`        | (required) Distance from `point` to surface.                                                          |
| `is_signed`       | `in`, opt.   | `logical`          | Default `.false.`.                                                                                    |
| `sign_algorithm`  | `in`, opt.   | `integer(I4P)`     | `SIGN_*`; default `SIGN_PSEUDO_NORMAL`.                                                               |
| `is_square_root`  | `in`, opt.   | `logical`          | Default `.false.` — d² is returned.                                                                  |
| `facet_index`     | `out`, opt.  | `integer(I4P)`     | 1-based index of the facet holding the closest point. 0 if the surface has no facets.                 |
| `edge_index`      | `out`, opt.  | `integer(I4P)`     | Local edge id (1, 2, or 3) when the closest point lies on an edge interior; 0 otherwise.              |
| `vertex_index`    | `out`, opt.  | `integer(I4P)`     | Local vertex id (1, 2, or 3) when the closest point coincides with a vertex; 0 otherwise.             |

**Example.** Identify the closest feature class on every query.

```fortran
program ex_compute_distance
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: p
real(R8P)                :: d
integer(I4P)             :: fid, eid, vid
character(len=8)         :: feature

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
p = 2._R8P * ex_R8P
call surface%compute_distance(point=p, distance=d, is_signed=.true., is_square_root=.true., &
                              facet_index=fid, edge_index=eid, vertex_index=vid)
if (vid /= 0) then
   feature = 'vertex'
else if (eid /= 0) then
   feature = 'edge'
else
   feature = 'face'
end if
print '(A,ES12.5,A,A,A,I0)', 'd=', d, ' closest feature=', trim(feature), &
                              ' on facet ', fid
end program ex_compute_distance
```

**Notes.**

- The `edge_index` / `vertex_index` outputs are mutually exclusive with
  "interior face": at most one of them is non-zero for a given query.
- When `sign_algorithm = SIGN_PSEUDO_NORMAL`, sign and closest-feature are
  computed in a single tree traversal; with `SIGN_RAY_INTERSECTIONS` or
  `SIGN_SOLID_ANGLE` a second O(N) pass classifies the sign.

**See also.** [`distance`](#distance),
[`is_point_inside`](#is_point_inside),
[`facet%compute_distance_with_region`](/guide/api/facet-object#compute_distance_with_region).

---

### `is_point_inside`

```fortran
inside = surface%is_point_inside(point=p)
```

**Purpose.** Return `.true.` if the query point lies strictly inside the
closed surface. A wrapper over [`distance`](#distance) with `is_signed=.true.`
that returns only the sign — slightly more readable at call sites that don't
care about the magnitude.

**Arguments.**

| Argument          | Intent       | Type               | Notes                                                              |
| ----------------- | ------------ | ------------------ | ------------------------------------------------------------------ |
| `point`           | `in`         | `type(vector_R8P)` | (required) Query point.                                            |
| `sign_algorithm`  | `in`, opt.   | `integer(I4P)`     | `SIGN_*`; default `SIGN_PSEUDO_NORMAL`.                            |

**Returns.** `logical` — `.true.` inside, `.false.` outside.

**Example.**

```fortran
program ex_is_point_inside
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
print '(A, L1)', 'centre inside? ', surface%is_point_inside(point=surface%get_centroid())
end program ex_is_point_inside
```

**Notes.**

- Behaviour on the boundary is **algorithm-dependent**: pseudo-normal returns
  whichever side wins the dot product (numerically unstable on the surface
  itself); ray-cast is technically undefined on the surface. Treat the result
  as meaningless when the query point is closer to the surface than the
  mesh's floating-point resolution.
- On open shells, `SIGN_PSEUDO_NORMAL` returns a well-defined value but the
  notion of "inside" is mesh-dependent. Prefer `SIGN_RAY_INTERSECTIONS` for
  open shells if you trust their orientation.

**See also.** [`distance`](#distance),
[`compute_distance`](#compute_distance),
[`is_watertight`](#getters-and-predicates).

---

## Manipulation

### `translate`

```fortran
call surface%translate(delta=vec)
call surface%translate(x=1._R8P, y=2._R8P, z=3._R8P)
```

**Purpose.** Rigid-body translation. Pass **either** a `vector_R8P` `delta`
**or** any subset of scalar `x`, `y`, `z` — mixing the two forms returns
`STATUS_AMBIGUOUS_ARGS`.

**Arguments.**

| Argument            | Intent       | Type                | Notes                                                                              |
| ------------------- | ------------ | ------------------- | ---------------------------------------------------------------------------------- |
| `x`, `y`, `z`       | `in`, opt.   | `real(R8P)`         | Per-axis increments. Missing axes are not translated.                              |
| `delta`             | `in`, opt.   | `type(vector_R8P)`  | Vectorial increment.                                                               |
| `recompute_metrix`  | `in`, opt.   | `logical`           | Default `.false.`. Pass `.true.` to refresh facet plane equations in place.        |
| `status`            | `out`, opt.  | `integer(I4P)`      | `STATUS_OK` or `STATUS_AMBIGUOUS_ARGS`.                                            |

**Example.**

```fortran
program ex_translate
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface

type(vector_R8P)         :: c

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%translate(delta=10._R8P * ex_R8P, recompute_metrix=.true.)
c = surface%get_centroid()
print '(A, 3ES12.5)', 'centroid = ', c%x, c%y, c%z
end program ex_translate
```

**Notes.**

- If you batch multiple transforms, pass `recompute_metrix=.false.` on every
  call but the last; the final one (or a closing [`analyze`](#analyze)) brings
  the cached state up to date.
- Translation does not change connectivity — `analyze` is not strictly
  required, but the bounding box, centroid, and AABB tree become stale.

**See also.** [`rotate`](#rotate),
[`mirror`](#mirror),
[`resize`](#resize).

---

### `rotate`

```fortran
call surface%rotate(axis=ax, angle=theta)
call surface%rotate(matrix=R)
call surface%rotate(axis=ax, angle=theta, center=surface%get_centroid())
```

**Purpose.** Rigid-body rotation. Either give an axis (a `vector_R8P` not
required to be normalised) and an angle in radians, or pass a 3×3 rotation
matrix directly. Pass `center=` to rotate about an arbitrary pivot — by
default the rotation pivots about the world origin.

**Arguments.**

| Argument            | Intent       | Type                | Notes                                                                                |
| ------------------- | ------------ | ------------------- | ------------------------------------------------------------------------------------ |
| `axis`              | `in`         | `type(vector_R8P)`  | (axis-angle form) Rotation axis.                                                     |
| `angle`             | `in`         | `real(R8P)`         | (axis-angle form) Angle in radians.                                                  |
| `matrix`            | `in`         | `real(R8P)(3,3)`    | (matrix form) Rotation matrix. Mutually exclusive with `axis`/`angle`.               |
| `center`            | `in`, opt.   | `type(vector_R8P)`  | Pivot. Defaults to the world origin. Pass `surface%get_centroid()` for body-frame.   |
| `recompute_metrix`  | `in`, opt.   | `logical`           | Default `.false.`.                                                                   |

**Example.** Body-frame rotation around the centroid (issue #6).

```fortran
program ex_rotate
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ez_R8P
implicit none
type(surface_stl_object) :: surface
real(R8P),    parameter  :: pi = acos(-1._R8P)

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%rotate(axis=ez_R8P, angle=0.25_R8P * pi, &
                    center=surface%get_centroid(), recompute_metrix=.true.)
end program ex_rotate
```

**Notes.**

- The axis is normalised internally; magnitude is irrelevant.
- World-origin rotation moves the centroid; body-frame rotation does not. If
  you want "rotate in place", always pass `center=surface%get_centroid()`.

**See also.** [`translate`](#translate),
[`mirror`](#mirror).

---

### `mirror`

```fortran
call surface%mirror(normal=plane_normal)
call surface%mirror(matrix=M)
```

**Purpose.** Mirror the surface through a plane. The plane is specified
either by its unit normal (passing through the origin) or by an explicit 3×3
reflection matrix. Vertex winding is reversed automatically so normals stay
geometrically meaningful.

**Arguments.**

| Argument            | Intent       | Type                | Notes                                                                |
| ------------------- | ------------ | ------------------- | -------------------------------------------------------------------- |
| `normal`            | `in`         | `type(vector_R8P)`  | (normal form) Plane normal. Plane passes through the origin.         |
| `matrix`            | `in`         | `real(R8P)(3,3)`    | (matrix form) Reflection matrix.                                     |
| `recompute_metrix`  | `in`, opt.   | `logical`           | Default `.false.`.                                                   |

**Example.**

```fortran
program ex_mirror
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%mirror(normal=ex_R8P, recompute_metrix=.true.)
end program ex_mirror
```

**Notes.**

- To mirror through a plane *not* through the origin, translate the surface
  so the plane passes through the origin, mirror, translate back.

**See also.** [`rotate`](#rotate),
[`translate`](#translate).

---

### `resize`

```fortran
call surface%resize(factor=2._R8P * ex_R8P + ey_R8P + ez_R8P)
call surface%resize(x=2._R8P, y=2._R8P, z=2._R8P)
call surface%resize(factor=2._R8P*ex_R8P+ey_R8P+ez_R8P, respect_centroid=.true.)
```

**Purpose.** Anisotropic scaling. As with [`translate`](#translate), pass
either a vector `factor` or scalar `x`/`y`/`z` — never both. By default the
scale is around the world origin; pass `respect_centroid=.true.` to pivot
about the current centroid (the scale is then "in place").

**Arguments.**

| Argument            | Intent       | Type                | Notes                                                                            |
| ------------------- | ------------ | ------------------- | -------------------------------------------------------------------------------- |
| `x`, `y`, `z`       | `in`, opt.   | `real(R8P)`         | Per-axis factors. Missing axes are not scaled (effective factor 1).              |
| `factor`            | `in`, opt.   | `type(vector_R8P)`  | Vectorial factor.                                                                |
| `respect_centroid`  | `in`, opt.   | `logical`           | Default `.false.`. Pivot at centroid when `.true.`.                              |
| `recompute_metrix`  | `in`, opt.   | `logical`           | Default `.false.`.                                                               |
| `status`            | `out`, opt.  | `integer(I4P)`      | `STATUS_OK` or `STATUS_AMBIGUOUS_ARGS`.                                          |

**Example.**

```fortran
program ex_resize
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : ex_R8P, ey_R8P, ez_R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%resize(factor=2._R8P*ex_R8P + 2._R8P*ey_R8P + 2._R8P*ez_R8P, &
                    respect_centroid=.true., recompute_metrix=.true.)
end program ex_resize
```

**Notes.**

- The method is called `resize` and not `scale` because `scale` is a Fortran
  intrinsic.
- Negative factors mirror across the corresponding axis. Use
  [`mirror`](#mirror) instead if that is your intent — it adjusts winding.

**See also.** [`translate`](#translate),
[`mirror`](#mirror).

---

### `clip`

```fortran
call surface%clip(bmin=lo, bmax=hi)
call surface%clip(bmin=lo, bmax=hi, remainder=outside)
```

**Purpose.** Drop facets that do not overlap the AABB `[bmin, bmax]`,
keeping only the part of the surface inside (or intersecting) the box. With
the optional `remainder` argument the discarded facets are moved into a fresh
surface so nothing is lost.

**Arguments.**

| Argument     | Intent       | Type                       | Notes                                                                |
| ------------ | ------------ | -------------------------- | -------------------------------------------------------------------- |
| `bmin`       | `in`         | `type(vector_R8P)`         | (required) Clip box minimum corner.                                  |
| `bmax`       | `in`         | `type(vector_R8P)`         | (required) Clip box maximum corner.                                  |
| `remainder`  | `out`, opt.  | `type(surface_stl_object)` | Receives the facets that were dropped.                               |
| `status`     | `out`, opt.  | `integer(I4P)`             | `STATUS_OK` or `STATUS_ALLOC_FAIL`.                                  |

**Example.** Split a surface into two halves on the xy-plane.

```fortran
program ex_clip
use fossil
use penf, only : I4P, R8P
use vecfor_R8P, only : vector_R8P
implicit none
type(surface_stl_object) :: surface, top
type(vector_R8P)         :: lo, hi

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
lo = surface%get_bmin()
hi = surface%get_bmax()
lo%z = 0.5_R8P * (lo%z + hi%z)
call surface%clip(bmin=lo, bmax=hi, remainder=top)
print '(A,I0,A,I0)', 'lower half facets=', surface%get_facets_number(), &
                     '  upper half facets=', top%get_facets_number()
end program ex_clip
```

**Notes.**

- Clip is **inclusion-by-overlap**, not by-containment: facets that straddle
  the box boundary go into the *inside* surface. There is no triangle
  subdivision along the clipping planes; if you need precise cuts, use a CSG
  boolean operation (planned — see roadmap issue #18).
- Both the kept and remainder surfaces are re-analysed automatically.

**See also.** [`merge_solids`](#merge_solids).

---

### `merge_solids`

```fortran
call surface%merge_solids(other=second_surface)
```

**Purpose.** Append every facet from `other` to `self`. The merged surface is
re-analysed afterwards. **Geometry is not de-duplicated** — overlapping
regions remain as overlapping triangles; the merge is concatenation, not a
boolean union.

**Arguments.**

| Argument | Intent       | Type                       | Notes                                                              |
| -------- | ------------ | -------------------------- | ------------------------------------------------------------------ |
| `other`  | `in`         | `type(surface_stl_object)` | (required) The donor surface. Unchanged on return.                 |
| `status` | `out`, opt.  | `integer(I4P)`             | `STATUS_OK` or `STATUS_ALLOC_FAIL`.                                |

**Example.** Reassemble a dragon from two halves.

```fortran
program ex_merge_solids
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: dragon, half2

call dragon%load_from_file(file_name='dragon_part_1.stl', guess_format=.true.)
call half2 %load_from_file(file_name='dragon_part_2.stl', guess_format=.true.)
call dragon%merge_solids(other=half2)
print '(A, I0)', 'merged facets = ', dragon%get_facets_number()
end program ex_merge_solids
```

**Notes.**

- A subsequent [`sanitize`](#sanitize) is recommended: the seam between the
  two donor meshes typically produces near-coincident vertices and duplicate
  facets.

**See also.** [`clip`](#clip),
[`sanitize`](#sanitize).

---

## Getters and predicates

The read-only accessors are listed here as a single table — they are
self-explanatory and do not benefit from the full per-method template. Every
one of them is `pure`.

| Method                                | Returns                                                                                |
| ------------------------------------- | -------------------------------------------------------------------------------------- |
| `get_facets_number()`                 | `integer(I4P)` — number of facets.                                                     |
| `get_bmin()`, `get_bmax()`            | `type(vector_R8P)` — min / max corners of the axis-aligned bounding box.               |
| `get_volume()`                        | `real(R8P)` — signed volume. Positive for closed bodies with outward normals.          |
| `get_area()`                          | `real(R8P)` — total surface area (sum of facet areas). Added in issue #7.              |
| `get_centroid()`                      | `type(vector_R8P)` — surface centroid.                                                 |
| `get_header()`                        | `character(FRLEN)` — STL header string.                                                |
| `set_header(header)`                  | Set the header for the next `save_into_file`.                                          |
| `get_non_manifold_edges_number()`     | `integer(I4P)` — edges shared by 3+ facets.                                            |
| `get_degenerate_facets_removed()`     | `integer(I4P)` — count from the most recent `remove_degenerate_facets` pass.           |
| `get_duplicate_facets_removed()`      | `integer(I4P)` — count from the most recent `remove_duplicate_facets` pass.            |
| `get_vertex_pool()`                   | Read-only pointer to the unique-vertex pool.                                           |
| `facet_at(i)`                         | Pointer to facet `i`; `null()` if `i` is out of range.                                 |
| `facets_ref()`                        | Pointer to the entire facet array.                                                     |

And the three composite validity predicates:

| Predicate          | True iff                                                                                          |
| ------------------ | ------------------------------------------------------------------------------------------------- |
| `is_watertight()`  | Every edge has exactly two incident facets (no boundary, no non-manifold edges).                  |
| `is_manifold()`    | No edge is shared by three or more facets (boundary edges are allowed — open shells qualify).     |
| `is_volume()`      | `is_watertight()` **and** signed volume > 0 **and** centroid is finite.                           |

**Idiomatic usage.** Read the predicates after [`sanitize`](#sanitize) to
classify the input:

```fortran
program ex_predicates
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
if (.not. surface%is_manifold()) then
   print '(A)', 'mesh has non-manifold edges; distance queries may be unreliable'
elseif (.not. surface%is_watertight()) then
   print '(A)', 'mesh is manifold but open; signed distance is mesh-dependent'
elseif (.not. surface%is_volume()) then
   print '(A)', 'mesh is closed but oriented inward; sanitize_normals should have fixed this'
else
   print '(A)', 'mesh is a valid solid'
end if
end program ex_predicates
```

**See also.** [`sanitize`](#sanitize),
[`statistics`](#statistics).
