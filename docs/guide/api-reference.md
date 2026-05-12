---
title: API Reference
---

# API Reference

FOSSIL exposes a single module:

```fortran
use fossil
```

The module re-exports two public types (`surface_stl_object`, `facet_object`), plus the named constants used to configure distance queries and the AABB tree. All geometry operations use [PENF](https://github.com/szaghi/PENF) numeric kinds (`I4P`, `R8P`) and [VecFor](https://github.com/szaghi/VecFor) vector types (`vector_R8P`).

## Public constants

### Sign algorithms (used by signed-distance queries)

| Constant | Description |
|---|---|
| `SIGN_PSEUDO_NORMAL` | **(default)** Bærentzen–Aanæs angle-weighted pseudo-normal test. Fused with the closest-facet traversal so signed distance costs the same as unsigned. Requires consistently outward-oriented normals — `sanitize` produces this state. |
| `SIGN_RAY_INTERSECTIONS` | Axis-aligned ray-cast: odd intersections = inside. Robust on meshes with mixed orientation; ~3× slower than pseudo-normal. |
| `SIGN_SOLID_ANGLE` | Sum of projected solid angles ≈ ±4π. Robust but unaccelerated (O(N) per query). |

`sign_algorithm_from_string(name)` parses `"pseudo_normal"`, `"ray_intersections"`, `"solid_angle"` (CLI/config use).

### AABB tree kind

| Constant | Description |
|---|---|
| `AABB_TREE_SAH_BVH` | **(default)** Binary BVH built top-down by partitioning triangles along the longest centroid-bbox axis, using the bucketed surface-area heuristic. Adapts to local triangle density. |
| `AABB_TREE_OCTREE` | Legacy 8-way space-partitioning octree, built to a uniform refinement depth (see `AABB_AUTO_REFINEMENT`). |

### Octree refinement levels

| Constant | Description |
|---|---|
| `AABB_AUTO_REFINEMENT` | Sentinel for the octree's auto-tune heuristic (depth picked from the facet count: `ceil(log8(N/64))` clamped to `[1, 6]`). Only meaningful for `AABB_TREE_OCTREE`. |

### Status codes

Returned via the optional `status` argument on mutating procedures (`load_from_file`, `clip`, `translate`, etc.):

| Constant | Returned when |
|---|---|
| `STATUS_OK` | Success (zero). |
| `STATUS_ALLOC_FAIL` | An internal allocation failed. |
| `STATUS_AMBIGUOUS_ARGS` | Conflicting optional arguments (e.g. both `delta=` and `x=`/`y=`/`z=` on `translate`). |
| `STATUS_FILE_NOT_FOUND` | `load_from_file` could not open the input path. |
| `STATUS_FILE_OPEN_FAIL` | `save_into_file` could not open the output path. |
| `STATUS_INVALID_INPUT` | Input contains NaN/Inf vertex coordinates; load refused. |

---

## `surface_stl_object`

The triangulated surface — holds the facet array plus connectivity, bounding box, volume, centroid, and AABB tree. All public state is accessed through getters; the components themselves are private.

### Read-only accessors

| Method | Returns |
|---|---|
| `get_facets_number()` | `integer(I4P)` — facet count |
| `get_bmin()`, `get_bmax()` | `type(vector_R8P)` — axis-aligned bounding-box corners |
| `get_volume()` | `real(R8P)` — signed volume (positive for outward-oriented closed bodies) |
| `get_centroid()` | `type(vector_R8P)` — surface centroid |
| `get_header()` | `character(FRLEN)` — STL file header |
| `get_non_manifold_edges_number()` | `integer(I4P)` — count of edges with 3+ incident facets |
| `get_degenerate_facets_removed()` | `integer(I4P)` — count removed by the last `remove_degenerate_facets` pass |
| `get_duplicate_facets_removed()` | `integer(I4P)` — count removed by the last `remove_duplicate_facets` pass |
| `facet_at(i)` | `pointer` to facet(i); `null()` if `i` is out of range |
| `facets_ref()` | `pointer` to the whole facet array |
| `aabb` (public type-bound) | the underlying `aabb_tree_object` (use `surface%aabb%get_tree_kind()`, etc.) |

### Validity predicates

| Method | True iff |
|---|---|
| `is_watertight()` | Every edge has exactly two incident facets (no boundary, no non-manifold) |
| `is_manifold()`   | No non-manifold edges (boundary edges allowed — open shells can be manifold) |
| `is_volume()`     | `is_watertight()` AND volume > 0 AND finite centroid |

### Loading / saving

#### [`load_from_file`](#load_from_file) {#load_from_file}

Load an STL file. Runs `analyze` internally — bounding box, connectivity, AABB tree, and pseudo-normals are all populated before return.

```fortran
call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
```

| Argument | Intent | Type | Description |
|---|---|---|---|
| `file_name` | `in` | `character(*)` | File path |
| `is_ascii` | `in`, optional | `logical` | Force ASCII format (ignored if `guess_format=.true.`) |
| `guess_format` | `in`, optional | `logical` | Auto-detect ASCII vs binary from file size |
| `clip_min`, `clip_max` | `in`, optional | `type(vector_R8P)` | Discard facets outside this AABB during load |
| `aabb_refinement_levels` | `in`, optional | `integer(I4P)` | Octree depth or `AABB_AUTO_REFINEMENT` (ignored by SAH BVH) |
| `aabb_tree_kind` | `in`, optional | `integer(I4P)` | `AABB_TREE_SAH_BVH` (default) or `AABB_TREE_OCTREE` |
| `status` | `out`, optional | `integer(I4P)` | `STATUS_OK` or one of `STATUS_FILE_NOT_FOUND` / `STATUS_INVALID_INPUT` / `STATUS_ALLOC_FAIL` |

#### `save_into_file` {#save_into_file}

```fortran
call surface%save_into_file(file_name='output.stl', is_ascii=.false., status=stat)
```

### Repair pipeline

| Method | Description |
|---|---|
| `sanitize` | Full repair orchestrator (see below) |
| `remove_degenerate_facets` | Drop zero-area / sliver triangles |
| `remove_duplicate_facets` | Drop literal duplicates (orientation-agnostic) |
| `connect_nearby_vertices` | Union-find vertex deduplication on a tolerance |
| `sanitize_normals` | BFS-propagate winding consistency, then flip globally so volume > 0 |
| `reverse_normals` | Flip every facet's normal (and winding) |

#### `sanitize`

The high-level entry point. Internally runs, in order:

1. `remove_degenerate_facets`
2. `connect_nearby_vertices` (only if disconnected edges remain)
3. `analyze`
4. `remove_duplicate_facets` (re-`analyze` if anything was removed)
5. `sanitize_normals`

Each pass emits a stderr warning if it detected any defects.

```fortran
call surface%sanitize
print '(A)', surface%statistics()
```

### Analysis

| Method | Description |
|---|---|
| `analyze` | Compute bounding box, connectivity, volume, centroid, AABB tree, pseudo-normals. Called automatically by `load_from_file`, `sanitize`, `clip`. |
| `compute_volume` | Recompute volume from the current facets. |
| `compute_centroid` | Recompute centroid from the current facets. |
| `build_connectivity` | Rebuild facet edge-adjacency only (sort-and-pair algorithm). |
| `statistics` | Return a multi-line formatted string with bounding extents, volume, centroid, facet count, disconnected/non-manifold/degenerate/duplicate counts, AABB refinement levels. |

`analyze` accepts the same optional `aabb_refinement_levels` and `aabb_tree_kind` arguments as `load_from_file`.

### Manipulation

| Method | Notes |
|---|---|
| `translate(delta=)` or `translate(x=, y=, z=)` | One form per call (mixing them sets `status=STATUS_AMBIGUOUS_ARGS`) |
| `rotate(axis=, angle=)` or `rotate(matrix=)` | Angle in radians |
| `mirror(normal=)` or `mirror(matrix=)` | |
| `resize(factor=)` or `resize(x=, y=, z=)`; `respect_centroid=.true.` to pivot at centroid | |
| `clip(bmin=, bmax=, remainder=)` | `remainder` is optional |
| `merge_solids(other=)` | Appends other's facets; call `analyze` afterwards |

### Distance / inside queries

#### `distance` {#distance}

```fortran
d  = surface%distance(point=p)                                            ! unsigned, squared
d  = surface%distance(point=p, is_square_root=.true.)                     ! unsigned, Euclidean
sd = surface%distance(point=p, is_signed=.true., is_square_root=.true.)   ! signed, Euclidean
sd = surface%distance(point=p, is_signed=.true.,                       &
                      sign_algorithm=SIGN_RAY_INTERSECTIONS)              ! override sign algo
```

| Argument | Intent | Type | Description |
|---|---|---|---|
| `point` | `in` | `type(vector_R8P)` | Query point |
| `is_signed` | `in`, optional | `logical` | Default `.false.` (unsigned) |
| `sign_algorithm` | `in`, optional | `integer(I4P)` | `SIGN_*` constant; default `SIGN_PSEUDO_NORMAL` |
| `is_square_root` | `in`, optional | `logical` | Default `.false.` (returns d²) |

Tree-accelerated; bit-exact vs brute force on the regression suite.

#### `compute_distance` {#compute_distance}

Same as `distance` but with an explicit `intent(out)` `distance` parameter plus optional outputs reporting which facet, edge, or vertex held the closest point:

```fortran
real(R8P)     :: d
integer(I4P)  :: fid, eid, vid
call surface%compute_distance(point=p, distance=d, is_signed=.true., is_square_root=.true., &
                              facet_index=fid, edge_index=eid, vertex_index=vid)
```

`edge_index` is the local edge (1, 2, or 3) when the closest point landed on an edge interior, else 0. `vertex_index` is the local vertex (1, 2, or 3) when the closest point coincided with a vertex, else 0.

#### `is_point_inside` {#is_point_inside}

```fortran
inside = surface%is_point_inside(point=p)
inside = surface%is_point_inside(point=p, sign_algorithm=SIGN_SOLID_ANGLE)
```

---

## `aabb_tree_object`

Embedded in `surface_stl_object%aabb`. You normally interact with it through the surface's analyze / load_from_file API; direct access is useful for benchmarking and diagnostics.

| Method | Description |
|---|---|
| `get_tree_kind()` | Returns `AABB_TREE_SAH_BVH` or `AABB_TREE_OCTREE` |
| `set_tree_kind(kind)` | Set the kind (invalidates the build; call `initialize` afterwards) |
| `get_refinement_levels()` | Octree depth (meaningless for BVH) |
| `set_refinement_levels(n)` | Set the depth or pass `AABB_AUTO_REFINEMENT` |
| `get_nodes_number()` | Total node count after build |
| `get_is_initialized()` | True if the tree has been built |
| `node_at(i)` | `pointer` to node `i` for tree-walking diagnostics |
| `set_use_index(flag)` | `AABB_USE_INDEX` (default, use tree) or `AABB_USE_BRUTE_FORCE` (force O(N) scan; for benchmarking) |

---

## `facet_object`

Represents a single triangular facet. User code accesses facets through `surface%facet_at(i)` rather than directly.

### Public members

| Member | Type | Description |
|---|---|---|
| `normal` | `type(vector_R8P)` | Outward unit normal |
| `vertex(3)` | `type(vector_R8P)` | Three vertices |
| `centroid` | `type(vector_R8P)` | Facet centroid |
| `id` | `integer(I4P)` | Global facet ID (1-based) |
| `fcon_edge(3)` | `integer(I4P)` | Connected-facet ids by edge. Conventions: `EDGE_12 = 1` (v1→v2), `EDGE_23 = 2` (v2→v3), `EDGE_31 = 3` (v3→v1). 0 means disconnected. |
| `bb(2)` | `type(vector_R8P)` | Facet AABB: `bb(1)` = min corner, `bb(2)` = max corner |
| `edge_pnormal(3)`, `vertex_pnormal(3)` | `type(vector_R8P)` | Angle-weighted pseudo-normals for the Bærentzen–Aanæs sign test (populated by `analyze`) |

### Selected methods

| Method | Description |
|---|---|
| `compute_normal` | Recompute the normal from vertex winding |
| `compute_metrix` | Recompute plane-equation coefficients (used by point-to-facet distance) |
| `compute_distance(point, distance)` | Closest squared distance from a point to this facet |
| `compute_distance_with_region(point, distance, closest, region)` | Same plus closest point and Voronoi region tag (face / edge / vertex) |
| `pseudo_normal_for_region(region)` | Pseudo-normal at the given Voronoi region |
| `solid_angle(point)` | Projected solid angle subtended by this facet |
| `do_ray_intersect(origin, direction)` | Ray-triangle intersection test |
| `translate`, `rotate`, `mirror`, `resize` | Per-facet geometry transforms |
| `load_from_file_ascii` / `load_from_file_binary` | Low-level facet I/O |
| `save_into_file_ascii` / `save_into_file_binary` | Low-level facet I/O |
