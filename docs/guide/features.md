---
title: Features
---

# Features

## STL File I/O

- Load ASCII or binary STL files from disk via `surface%load_from_file`
- Automatic format detection (`guess_format=.true.`) — no need to know the file format upfront
- Load with on-the-fly clipping: `clip_min` / `clip_max` bounding box arguments discard facets outside the box during load
- NaN/Inf coordinates are rejected at load time via the optional `status` argument (`STATUS_INVALID_INPUT`); garbage input never enters the geometry pipeline
- Save surfaces to ASCII or binary STL via `surface%save_into_file`

## Surface Analysis

After calling `surface%analyze` (or `surface%sanitize`, which calls it internally), the following data are computed and stored:

- Axis-aligned bounding box (`get_bmin`, `get_bmax`)
- Enclosed signed volume (`get_volume`) — positive for outward-oriented closed bodies (standard divergence-theorem convention) and centroid (`get_centroid`)
- Facet connectivity — symmetric edge-adjacency via the sort-and-pair algorithm (matches trimesh / Open3D / libigl). `f1.fcon_edge(e) = f2` always implies a back-link from `f2` to `f1`
- Non-manifold edge detection (3+ incident facets are explicitly counted, never silently collapsed to a single neighbour)
- Disconnected-edge detection: facets with 1, 2, or 3 disconnected edges are catalogued separately
- AABB acceleration structure (SAH BVH by default; octree available)
- Vertex / edge angle-weighted pseudo-normals (Bærentzen–Aanæs) used by signed-distance queries

### Validity predicates

| Predicate | Returns `.true.` iff |
|---|---|
| `is_watertight()` | Every edge has exactly two incident facets (no boundary, no non-manifold) |
| `is_manifold()`   | No non-manifold edges (boundary edges allowed — open shells can be manifold) |
| `is_volume()`     | `is_watertight()` AND volume > 0 AND finite centroid |

## Surface Manipulation

All manipulation methods operate directly on `surface_stl_object`:

- **Translate** — by a 3D vector (`delta=`) or by scalar components (`x=`, `y=`, `z=`)
- **Rotate** — around an arbitrary axis by angle in radians, or by a given rotation matrix
- **Mirror** — with respect to a plane defined by its normal, or by a given mirror matrix
- **Resize (scale)** — by a 3D vector factor or by scalar per-axis factors; optionally scale about the surface centroid
- **Clip** — discard facets outside an axis-aligned bounding box; the cut-off remainder is returned as a separate surface
- **Merge** — combine two STL surfaces into one (`merge_solids`)

## Mesh repair pipeline

`surface%sanitize` runs the full repair sequence in the right order, with one stderr warning per defect class detected:

1. `remove_degenerate_facets` — drop zero-area / sliver triangles before their NaN normals contaminate downstream queries
2. `connect_nearby_vertices` — union-find vertex deduplication on a tolerance
3. `analyze` — rebuild metrix, connectivity, vertex occurrences
4. `remove_duplicate_facets` — drop literal duplicates (orientation-agnostic; matches the typical CAD-export defect)
5. `analyze` — rebuild after duplicate removal
6. `sanitize_normals` — BFS-propagate winding consistency across the connectivity graph, then flip globally so volume > 0 (outward orientation)

Individual passes (`remove_degenerate_facets`, `remove_duplicate_facets`, `connect_nearby_vertices`, `sanitize_normals`, `reverse_normals`) remain callable for fine-grained control.

## Distance and point-in-polyhedron queries

- `surface%distance(point, is_signed, is_square_root, sign_algorithm)` returns the minimum distance from a point to the surface
  - Default returns squared distance; pass `is_square_root=.true.` for Euclidean
  - Pass `is_signed=.true.` for the signed distance (negative inside, positive outside)
- Sign-determination algorithms (`sign_algorithm` argument):
  - `SIGN_PSEUDO_NORMAL` (default) — Bærentzen–Aanæs pseudo-normal test. Fused with the closest-facet traversal so signed distance costs the same as unsigned. Requires consistently oriented outward normals (which `sanitize` produces).
  - `SIGN_RAY_INTERSECTIONS` — odd/even axis-aligned ray-cast count. Robust on meshes with mixed orientation but ~3× slower than pseudo-normal.
  - `SIGN_SOLID_ANGLE` — sum of projected solid angles ≈ ±4π. Robust but unaccelerated; O(N) per query.
- `surface%is_point_inside(point, sign_algorithm)` is the predicate form of the same test.

All tree-accelerated queries are **bit-exact** vs. brute force on the regression suite — best-first traversal with d² pruning never prunes the true closest facet.

## AABB acceleration

Two tree structures are available; both share the same traversal code and produce identical query results:

- `AABB_TREE_SAH_BVH` **(default)** — binary BVH built top-down by partitioning **triangles** along the longest centroid-bbox axis, choosing splits via the bucketed surface-area heuristic (16 buckets, 3 axes, evaluate at each bucket boundary). Adapts to local triangle density.
- `AABB_TREE_OCTREE` — 8-way space-partitioning octree built to a uniform refinement depth (auto-tuned from facet count via `AABB_AUTO_REFINEMENT`).

Selection is per-surface via the optional `aabb_tree_kind` argument on `load_from_file` / `adopt_facets` / `analyze`. Benchmark numbers (`gfortran -O2`, dragon-fine, 24318 facets, 32³ query points):

| Tree | Build | Query |
|---|---|---|
| Octree | 0.53 s | 27.32 s |
| SAH BVH | 0.32 s | **0.25 s** |

## `fossilizer` CLI

A command-line utility (`src/app/fossilizer.f90`) wraps the library for interactive STL processing:

- Load one or more STL files
- Apply clip, merge, mirror, rotate, translate, resize, sanitize, and connectivity repair operations
- Print per-file statistics to stdout
- Save the result to a named output file

## Compiler Support

| Compiler | Status |
|----------|--------|
| GNU gfortran ≥ 5.3 | Supported |
| Intel Fortran ≥ 16.x | Supported |
| NVIDIA HPC Fortran (nvfortran) | Supported |
| IBM XL Fortran | Not tested |
| g95 | Not tested |
| NAG Fortran | Not tested |

## Design Principles

- **Pure Fortran** — no C extensions, no system calls beyond standard I/O
- **OOP** — all functionality exposed as type-bound procedures; the public API is just `surface_stl_object` and `facet_object`
- **TDD** — every public operation is exercised by automated tests in `src/tests/`. The signed-distance regression test asserts bit-exact agreement between tree-accelerated and brute-force queries
- **KISS** — simple, focused API without unnecessary abstractions
- **Free & Open Source** — multi-licensed for both FOSS and commercial use
