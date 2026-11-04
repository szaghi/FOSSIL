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
- **Rotate** — around an arbitrary axis by angle in radians, or by a given rotation matrix (an improper matrix keeps the outward orientation, as mirror does)
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
  - `SIGN_RAY_INTERSECTIONS` — odd/even crossing count along three generic (non-axis-aligned) rays, majority vote. Robust on meshes with mixed orientation but ~3× slower than pseudo-normal.
  - `SIGN_SOLID_ANGLE` — sum of projected solid angles ≈ ±4π. Robust but unaccelerated; O(N) per query.
  - `SIGN_FAST_WINDING_NUMBER` — inside when the fast (hierarchical, second-order) winding number exceeds 0.5. The most robust choice on surfaces with holes or gaps; about 13 µs per query on a 70k-facet mesh on top of the distance search.
- `surface%is_point_inside(point, sign_algorithm)` is the predicate form of the same test.
- Measured at sharp features (issue #39): on watertight meshes all four are right down to 1e-12 of the body size (the pseudo-normal down to ~1e-9 on needle triangles); on dirty meshes the fast winding number is the most robust. See [Sign algorithms at sharp features](/guide/advanced/sign-benchmark).
- `surface%distance_many(points, distances, ...)` is the batch form: one distance per point of an array, with the same options. The loop is OpenMP-parallel when built with OpenMP (`tests-gnu-openmp`) and serial otherwise, with identical results.
- `surface%distance_many_device(points, distances, ...)` runs the same batch on the GPU (OpenACC, nvfortran), between `enter_device` and `exit_device`. SAH BVH only; `SIGN_PSEUDO_NORMAL` or `SIGN_FAST_WINDING_NUMBER`; bit-exact vs the CPU path where validated. See [GPU distance queries](/guide/advanced/device-distance).
- `distance_field_object` computes the signed distance of every cell of a rectilinear, possibly stretched grid: exact distances in a narrow band, signs by flood fill, the far field by fast sweeping (OpenMP, bit-exact with serial; band and point queries optionally on the GPU). `fossil_ib_generator` uses it. See [Signed-distance field](/guide/advanced/distance-field).
- `curve_object` is the 2-D geometry for quadtree solvers: closed (or open) polylines in an axis-aligned plane with a segment BVH, signed distance, a 2-D fast winding number, exact self-intersection and nesting tests, repair, and batched queries on the CPU or the device. `slice_surface` cuts an STL by an axis-aligned plane into a curve, and `distance_field_object` builds 2-D fields from curves. See [2-D geometry](/guide/advanced/2d-geometry).
- `query_many` / `query_many_device` return, per point, the closest point, the closest facet (segment), the feature, the pseudo-normal and the body; `probe_many` and `image_points` give probe and image points with a visibility test; `query_many(ambiguous=)` flags points with several foot points; `distance_field_object%compute(store_closest=.true.)` stores closest point, id and normal on the band (vector level set). Shells and bodies come from `analyze`. 3-D and 2-D. See [Immersed-boundary queries](/guide/advanced/ib-queries).
- `wall_distance_many` (CPU/GPU) gives the wall distance with the body id and an early cutoff; `distance_field_object%gradient_diagnostic` reports q = |1 - |grad phi|| per node; `query_many(curvatures=)` gives the wall curvature at the foot point; `surface_markers` / `curve_markers` give grid-adapted markers (spacing near Δx, with areas and normals). See [Immersed-boundary queries](/guide/advanced/ib-queries).
- Rigidly moving bodies: `rigid_pose_t` poses, an optional `pose=` on the point queries (body-frame queries, no rebuild, CPU and device), `move(pose)` with an O(N) tree refit (queries equal a rebuild bit for bit), `state_changes` for fresh and dead cells, and `swept_cells_2d` / `swept_cells_3d` for the swept geometry of a step (exact swept volumes, time-averaged face apertures, projected walls; exact for translations). See [Moving bodies](/guide/advanced/moving-bodies).
- `line_crossings` / `cell_wall_fractions` give the exact crossings of grid lines with the wall (each counted once, also through edges and vertices) and the wall fractions at cell centres; `curve%cut_cells` gives 2-D cut cells (fractions, apertures, centroids, wetted boundary, pieces) and `cut_cells_3d` 3-D face fractions, cell volumes and fluid pieces. See [Cut-cell geometry](/guide/advanced/cut-cells).
- `box_classify_many` / `box_classify_many_device` classify AMR blocks or cells as cut, near (within a radius), inside or outside from the BVH alone (3-D conservative separating-axis test, 2-D exact); `tag_block` gives signed tag bands (Chebyshev cell distance to the nearest cut cell) on a stretched block. 3-D and 2-D. See [AMR tagging](/guide/advanced/amr-tagging).

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

## Advanced features

Beyond the core load / sanitize / distance API above, FOSSIL ships
**thirteen** algorithmic primitives drawn from the libigl / CGAL geometry-
processing toolkit. Each is a single TBP on `surface_stl_object`.
Conceptual pages (what the algorithm is, when to use it, when not to,
code examples, known limitations) live under
[Advanced Features](/guide/advanced/).

| Feature | API entry point | Page |
|---|---|---|
| Boolean operations (∪ ∩ ∖ △) | `surface%boolean(other, op, status)` | [Booleans](/guide/advanced/booleans) |
| Self-intersection detection | `surface%find_self_intersections(pairs, status)` | [Self-intersection](/guide/advanced/self-intersection) |
| Self-intersection resolution | `surface%resolve_self_intersections(status)` | [Self-intersection](/guide/advanced/self-intersection) |
| Mesh decimation (QEM) | `surface%decimate(target_facets, status)` | [Decimation](/guide/advanced/decimate) |
| Generalized winding number | `surface%winding_number(point, beta)` | [Winding number](/guide/advanced/winding-number) |
| Marching cubes | `extract_isosurface(values, bmin, bmax, iso, ...)` | [Marching cubes](/guide/advanced/marching-cubes) |
| Alpha wrapping | `surface%alpha_wrap(alpha, offset, ...)` | [Alpha wrap](/guide/advanced/alpha-wrap) |
| Isotropic remeshing | `surface%isotropic_remesh(target_length, ...)` | [Isotropic remesh](/guide/advanced/isotropic-remesh) |
| SDF segmentation | `surface%segment_sdf(facet_labels, ...)` | [SDF segmentation](/guide/advanced/sdf-segmentation) |
| Ray-mesh queries | `surface%intersect_ray_first/_all/_any(...)` | [Ray queries](/guide/advanced/ray-queries) |
| Cotangent Laplacian + barycentric mass | `surface%cotangent_laplacian(L, M, status)` | [Cotangent Laplacian](/guide/advanced/cotangent-laplacian) |
| Per-vertex Gaussian curvature | `surface%gaussian_curvature(K, status)` | [Curvature](/guide/advanced/curvature) |
| Per-vertex signed mean curvature | `surface%mean_curvature(H, status)` | [Curvature](/guide/advanced/curvature) |
| Per-vertex H and K by a cubic jet fit (pointwise O(h²)) | `surface%jet_curvature(H, K, status)` | [Curvature](/guide/advanced/curvature) |
| Mesh smoothing (explicit + Taubin) | `surface%smooth(method, lambda, mu, iterations, status)` | [Smoothing](/guide/advanced/smoothing) |
| Hole filling (Liepa minimum-weight patches, optional refinement and thin-plate fairing) | `surface%fill_holes(max_hole_edges, max_hole_size, ..., refine, fair)` | [Hole filling](/guide/advanced/hole-filling) |

RANSAC shape detection was deliberately dropped from the roadmap; a few
solver-dependent items (heat-method geodesics, implicit smoothing)
remain deferred — see the [implementation
issue](https://github.com/szaghi/FOSSIL/issues/18) for the full
rationale.

## `fossilizer` CLI

A command-line utility (`src/app/fossilizer.f90`) wraps the library for interactive STL processing:

- Load one or more STL files
- Apply clip, merge, mirror, rotate, translate, resize, sanitize, and connectivity repair operations
- Print per-file statistics to stdout
- Save the result to a named output file

## Compiler Support

| Compiler | Status |
|----------|--------|
| GNU gfortran 14.3+ (14 series) or ≥ 15.1 | Supported (13.x and 14.0–14.2 untested, see [Installation](installation.md)) |
| Intel Fortran ≥ 16.x | Supported |
| NVIDIA HPC Fortran (nvfortran) ≥ 26.9 | GPU device build of the BVH distance path only |
| IBM XL Fortran | Not tested |
| g95 | Not tested |
| NAG Fortran | Not tested |

## Design Principles

- **Pure Fortran** — no C extensions, no system calls beyond standard I/O
- **OOP** — all functionality exposed as type-bound procedures; the public API is just `surface_stl_object` and `facet_object`
- **TDD** — every public operation is exercised by automated tests in `src/tests/`. The signed-distance regression test asserts bit-exact agreement between tree-accelerated and brute-force queries
- **KISS** — simple, focused API without unnecessary abstractions
- **Free & Open Source** — multi-licensed for both FOSS and commercial use
