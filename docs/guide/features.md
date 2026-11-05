---
title: Feature map
---

# Feature map

Everything FOSSIL does, grouped as the reference is grouped. Each row names the entry point, the reference page
and, where there is one, the tutorial chapter that uses it. Every public operation is a type-bound procedure of one
of four types, `surface_stl_object`, `curve_object`, `distance_field_object` and `facet_object`, or a free procedure
of a module of its own (named in the row); `use fossil` re-exports the types and the named constants.

## Surfaces

| Capability | Entry point | Reference | Tutorial |
|---|---|---|---|
| Load ASCII or binary STL, format detected, NaN/Inf rejected, optional clip box | `surface%load_from_file` | [I/O, analysis, manipulation](/guide/usage) | [1](/manual/tutorial/01-first-surface) |
| Save ASCII or binary STL | `surface%save_into_file` | [I/O, analysis, manipulation](/guide/usage) | [1](/manual/tutorial/01-first-surface) |
| Bounding box, volume, centroid, connectivity, shells and bodies, pseudo-normals, AABB tree | `surface%analyze` (run by load, sanitize, clip) | [I/O, analysis, manipulation](/guide/usage) | [1](/manual/tutorial/01-first-surface) |
| Validity: watertight, manifold, closed volume | `is_watertight`, `is_manifold`, `is_volume` | [I/O, analysis, manipulation](/guide/usage) | [2](/manual/tutorial/02-repair) |
| Repair: slivers, duplicates, gaps, orientation, one warning per defect class | `surface%sanitize` | [I/O, analysis, manipulation](/guide/usage) | [2](/manual/tutorial/02-repair) |
| Hole filling (Liepa minimum-weight patches, refinement, thin-plate fairing) | `surface%fill_holes` | [Hole filling](/guide/advanced/hole-filling) | [2](/manual/tutorial/02-repair) |
| Translate, rotate, mirror, resize, clip, merge | `translate`, `rotate`, `mirror`, `resize`, `clip`, `merge_solids` | [I/O, analysis, manipulation](/guide/usage) | [3](/manual/tutorial/03-manipulation) |
| Self-intersection detection and resolution (Möller tri-tri, exact predicates) | `find_self_intersections`, `resolve_self_intersections` | [Self-intersection](/guide/advanced/self-intersection) | |
| Booleans: union, intersection, difference, symmetric difference (Zhou et al. 2016) | `surface%boolean` | [Boolean operations](/guide/advanced/booleans) | |
| Quadric decimation (Garland and Heckbert 1997) | `surface%decimate` | [Mesh decimation](/guide/advanced/decimate) | |
| Isotropic remeshing, optional feature preservation (Botsch and Kobbelt 2004) | `surface%isotropic_remesh` | [Isotropic remeshing](/guide/advanced/isotropic-remesh) | |
| Explicit and Taubin smoothing | `surface%smooth` | [Mesh smoothing](/guide/advanced/smoothing) | |
| Alpha wrapping: a watertight surrogate of any triangle soup (Portaneri et al. 2022) | `surface%alpha_wrap` | [Alpha wrapping](/guide/advanced/alpha-wrap) | |
| Isosurface extraction (Lorensen and Cline 1987) | `extract_isosurface` | [Marching cubes](/guide/advanced/marching-cubes) | |
| Shape-diameter segmentation (Shapira et al. 2008) | `surface%segment_sdf` | [SDF segmentation](/guide/advanced/sdf-segmentation) | |
| Cotangent Laplacian and barycentric mass | `surface%cotangent_laplacian` | [Cotangent Laplacian](/guide/advanced/cotangent-laplacian) | |
| Gaussian and mean curvature per vertex; pointwise O(h²) by a cubic jet fit | `gaussian_curvature`, `mean_curvature`, `jet_curvature` | [Per-vertex curvature](/guide/advanced/curvature) | |
| Mass, centre of mass, inertia tensor, principal axes (exact polyhedral integrals); 2-D second moments | `surface%mass_properties`, `curve%mass_properties`, `inertia_to_world`, `principal_frame` | [Mass properties](/guide/advanced/mass-properties) | |
| Surface deviation: one-sided and symmetric Hausdorff distance, mean, rms, per-vertex distance | `surface%deviation`, `surface%hausdorff` | [Surface deviation](/guide/advanced/surface-deviation) | |

## Distance and signs

| Capability | Entry point | Reference | Tutorial |
|---|---|---|---|
| Signed or unsigned distance of a point, bit-exact against brute force | `surface%distance` | [I/O, analysis, manipulation](/guide/usage) | [4](/manual/tutorial/04-distance) |
| Closest facet, edge, vertex, point and normal | `surface%compute_distance` | [I/O, analysis, manipulation](/guide/usage) | [4](/manual/tutorial/04-distance) |
| Inside or outside | `surface%is_point_inside` | [Sign algorithms benchmark](/guide/advanced/sign-benchmark) | [4](/manual/tutorial/04-distance) |
| Four sign algorithms: pseudo-normal (default), ray parity, solid angle, fast winding number | `sign_algorithm=` | [Sign algorithms benchmark](/guide/advanced/sign-benchmark) | [4](/manual/tutorial/04-distance) |
| Generalized and fast winding number (Jacobson 2013, Barill 2018) | `surface%winding_number` | [Winding number](/guide/advanced/winding-number) | |
| Batched distances, OpenMP-parallel and bit-exact with serial | `surface%distance_many` | [I/O, analysis, manipulation](/guide/usage) | [5](/manual/tutorial/05-batched) |
| SAH BVH (default) or octree acceleration | `aabb_tree_kind=` | [I/O, analysis, manipulation](/guide/usage) | [5](/manual/tutorial/05-batched) |
| Ray queries: first hit, all hits, any hit | `intersect_ray_first`, `_all`, `_any` | [Ray queries](/guide/advanced/ray-queries) | |
| Signed-distance field on a stretched grid: exact band, flood-fill signs, fast sweeping | `distance_field_object%compute` | [Signed-distance fields](/guide/advanced/distance-field) | [6](/manual/tutorial/06-distance-field) |
| Eikonal diagnostic of a field | `distance_field_object%gradient_diagnostic` | [Signed-distance fields](/guide/advanced/distance-field) | |

## Immersed boundaries

| Capability | Entry point | Reference | Tutorial |
|---|---|---|---|
| 2-D curves: closed and open polylines, segment BVH, signed distance, 2-D winding number, repair | `curve_object` | [2-D geometry](/guide/advanced/2d-geometry) | [7](/manual/tutorial/07-curves) |
| Axis-aligned section of an STL as a curve | `slice_surface` | [2-D geometry](/guide/advanced/2d-geometry) | [7](/manual/tutorial/07-curves) |
| Rich queries: closest point, feature, pseudo-normal, body; ambiguous feet; wall curvature | `query_point`, `query_many` | [Immersed-boundary queries](/guide/advanced/ib-queries) | [8](/manual/tutorial/08-ib-queries) |
| Image points and probes with a visibility test | `image_points`, `probe_many` | [Immersed-boundary queries](/guide/advanced/ib-queries) | [8](/manual/tutorial/08-ib-queries) |
| Wall distance with body id and early cutoff | `wall_distance_many` | [Immersed-boundary queries](/guide/advanced/ib-queries) | |
| Grid-adapted markers | `surface_markers`, `curve_markers` | [Immersed-boundary queries](/guide/advanced/ib-queries) | |
| AMR block and cell classification: cut, near, inside, outside | `box_classify_many` | [AMR tagging](/guide/advanced/amr-tagging) | [9](/manual/tutorial/09-amr-tagging) |
| Signed tag bands on a block | `tag_block` | [AMR tagging](/guide/advanced/amr-tagging) | [9](/manual/tutorial/09-amr-tagging) |
| Exact grid-line crossings and wall fractions | `line_crossings`, `cell_wall_fractions` | [Cut-cell geometry](/guide/advanced/cut-cells) | |
| 2-D cut cells: fractions, apertures, centroids, wetted boundary, pieces | `curve%cut_cells` | [Cut-cell geometry](/guide/advanced/cut-cells) | [10](/manual/tutorial/10-cut-cells) |
| 3-D cut cells: face fractions, volumes, fluid pieces with their own measures | `cut_cells_3d` (module `fossil_cut_cells_3d`) | [Cut-cell geometry](/guide/advanced/cut-cells) | [10](/manual/tutorial/10-cut-cells) |
| Rigid poses, body-frame queries, O(N) tree refit | `rigid_pose_t`, `pose=`, `move` | [Moving bodies](/guide/advanced/moving-bodies) | [11](/manual/tutorial/11-moving-bodies) |
| Swept geometry of a time step: exact swept volumes, time-averaged apertures | `swept_cells_2d`, `swept_cells_3d` | [Moving bodies](/guide/advanced/moving-bodies) | [11](/manual/tutorial/11-moving-bodies) |

## GPU

| Capability | Entry point | Reference | Tutorial |
|---|---|---|---|
| Distance, rich queries, wall distance, box classification, posed queries on the device (OpenACC), bit-exact against the CPU path | `enter_device`, `*_device`, `exit_device` | [GPU queries](/guide/advanced/device-distance) | [12](/manual/tutorial/12-gpu) |

## `fossilizer` CLI

`src/app/fossilizer.f90` wraps the library for the command line: load one or more STL files, clip, merge, mirror,
rotate, translate, resize, sanitize and repair connectivity, print statistics and, with `--mass_properties
[--density ρ]`, the mass properties of a closed body; save the result.

## Compilers

| Compiler | Status |
|---|---|
| GNU gfortran 14.3+ (14 series) or ≥ 15.1 | Supported (13.x and 14.0–14.2 untested, see [Installation](/guide/installation)) |
| Intel Fortran ≥ 16.x | Supported |
| NVIDIA HPC Fortran (nvfortran) ≥ 26.9 | The device build of the GPU queries, compiled without FMA contraction |

## Design principles

- **Pure Fortran**: no C extensions, no system calls beyond standard I/O.
- **Exact where it decides**: orientation predicates are exact (`orient2d`, `orient3d`); tree-accelerated queries
  equal brute force bit for bit; the GPU path equals the CPU path bit for bit.
- **OOP and TDD**: every operation is a type-bound procedure, every public operation is exercised by a test in
  `src/tests/`.
- **Free and open source**: multi-licensed for FOSS and commercial use.
