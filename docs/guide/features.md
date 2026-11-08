---
title: Feature map
---

# Feature map

Everything FOSSIL does, grouped as the reference is grouped. Each row names the entry point, the reference page
and, where there is one, the tutorial chapter that uses it. Every public operation is a type-bound procedure of one
of four types, `surface_stl_object`, `curve_object`, `distance_field_object` and `facet_object`, or a free procedure
of a module of its own (named in the row); `use fossil` re-exports the types and the named constants.

## Surfaces

| Capability | Entry point | Reference | Tutorial | Cookbook |
|---|---|---|---|---|
| Load ASCII or binary STL, format detected, NaN/Inf rejected, optional clip box | `surface%load_from_file` | [I/O, analysis, manipulation](/guide/usage) | [1](/manual/tutorial/01-first-surface) | [1](/manual/cookbook#read-an-stl-file-ascii-or-binary), [2](/manual/cookbook#reject-a-file-with-nan-or-inf-coordinates), [3](/manual/cookbook#load-only-a-region-of-a-large-file), [4](/manual/cookbook#load-several-files-as-one-surface) |
| Save ASCII or binary STL | `surface%save_into_file` | [I/O, analysis, manipulation](/guide/usage) | [1](/manual/tutorial/01-first-surface) | [1](/manual/cookbook#save-as-ascii-or-binary) |
| Bounding box, volume, centroid, connectivity, shells and bodies, pseudo-normals, AABB tree | `surface%analyze` (run by load, sanitize, clip) | [I/O, analysis, manipulation](/guide/usage) | [1](/manual/tutorial/01-first-surface) | [1](/manual/cookbook#check-that-a-surface-is-valid) |
| Validity: watertight, manifold, closed volume | `is_watertight`, `is_manifold`, `is_volume` | [I/O, analysis, manipulation](/guide/usage) | [2](/manual/tutorial/02-repair) | [1](/manual/cookbook#check-that-a-surface-is-valid) |
| Repair: slivers, duplicates, gaps, orientation, one warning per defect class | `surface%sanitize` | [I/O, analysis, manipulation](/guide/usage) | [2](/manual/tutorial/02-repair) | [1](/manual/cookbook#fix-an-inside-out-surface) |
| Hole filling (Liepa minimum-weight patches, refinement, thin-plate fairing) | `surface%fill_holes` | [Hole filling](/guide/advanced/hole-filling) | [2](/manual/tutorial/02-repair) | [1](/manual/cookbook#close-the-holes-of-a-surface) |
| Translate, rotate, mirror, resize, clip, merge | `translate`, `rotate`, `mirror`, `resize`, `clip`, `merge_solids` | [I/O, analysis, manipulation](/guide/usage) | [3](/manual/tutorial/03-manipulation) | [1](/manual/cookbook#load-only-a-region-of-a-large-file), [2](/manual/cookbook#rotate-a-body), [3](/manual/cookbook#split-a-surface-by-a-box) |
| Self-intersection detection and resolution (exact `pair_cut`, exact predicates) | `find_self_intersections`, `resolve_self_intersections` | [Self-intersection](/guide/advanced/self-intersection) | [12](/manual/tutorial/12-geometry) | [1](/manual/cookbook#find-self-intersections-and-remove-them), [2](/manual/cookbook#resolve-the-self-intersections-of-one-surface) |
| Booleans: union, intersection, difference, symmetric difference (Zhou et al. 2016) | `surface%boolean` | [Boolean operations](/guide/advanced/booleans) | [12](/manual/tutorial/12-geometry) | [1](/manual/cookbook#unite-two-bodies) |
| Quadric decimation (Garland and Heckbert 1997) | `surface%decimate` | [Mesh decimation](/guide/advanced/decimate) | | [1](/manual/cookbook#decimate-a-mesh), [2](/manual/cookbook#decimate-to-a-tolerance) |
| Isotropic remeshing, optional feature preservation (Botsch and Kobbelt 2004) | `surface%isotropic_remesh` | [Isotropic remeshing](/guide/advanced/isotropic-remesh) | [12](/manual/tutorial/12-geometry) | [1](/manual/cookbook#remesh-to-a-uniform-edge-length) |
| Explicit and Taubin smoothing | `surface%smooth` | [Mesh smoothing](/guide/advanced/smoothing) | | [1](/manual/cookbook#smooth-a-noisy-surface) |
| Alpha wrapping: a watertight surrogate of any triangle soup (Portaneri et al. 2022) | `surface%alpha_wrap` | [Alpha wrapping](/guide/advanced/alpha-wrap) | | [1](/manual/cookbook#wrap-a-triangle-soup-into-a-closed-surface) |
| Isosurface extraction (Lorensen and Cline 1987) | `extract_isosurface` | [Marching cubes](/guide/advanced/marching-cubes) | | [1](/manual/cookbook#extract-an-isosurface-from-a-sampled-field) |
| Shape-diameter segmentation (Shapira et al. 2008) | `surface%segment_sdf` | [SDF segmentation](/guide/advanced/sdf-segmentation) | | [1](/manual/cookbook#segment-a-part-by-its-shape-diameter) |
| Cotangent Laplacian and barycentric mass | `surface%cotangent_laplacian` | [Cotangent Laplacian](/guide/advanced/cotangent-laplacian) |  | [1](/manual/cookbook#assemble-the-cotangent-laplacian-of-a-surface) |
| Gaussian and mean curvature per vertex; pointwise O(h²) by a cubic jet fit | `gaussian_curvature`, `mean_curvature`, `jet_curvature` | [Per-vertex curvature](/guide/advanced/curvature) | | [1](/manual/cookbook#measure-the-curvature-at-a-point), [2](/manual/cookbook#curvature-at-every-vertex) |
| Mass, centre of mass, inertia tensor, principal axes (exact polyhedral integrals); 2-D second moments | `surface%mass_properties`, `curve%mass_properties`, `inertia_to_world`, `principal_frame` | [Mass properties](/guide/advanced/mass-properties) | [12](/manual/tutorial/12-geometry) | [1](/manual/cookbook#mass-and-inertia-of-a-body), [2](/manual/cookbook#second-moments-of-area-of-a-section), [3](/manual/cookbook#inertia-of-a-turned-body) |
| Surface deviation: one-sided and symmetric Hausdorff distance, mean, rms, per-vertex distance | `surface%deviation`, `surface%hausdorff`, `curve%deviation`, `curve%hausdorff` | [Surface deviation](/guide/advanced/surface-deviation) | [12](/manual/tutorial/12-geometry) | [1](/manual/cookbook#measure-how-far-a-repair-moved-the-surface), [2](/manual/cookbook#compare-two-versions-of-a-part) |
| Feature edges: sharp, boundary and non-manifold edges, corners, polylines; refinement and distances along them; 2-D corners | `surface%feature_edges`, `build_features`, `box_feature_many`, `tag_block(feature_tags=)`, `query_many(feature_distances=)`, `curve%feature_corners` | [Feature edges](/guide/advanced/feature-edges) | [12](/manual/tutorial/12-geometry) | [1](/manual/cookbook#find-the-sharp-edges-of-a-part), [2](/manual/cookbook#refine-the-cells-along-sharp-edges), [3](/manual/cookbook#find-the-corners-of-a-section) |
| Local thickness and gaps per facet: ray (shape diameter) and medial (shrinking ball); thin regions under a cell size; refinement of thin walls | `surface%thickness`, `thin_regions`, `build_thin`, `box_thin_many`, `tag_block(thin_tags=)`, `curve%thickness` | [Thickness and gaps](/guide/advanced/thickness) | [12](/manual/tutorial/12-geometry) | [1](/manual/cookbook#find-walls-too-thin-for-my-grid), [2](/manual/cookbook#refine-the-cells-across-thin-walls) |
| Named patches: one per `solid` block or per file, kept on every facet; multi-solid ASCII output | `load_from_files`, `get_patches_number`, `patch_name`, `facet_patch`, `set_patch` | [Named patches](/guide/advanced/patches) | [12](/manual/tutorial/12-geometry) | [1](/manual/cookbook#wall-inlet-and-outlet-from-one-stl), [2](/manual/cookbook#load-several-files-as-one-surface) |
| Distance between two bodies: closest facets and points, contact; bit-identical to brute force | `surface%distance_to` | [Body–body distance](/guide/advanced/body-distance) | [14](/manual/tutorial/14-contact) | [1](/manual/cookbook#how-close-are-my-two-bodies) |
| Time of first contact of moving bodies over a step: certified (never late, no tunnelling), several bodies, 2-D twins | `surface%time_of_impact`, `curve%time_of_impact`, `first_impact`, `first_impact_2d` | [Continuous collision detection](/guide/advanced/collision) | [14](/manual/tutorial/14-contact) | [1](/manual/cookbook#when-do-my-moving-bodies-first-touch), [2](/manual/cookbook#which-of-my-moving-bodies-touch-first) |

## Distance and signs

| Capability | Entry point | Reference | Tutorial | Cookbook |
|---|---|---|---|---|
| Signed or unsigned distance of a point, bit-exact against brute force | `surface%distance` | [I/O, analysis, manipulation](/guide/usage) | [4](/manual/tutorial/04-distance) | [1](/manual/cookbook#tell-whether-a-point-is-inside) |
| Closest facet, edge, vertex, point and normal | `surface%compute_distance` | [I/O, analysis, manipulation](/guide/usage) | [4](/manual/tutorial/04-distance) | [1](/manual/cookbook#find-the-closest-point-and-its-normal) |
| Inside or outside | `surface%is_point_inside` | [Sign algorithms benchmark](/guide/advanced/sign-benchmark) | [4](/manual/tutorial/04-distance) | [1](/manual/cookbook#tell-whether-a-point-is-inside) |
| Four sign algorithms: pseudo-normal (default), ray parity, solid angle, fast winding number | `sign_algorithm=` | [Sign algorithms benchmark](/guide/advanced/sign-benchmark) | [4](/manual/tutorial/04-distance) | [1](/manual/cookbook#get-a-reliable-sign-on-a-mesh-with-holes) |
| Generalized and fast winding number (Jacobson 2013, Barill 2018) | `surface%winding_number` | [Winding number](/guide/advanced/winding-number) | | [1](/manual/cookbook#measure-how-much-of-a-surface-surrounds-a-point) |
| Batched distances, OpenMP-parallel and bit-exact with serial | `surface%distance_many` | [I/O, analysis, manipulation](/guide/usage) | [5](/manual/tutorial/05-batched) | [1](/manual/cookbook#compute-the-signed-distance-on-a-whole-grid) |
| SAH BVH (default) or octree acceleration | `aabb_tree_kind=` | [I/O, analysis, manipulation](/guide/usage) | [5](/manual/tutorial/05-batched) | [1](/manual/cookbook#choose-the-acceleration-tree) |
| Read-only access to the acceleration tree: node boxes down to a depth, leaves and their facets (both kinds; curve twin) | `surface%tree_boxes`, `curve%tree_boxes` | [I/O, analysis, manipulation](/guide/usage#reading-the-tree) | [5](/manual/tutorial/05-batched) | [1](/manual/cookbook#read-the-acceleration-tree) |
| Ray queries: first hit, all hits (one per crossing, exact), any hit | `intersect_ray_first`, `_all`, `_any` | [Ray queries](/guide/advanced/ray-queries) | | [1](/manual/cookbook#cast-a-ray) |
| Signed-distance field on a stretched grid: exact band, flood-fill signs, fast sweeping | `distance_field_object%compute` | [Signed-distance fields](/guide/advanced/distance-field) | [6](/manual/tutorial/06-distance-field) | [1](/manual/cookbook#compute-the-signed-distance-on-a-whole-grid) |
| Eikonal diagnostic of a field | `distance_field_object%gradient_diagnostic` | [Signed-distance fields](/guide/advanced/distance-field) | | [1](/manual/cookbook#check-a-distance-field-by-the-eikonal-diagnostic) |

## Immersed boundaries

| Capability | Entry point | Reference | Tutorial | Cookbook |
|---|---|---|---|---|
| 2-D curves: closed and open polylines, segment BVH, signed distance, 2-D winding number, repair | `curve_object` | [2-D geometry](/guide/advanced/2d-geometry) | [7](/manual/tutorial/07-curves) | [1](/manual/cookbook#compute-the-signed-distance-to-a-polygon) |
| Axis-aligned section of an STL as a curve | `slice_surface` | [2-D geometry](/guide/advanced/2d-geometry) | [7](/manual/tutorial/07-curves) | [1](/manual/cookbook#cut-a-2-d-section-of-an-stl) |
| Rich queries: closest point, feature, pseudo-normal, body; ambiguous feet; wall curvature | `query_point`, `query_many` | [Immersed-boundary queries](/guide/advanced/ib-queries) | [8](/manual/tutorial/08-ib-queries) | [1](/manual/cookbook#find-the-closest-point-and-its-normal), [2](/manual/cookbook#tell-whether-the-closest-point-is-on-a-face-an-edge-or-a-vertex) |
| Image points and probes with a visibility test | `image_points`, `probe_many` | [Immersed-boundary queries](/guide/advanced/ib-queries) | [8](/manual/tutorial/08-ib-queries) | [1](/manual/cookbook#image-points-and-probes-for-a-ghost-cell-method) |
| Wall distance with body id and early cutoff | `wall_distance_many` | [Immersed-boundary queries](/guide/advanced/ib-queries) | | [1](/manual/cookbook#compute-the-wall-distance-up-to-a-cutoff) |
| Grid-adapted markers | `surface_markers`, `curve_markers` | [Immersed-boundary queries](/guide/advanced/ib-queries) | | [1](/manual/cookbook#place-lagrangian-markers-for-an-ib-method), [2](/manual/cookbook#place-markers-on-a-curve) |
| AMR block and cell classification: cut, near, inside, outside | `box_classify_many` | [AMR tagging](/guide/advanced/amr-tagging) | [9](/manual/tutorial/09-amr-tagging) | [1](/manual/cookbook#tag-the-cells-near-the-wall) |
| Signed tag bands on a block | `tag_block` | [AMR tagging](/guide/advanced/amr-tagging) | [9](/manual/tutorial/09-amr-tagging) | [1](/manual/cookbook#tag-refinement-bands-on-a-block) |
| Exact grid-line crossings and wall fractions | `line_crossings`, `cell_wall_fractions` | [Cut-cell geometry](/guide/advanced/cut-cells) |  | [1](/manual/cookbook#find-where-a-grid-line-crosses-a-curve), [2](/manual/cookbook#find-where-the-wall-cuts-the-links-between-cell-centres) |
| 2-D cut cells: fractions, apertures, centroids, wetted boundary, pieces | `curve%cut_cells` | [Cut-cell geometry](/guide/advanced/cut-cells) | [10](/manual/tutorial/10-cut-cells) | [1](/manual/cookbook#find-where-the-wall-cuts-the-links-between-cell-centres), [2](/manual/cookbook#moments-and-quadrature-of-cut-cells) |
| 3-D cut cells: face fractions, volumes, fluid pieces with their own measures | `cut_cells_3d` (module `fossil_cut_cells_3d`) | [Cut-cell geometry](/guide/advanced/cut-cells) | [10](/manual/tutorial/10-cut-cells) | [1](/manual/cookbook#cut-a-3-d-grid-by-a-body) |
| Cut-cell solver support: piece links, state- and flux-redistribution neighbourhoods, merge targets | `cut_links_2d`, `cut_cells_3d(links=)`, `srd_neighbourhoods`, `frd_neighbourhoods`, `merge_targets` | [Cut cells: solver support](/guide/advanced/cut-solver) | | [1](/manual/cookbook#neighbourhoods-for-a-small-cell-cure), [2](/manual/cookbook#flux-redistribution-weights-of-a-small-cell) |
| Cut-cell moments up to order 4, quadrature rules, edge fractions, exact coarsening | `cut_cells_3d(moments_order=, quadrature=, edges=)`, `curve%cut_cells(moments_order=, quadrature=)`, `cut_coarsen_3d`, `cut_coarsen_2d` | [Cut cells: solver support](/guide/advanced/cut-solver) | | [1](/manual/cookbook#moments-and-quadrature-of-cut-cells), [2](/manual/cookbook#coarsen-cut-cells), [3](/manual/cookbook#coarsen-3-d-cut-cells) |
| Curvilinear structured blocks: topology records (seams, cuts, poles, periodic, interfaces) and their detection | `curvilinear_block_t`, `curvilinear_block_2d_t`, `detect_topology` | [Curvilinear blocks](/guide/advanced/curvilinear-blocks) | [13](/manual/tutorial/13-blocks-overset) | [1](/manual/cookbook#find-how-blocks-connect), [2](/manual/cookbook#find-the-seam-and-the-poles-of-a-3-d-block) |
| Exact cell geometry, certified cell validity, quality metrics | `geometry`, `validity`, `quality` | [Curvilinear blocks](/guide/advanced/curvilinear-blocks) | [13](/manual/tutorial/13-blocks-overset) | [1](/manual/cookbook#check-a-block-s-validity-and-quality) |
| Point location in curvilinear blocks (robust, walking, GPU) and Lagrange stencils of width q | `locate_many`, `locate_many_device`, `locate_in_blocks`, `stencils` | [Curvilinear blocks](/guide/advanced/curvilinear-blocks) | [13](/manual/tutorial/13-blocks-overset) | [1](/manual/cookbook#locate-a-point-in-a-curvilinear-block), [2](/manual/cookbook#find-which-block-holds-a-point), [3](/manual/cookbook#locate-many-points-on-the-gpu) |
| Overset cutters and hole cutting: wall cutters with trimming, exact hole flags (CPU and GPU), hole maps, X-rays | `wall_cutter`, `hole_flags`, `hole_map`, `xray_flags` | [Overset assembly](/guide/advanced/overset) | [13](/manual/tutorial/13-blocks-overset) | [1](/manual/cookbook#cut-holes-with-a-hole-map-or-x-rays), [2](/manual/cookbook#cut-holes-around-a-3-d-body) |
| Overset donors and classification: candidates with quality (CPU and GPU), pair matrix, explicit and implicit selection, fringes, orphans by cause, overlap width, wall projection, moving grids, connectivity output | `donor_candidates`, `classify`, `fringe_levels`, `wall_projection`, `exposed_points`, `build_connectivity` | [Overset assembly](/guide/advanced/overset) | [13](/manual/tutorial/13-blocks-overset) | [1](/manual/cookbook#assemble-an-overset-grid), [2](/manual/cookbook#write-the-overset-connectivity-for-a-flow-solver), [3](/manual/cookbook#assemble-a-3-d-overset-grid), [4](/manual/cookbook#find-donors-next-to-a-curved-wall-that-two-blocks-discretize-differently), [5](/manual/cookbook#find-the-points-a-moving-body-uncovers) |
| Body-fitted grid support: constrained projection onto bodies (patches, feature and patch-boundary curves, corners), block faces, wall-spacing and quality report (y1, y1+, growth), 2-D hyperbolic extrusion | `projection_target`, `project_points`, `project_block_face`, `block_report`, `hyperbolic_grid_2d` | [Body-fitted grid support](/guide/advanced/body-fitted) | [15](/manual/tutorial/15-body-fitted) | [1](/manual/cookbook#put-a-block-face-on-a-body-with-sharp-edges), [2](/manual/cookbook#extrude-an-o-grid-and-check-its-wall-spacing), [3](/manual/cookbook#save-the-wall-spacing-report-of-a-grid), [4](/manual/cookbook#snap-points-onto-a-body-with-sharp-edges), [5](/manual/cookbook#check-the-first-cell-height-of-a-3-d-wall-block) |
| Rigid poses, body-frame queries, O(N) tree refit | `rigid_pose_t`, `pose=`, `move` | [Moving bodies](/guide/advanced/moving-bodies) | [11](/manual/tutorial/11-moving-bodies) | [1](/manual/cookbook#query-a-rotated-body-without-rebuilding-it) |
| Swept geometry of a time step: exact swept volumes, time-averaged apertures | `swept_cells_2d`, `swept_cells_3d` | [Moving bodies](/guide/advanced/moving-bodies) | [11](/manual/tutorial/11-moving-bodies) | [1](/manual/cookbook#swept-volume-of-a-moving-body-on-a-grid), [2](/manual/cookbook#swept-volume-of-a-moving-body-on-a-3-d-grid) |

## GPU

| Capability | Entry point | Reference | Tutorial | Cookbook |
|---|---|---|---|---|
| Distance, rich queries, wall distance, box classification, posed queries on the device (OpenACC), bit-exact against the CPU path | `enter_device`, `*_device`, `exit_device` | [GPU queries](/guide/advanced/device-distance) | [16](/manual/tutorial/16-gpu) | [1](/manual/cookbook#run-a-batch-of-queries-on-the-gpu), [2](/manual/cookbook#locate-many-points-on-the-gpu) |

## `fossilizer` CLI

`src/app/fossilizer.f90` wraps the library for the command line: load one or more STL files, clip, merge, mirror,
rotate, translate, resize, sanitize and repair connectivity, print statistics and, with `--mass_properties
[--density ρ]`, the mass properties of a closed body, or with `--compare other.stl`, the certified Hausdorff
distance to another surface, or with `--distance other.stl`, the body–body distance, or with `--patches prefix`, the named patches (one STL per patch), or with `--features file.vtk`, the feature edges and corners, or with `--thickness
file.vtk [--thickness_method ray|medial] [--h h --thin_factor f]`, the thinnest wall, the narrowest gap and the facets
thinner than the cell; save the result.

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
