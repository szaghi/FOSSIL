# Changelog

All notable changes to this project are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

## [2.0.0] — 2026-11-04
### Added
- **tests**: Add concurrency test for fossil distance calculations

- **ib_generator**: Print CLI options and rework xall export bounds

- **distance**: GPU-offloadable distance_many_device path (issue #20 Step 5)

- **geometry**: Add a plain-real, device-callable solid-angle core

- **winding**: Persistent per-node fast-winding-number cache

- **winding**: Iterative second-order fast-winding-number traversal

- **distance**: Add the SIGN_FAST_WINDING_NUMBER sign algorithm

- **ib_generator**: List fast_winding_number among the sign algorithms

- **repair**: Extract boundary loops for hole filling

- **repair**: Triangulate boundary loops by Liepa minimum-weight DP

- **repair**: Add surface%fill_holes

- **repair**: Optional hole filling in sanitize

- **fossilizer**: Add --fill_holes, --max_hole_edges, --max_hole_size

- **distance**: Fast-winding-number sign in distance_many_device

- **distance-field**: Narrow-band exact distances on rectilinear grids

- **distance-field**: Signs outside the band by flood fill

- **distance-field**: Far field by fast sweeping on stretched grids

- **distance-field**: Narrow band on the GPU (use_device)

- **ib_generator**: Far field from the distance field, not copied values

- **ib_generator**: --squared writes the former output quantity

- **arrangement**: Exact pair classification and canonical cut points

- **boolean**: Cut and classify coplanar overlaps

- **arrangement**: Close T-junctions; coplanar boolean matrix passes

- **csr**: Conjugate gradient solver for sparse SPD systems

- **repair**: Refine hole patches to the surrounding density

- **repair**: Fair hole patches with a bi-Laplacian solve

- **repair**: Refine and fair holes from sanitize and fossilizer

- **predicates**: Exact orient2d

- **curve**: Closed planar curves for 2-D geometry

- **curve**: Segment BVH, signed distance and self-intersections

- **curve**: 2-D winding number, inside test and open loops

- **curve**: Slice a surface by an axis-aligned plane

- **curve**: Repair loops (join, close, merge, orient)

- **curve**: Batched distances on the CPU and the device

- **distance-field**: 2-D fields from curves

- **surface**: Shells and bodies for surfaces and curves

- **surface**: Batched rich queries

- **surface**: Rich queries on the device

- **curve**: Batched rich queries on curves

- **distance-field**: Vector level set

- **query**: Probe and image points with visibility

- **query**: Report multiple foot points

- **box**: Box predicates for AMR tagging

- **box**: Classify AMR boxes against surfaces and curves

- **box**: Box classification on the device

- **box**: Signed tag bands on a block of cells

- **line**: Exact crossings of grid lines with the wall

- **line**: Grid-line crossings on the device

- **cut**: Exact 2-D cut cells

- **cut**: 3-D cut-cell faces from exact plane sections

- **cut**: 3-D cut-cell volumes by the divergence theorem

- **pose**: Rigid poses and body-frame point queries

- **pose**: Body-frame queries on the device

- **pose**: Rigid move with tree refit

- **pose**: Fresh and dead points between two poses

- **field**: Eikonal diagnostic of a distance field

- **query**: Wall distance with body id and early cutoff

- **query**: Curvature at the foot point

- **query**: Grid-adapted surface and curve markers

- **box**: Fast winding number as the default sign of box tagging ⚠ BREAKING CHANGE

- **curvature**: Pointwise-convergent curvature by a local jet fit

- **pose**: Swept geometry of a translating curve on a fixed grid

- **pose**: Exact swept geometry of a translating surface

- **predicates**: Exact orient3d

- **cut**: Fluid pieces of the 3-D cut cells

- **cut**: Measures of the fluid pieces of the 3-D cut cells

- **curve**: Optional self-intersection checks of curves and swept inputs


### Changed
- **aabb**: Inline aabb_object into aabb_node_object (issue #20 Step 2)

- **distance**: Iterative explicit-stack bvh traversal (issue #20 Step 3)

- **distance**: Device-callable bvh kernel chain (issue #20 Step 4)

- **distance**: Drop obsolete nvfortran 26.1 device workarounds (issue #20)


### Documentation
- **curvature**: Make curvature example snippets compile standalone

- **nvidia**: Correct manual-validation workflow + narrow tests-nvidia target

- **install**: Correct the minimum gfortran version (issue #24)

- **install**: Replace stale submodule and CMake instructions

- **advanced**: Add the GPU distance queries page

- **features**: List the batch and GPU distance queries

- **distance**: Correct the exit_device docstring

- **advanced**: Describe the BVH depth cap instead of the overflow risk

- **api**: State what the transforms leave stale and how to refresh it

- **api**: Document the re-analyze after sanitize_normals flips

- **api**: Reflections now preserve orientation

- Describe the generic sign rays, per-shell orientation and rotate

- Document the fast-winding-number sign and the rebuilt winding number

- Add the hole-filling page and list fill_holes

- Document the fast-winding-number sign on the device

- Document the signed-distance field

- Coplanar booleans

- Hole refinement and fairing

- **boolean**: Update stale module header comments

- 2-D geometry (curves)

- Immersed-boundary queries

- AMR tagging

- Sign algorithms at sharp features

- Cut-cell geometry

- Moving bodies

- Wall distance and field diagnostic

- Foot curvature and grid-adapted markers

- Swept geometry of moving bodies

- Exact swept geometry of a translating surface

- Exact orient3d and 3-D fluid pieces

- Measures of the 3-D cut-cell pieces

- **surface**: Correct the is_manifold and clip-box comments

- Refactor the documentation on the sibling skeleton

- A self-intersecting wall breaks the swept-geometry balance

- The regression and validation test tiers

- Self-intersection pairs sharing one vertex


### Fixed
- **distance**: Work around nvfortran 26.1 codegen bugs in device build (issue #20 Step 4)

- **docs**: Regenerate package-lock.json for vitepress 1.6 downgrade

- **docs**: Pin esbuild override to 0.25.x for vite 5.4 compatibility

- **ci**: Run install smoke test from the release workflow

- **build**: Exclude third-party docs/ and scripts/ sources from builds

- **build**: Exclude the FoXy benchmarks from the builds

- **build**: Exclude the StringiFor scripts from the builds

- **tests**: Avoid gfortran ICE on associate with an operator selector

- **aabb**: Cap the SAH BVH depth to bound the traversal stack

- **sanitize**: Rebuild derived data after sanitize_normals flips facets

- **facet**: Keep the outward orientation through reflections

- **facet**: Keep the orientation through an improper rotation matrix

- **distance**: Cast the ray-intersection sign rays off the axes

- **sanitize**: Orient each closed shell by its nesting parity

- **winding**: Use the persistent cache; correct sign and octree counts

- **distance-field**: Build the module into libfossil and re-export it

- **ib_generator**: Read list options with get_varying under nvfortran

- **repair**: Keep sharp edges when fairing hole patches

- **build**: Split code lines over 132 columns

- **distance**: Accurate point-triangle distance near the surface

- **remesh**: Keep sharp edges split during remeshing

- **build**: Split a test line over 132 columns

- **build**: Follow the VecFor _dev rename of device procedures

- **distance**: Robust closest point on needle triangles

- **device**: Restore the bit-exact nvfortran device build

- **aabb**: Drop the unused octree nodes when the depth is lowered

- **self-intersection**: Report crossings of facets sharing one vertex


### Performance
- **distance-field**: OpenMP fast sweeping, bit-exact with serial

- **distance-field**: Sign point queries on the device too


### Security
- **docs**: Untrack docs/package-lock.json and gitignore it


## [1.3.3] — 2026-05-14
### Performance
- **distance**: Split point-to-triangle kernel into lean and rich faces


## [1.3.2] — 2026-05-14
### Added
- **surface**: Add OpenMP-parallel distance_many batch query


## [1.3.1] — 2026-05-14
### Performance
- **distance**: Flatten BVH leaf facets into a contiguous distance SoA


## [1.3.0] — 2026-05-14
### Documentation
- **advanced**: Document cotangent Laplacian, curvature, smoothing + sync hero


## [1.2.13] — 2026-05-14
### Added
- **laplacian,surface**: Cotangent Laplacian + barycentric mass (issue #18 §2.1)

- **curvature,surface**: Per-vertex Gaussian curvature (issue #18 §2.4)

- **curvature,surface**: Per-vertex signed mean curvature (issue #18 §2.4)

- **smoothing,surface**: Explicit Laplacian + Taubin lambda|mu (issue #18 §2.3)


### Documentation
- **advanced**: Add user-facing pages for the 9 advanced primitives


## [1.2.12] — 2026-05-13
### Added
- **alpha-wrap**: Step 1 — leaf-size-α octree (issue #18 §1.6)

- **alpha-wrap**: Step 2 — inside/outside flood fill (issue #18 §1.6)

- **alpha-wrap**: Step 3 — boundary surface extraction (issue #18 §1.6)

- **alpha-wrap**: Step 4 — vertex projection toward input (issue #18 §1.6)

- **alpha-wrap**: Step 5 — surface%alpha_wrap capstone TBP (issue #18 §1.6)


## [1.2.11] — 2026-05-13
### Added
- **fossil**: Implement SDF functionality with AABB tree ray query and test coverage

- **sdf**: Implement Laplacian smoothing for SDF values with test coverage

- **sdf,surface**: Step 3 — GMM clustering + segment_sdf TBP (issue #18 §1.9)


## [1.2.10] — 2026-05-13
### Added
- **fossil**: Implement ray query functionality with test coverage

- **fossil**: Implement AABB tree ray query with test coverage

- **ray-query,surface**: Step 3 — any-hit early-exit (issue #18 §2.5)


### Fixed
- **fossil_decimate**: Prevent out-of-bounds heap access with nested guards


## [1.2.9] — 2026-05-13
### Added
- **remesh**: Isotropic remeshing foundation, steps 1-4 (issue #18 §1.7)

- **remesh**: Step 5 — edge flip pass for valence balancing (issue #18 §1.7)

- **remesh**: Step 6 — relax + projection (issue #18 §1.7)

- **remesh,surface**: Step 7 — feature lock + surface TBP (issue #18 §1.7)


## [1.2.8] — 2026-05-13
### Added
- **decimate,surface**: QEM mesh decimation (issue #18 §1.3)


## [1.2.7] — 2026-05-13
### Added
- **mc,surface**: Marching cubes isosurface extraction (issue #18 §1.5)


## [1.2.6] — 2026-05-13
### Added
- **predicates,arrangement**: Boolean foundation (issue #18 §1.1 stage 1)

- **dt**: 2D Delaunay triangulation, Phase 1 (issue #18 §1.1 stage 2)

- **dt**: Constrained Delaunay triangulation, Phase 2 (issue #18 §1.1 stage 2)

- **boolean**: Mesh-mesh booleans, BOOL_DIFFERENCE MVP (issue #18 §1.1)

- **facet**: Clip tri-tri intersection segments to both triangles (issue #18 §1.1)

- **boolean**: Coplanar-aware shared-boundary tagging (issue #18 §1.1)

- **boolean**: Wire UNION/INTERSECT/SYMDIFF + complete §1.1 (issue #18 §1.1)


## [1.2.5] — 2026-05-13
### Added
- **surface**: Add self-intersection detection (issue #18 §1.2)


## [1.2.4] — 2026-05-12
### Added
- **surface**: Add generalized / fast winding number (issue #18 §1.4)


### Documentation
- **guide**: Add API companion guide with CI-verified snippets


## [1.2.3] — 2026-05-12
### Added
- **facet,surface**: Expose facet%area() and surface%get_area() (#7)


## [1.2.2] — 2026-05-12
### Added
- **rotate**: Add optional `center` argument for body-frame rotation (#6)


## [1.2.1] — 2026-05-12
### Added
- **lib**: Implement vertex pool object for fossil surface handling


### Changed
- **surface**: Make vertex pool the source of truth for coincidence (#5)


### Fixed
- **clip**: Use full SAT for triangle-AABB overlap (#1)


## [1.2.0] — 2026-05-12
### Added
- **distance**: Correct AABB traversal and add pseudo-normal sign

- **surface**: Implement connectivity symmetry checks for fossil surfaces

- **surface**: Expand sanitization pipeline and add validity predicates

- **aabb**: Auto-tune octree refinement depth from facet count

- **aabb**: Scaffold tree-kind selector for upcoming SAH BVH

- **aabb**: Implement SAH BVH builder

- **aabb**: Make SAH BVH queryable via tree-kind dispatch


### Documentation
- Refresh user-facing docs after the signed-distance arc


### Fixed
- **block**: Avoid associate alias as allocate bound; wrap long lines

- **surface**: Orient sanitize_normals outward, populate pseudo-normals ⚠ BREAKING CHANGE


### Performance
- **aabb**: Switch default tree to SAH BVH


## [1.1.0] — 2026-05-11
### Added
- **api**: Add optional status argument to mutating procedures (audit #14 S1)

- **api**: Add optional status argument to mutating procedures (audit #14 S1)


### Changed
- **aabb-tree**: Encapsulate aabb_tree_object state via accessors

- **surface,file**: Encapsulate components and remove stale API docs ⚠ BREAKING CHANGE

- **facet**: Collapse edge fields and replace edge_dir with integer ⚠ BREAKING CHANGE

- **surface**: Remove unimplemented compute_mesh_distance ⚠ BREAKING CHANGE

- **api**: Replace sign_algorithm string with integer enum ⚠ BREAKING CHANGE

- **surface**: Encapsulate surface fully, delete file_stl_object ⚠ BREAKING CHANGE


### Fixed
- **api**: Define intent(out) and stop on invalid enums ⚠ BREAKING CHANGE

- **api**: Guard centroid, strict ascii prefix, add finalisers

- **io**: Size-based STL probe and EOF/iostat guards on facet reads

- **surface**: Replace order-dependent vertex merge with union-find

- Drop custom assignment overloads, rely on intrinsic assignment

- **surface**: Error stop on ambiguous translate/resize arguments

- **api,surface**: Rename analize, warn on residual disconnected facets after sanitize ⚠ BREAKING CHANGE


## [1.0.12] — 2026-05-09
### Documentation
- Fix src_dir path and update package dependencies

- **cliff**: Remove unused github remote section from cliff.toml


### Fixed
- **fobos**: Correct gcov_analyzer flag syntax in makecoverage-analysis rule

- **scaffold**: Update release script to handle fobos directory changes

- **scaffold**: Improve trunk branch detection in release script

- **scaffold**: Handle missing trunk branch detection in release script

- **scaffold**: Improve git fetch error handling in release script


### ci
- Modernize CI/CD workflows and actions ⚠ BREAKING CHANGE


## [1.0.7] — 2026-02-21
### Added
- Add full VitePress guide and modernise project tooling


## [1.0.6] — 2024-09-24
### Fixed
- **update**: Update third party submodules



