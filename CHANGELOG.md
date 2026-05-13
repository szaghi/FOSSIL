# Changelog

All notable changes to this project are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

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



