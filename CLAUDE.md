# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is FOSSIL?

FOSSIL (FOrtran Surfaces, Solids, Immersed boundaries and Level sets; until 2026-11 "FOrtran Stereo Litography parser") is a pure Fortran 2003+ geometry library for simulation codes: STL surfaces and 2-D curves (I/O, repair, processing), signed distance and level sets, immersed-boundary queries, AMR tagging, exact cut cells and moving bodies, on the CPU or bit-exact on the GPU. It is OOP-designed and TDD-designed. The companion CLI app `fossilizer` (`src/app/fossilizer.f90`) wraps the library for command-line STL processing.

## Build System

The project uses **FoBiS.py** (Fortran Build System) with a `fobos` configuration file. Install once via pip (`pip install FoBiS.py` — the PyPI package name is unchanged), then invoke as `fobis` (the modern CLI binary).

### Common build commands

```bash
# Build static library (GNU compiler)
fobis build --mode static-gnu

# Build shared library (GNU compiler)
fobis build --mode shared-gnu

# Build all tests (GNU compiler)
fobis build --mode tests-gnu

# Build with debug flags (GNU compiler)
fobis build --mode tests-gnu-debug

# Clean build artifacts
fobis rule --ex clean
```

Intel compiler variants are available by replacing `-gnu` with `-intel`.

### NVIDIA / device build (OpenACC, issue #20)

The `static-nvidia` and `tests-nvidia` modes target the GPU via OpenACC with
`nvfortran`. They are **manual-only** — not built by CI, but **validated
locally** by loading the nvhpc module before the build:

```bash
module load nvhpc/26.9
fobis build --mode tests-nvidia --varset local_nvf
```

**Requires nvfortran >= 26.9.** Validated on 26.9 (RTX 4070, `cc89`; last on 2026-11-05 at 8c78a46 with
bunny, cube-01, dragon-fine at 32³, every device-vs-serial line 0); 26.1 is no
longer supported, since the query-point mapping workaround it needed was removed.
In a non-interactive shell `module load` may not populate `PATH`; export the SDK's
`compilers/bin` and `compilers/lib` (+ `cuda/lib64`) directly instead.

Device-code rules (issue #20; see the comments in the named files):

- No type-bound procedure calls or `class(...)` dummies in device code: use the
  free procedures: FOSSIL's `_oac` ones and VecFor's `_dev` ones (VecFor >= 1.9.0
  renamed its `_oac` suffix to `_dev`, no aliases). This includes VecFor's overloaded
  `=`, so set `vector_R8P` components directly (or call `assign_vector_R8P_dev`).
- A `!$acc routine seq` procedure that uses `vector_R8P` must not be added to
  `fossil_facet_object`: it triggers an nvfortran NVVM duplicate-type bug (still
  present in 26.9). Put it in its own module, like `fossil_triangle_closest`, and keep ONE such
  routine per module: two in the same module trip it too (helpers take plain reals).
- `enter_device` maps member arrays, never `self`: inside a compute region,
  reference them through top-level pointer aliases, never as
  `self%aabb%payload(i)`. The mapped arrays are `aabb%node`, `aabb%payload`,
  `pnormal_payload` and `wn_payload` (aliases `node_p`, `payload_p`,
  `pnormal_p`, `wn_p` in `distance_many_device`).
- The device sign algorithms are `SIGN_PSEUDO_NORMAL` and
  `SIGN_FAST_WINDING_NUMBER` (issue #30). The FWN branch calls
  `bvh_winding_number_oac` with `FOUR_PI` / `DEFAULT_BETA` exported by
  `fossil_winding_number`, so the CPU and device paths share routine and
  constants: keep it that way, or bit-exactness breaks.

The `local_nvf` varset supplies `$NVF_CC` (GPU compute capability, default
`cc89`); override for a different GPU (`cc80` for A100/V100). The default
varset is `local_gnu`, a no-op placeholder, so every GNU/Intel mode is
unaffected.

The `tests-nvidia` mode builds a **single targeted test** (`src/tests/validation/fossil_test_distance_bench.f90`)
that exercises the distance kernel chain — not the full suite. Host-only tests
(file I/O, mesh-quality invariants) stay in `tests-gnu` / `tests-intel`; some
tripped nvfortran 26.1 parser quirks unrelated to the device port (e.g. parser
rejects `(expr) .cross. (expr)` in one strict-normal test).

The `tests-nvidia-curve` mode is its 2-D twin (issue #36): it builds
`fossil_test_curve_bench` only (segment BVH closest point, 2-D fast winding number on
the device). Same clean-build rule (shared `exe/`), pass criterion
`device vs serial max |d| err 0.000E+00`; the `rich`, `box`, `posed`, `sign bench`, `line crossings`, `deviation`, `feature` and `thin` mismatch lines must be 0 too (the 3-D bench prints the same lines). Validated on the GPU 2026-11-05 (26.9, RTX 4070, 512²).

The NVIDIA modes compile with `-gpu=nofma -Mnofma` (`$CSTATIC_NV`): with FMA contraction the
device lost host bit-exactness (1-ulp distances, sign and facet-id flips at near-ties); the device
kernels are 15-30 % slower for it. Do not remove it unless the gates stay at 0.

The device build scope is **BVH-only** (the octree stays CPU-only).

### Running tests

```bash
# Build tests first, then run all
fobis build --mode tests-gnu
./scripts/run_tests.sh

# Run a single test executable
./exe/fossil_test_clip

# Full coverage analysis (build debug + run tests + gcov)
fobis rule --ex makecoverage
```

Test programs output `Are all tests passed? T` on success or `Are all tests passed? F` on failure.

Two tiers (issue #45). **Regression** (`src/tests/*.f90`, every push, CI builds `tests-gnu-debug --coverage`): each
test <= 2 s in that debug coverage build (`boolean_coplanar` 3 s is the one exception), enforced by
`scripts/check_test_budgets.sh` in the `makecoverage` rule (exceptions in `src/tests/budgets.txt`; factor 2 on CI).
Re-time a test in that build after the last change of a step, not mid-step. Choose inputs for it:
`dragon.stl` for self-consistency checks, a clean closed body (`cube-0N.stl`, `make_apex_body`, a marching-cubes
sphere) for exact-sign checks, not `bunny.stl` (1.4 s to load in debug since issue #57, 19 s before; its tests
are in validation); never loosen a tolerance to save time.
**Validation** (`src/tests/validation/*.f90`: benches, large meshes, sweeps, the full-size `fossil_val_*` copies,
`fossil_test_distance_bench`), nightly (`.github/workflows/validation.yml`) and before a release:

```bash
fobis build --mode validation-gnu    # -O2, into exe-validation/
./scripts/run_validation.sh          # or: ./scripts/run_validation.sh fossil_val_sdf
```

Every `tests-*` mode excludes `src/tests/validation/`. `fobis clean` leaves stale executables of
moved sources in `exe/`: delete them by hand, or `run_tests.sh` runs them.

### Documentation

```bash
fobis rule --ex makedoc   # Requires 'formal' and npm/VitePress
fobis rule --ex deldoc    # Delete doc build artifacts
bash scripts/docs_examples.sh   # rebuild static-gnu, compile and run docs/examples/src/*.f90, regenerate
                                # docs/examples/{snippets,output,images} (~5 min); never edit those by hand
cd docs && npm run docs:build   # API pages (formal) + VitePress site
```

The site follows the FLAP/VecFor/StringiFor/FoXy skeleton (issue #43): Start here (`docs/guide/index.md`,
`installation.md`) → Tutorial (`docs/manual/tutorial/NN-*.md`, 16 chapters, the GPU last) → Cookbook (`docs/manual/cookbook.md`) →
Reference (`docs/guide/features.md` is the feature map, then `usage.md`, `advanced/*.md`, `api/*.md`) → Project. One
sidebar in reading order in `docs/.vitepress/config.mts`. Every tutorial/cookbook snippet and output comes from a
program in `docs/examples/src/` (markers `!run`, `!region`, `!cast`, `!anim`; see the script header); explicit output
formats only. `scripts/anim2svg.py` draws the SMIL geometry animations (`!anim`), `scripts/ansi2svg.py` the terminal
casts, `scripts/mesh2svg.py` the figures (`!fig NAME`: the program writes NAME.fig through
`docs/examples/lib/fossil_fig.f90`; issue #47). The landing gallery and the README's are generated by
`scripts/gallery.py` (regions `<!-- gallery:start/end -->`; never edit them by hand; `--check` runs in
`check_doc_snippets.sh`). Brand: `docs/brand/marks.py` generates `docs/public/{logo,logo-dark,favicon,favicon-dark}.svg` (the Sumi ammonite; PNG favicons rasterised from `favicon.svg`; never hand-edit the SVGs);
palette in `docs/.vitepress/theme/custom.css`. `scripts/check_doc_snippets.sh` still compiles the fortran blocks of
the hand-written reference pages, and runs `scripts/check_doc_coverage.py --strict` (issue #77): every reference page
and tutorial chapter shows a generated figure and a run example (`ref_<page>.f90`, `fig_ref_<page>_*.f90`), and
(issue #85) a figure in every explaining `##` section (exempt titles in the script), an animation on the pages of
iterative or time-dependent features (`ANIMATED`) and in every tutorial chapter (with at least two runs), every
feature-map entry point in a cookbook recipe's snippet; every feature-map row links a tutorial chapter or a recipe
(Tutorial and Cookbook columns), the landing cards equal the gallery tiles, every example program is shown somewhere,
every included file exists. A new feature needs all of them.

### Alternative build (fpm)

```bash
fpm build --profile release   # dependencies resolved from fpm.toml
```

There are no CMake/Make build files: `scripts/install.sh --build cmake|make` skips.

## Source Architecture

```
src/
├── lib/           # Core library modules
├── tests/         # Test programs (one per feature)
├── app/           # fossilizer CLI application
└── third_party/   # fobos [dependencies], fetched by `fobis fetch`
```

### Core object hierarchy

The public API is accessed via `use fossil`, which re-exports four types plus the named constants:

| Type | Module file | Purpose |
|------|-------------|---------|
| `surface_stl_object` | `fossil_surface_stl.f90` | Surface geometry: I/O, analysis, sanitization, manipulation, distance queries |
| `facet_object` | `fossil_facet_object.f90` | Individual triangular facet — vertices, normal, connectivity, pseudo-normals |
| `distance_field_object` | `fossil_distance_field.f90` | Signed-distance field on a rectilinear (stretched) grid: narrow band, flood-fill signs, fast sweeping; 2-D fields from curves |
| `curve_object` | `fossil_curve_object.f90` | 2-D geometry (issue #36): closed/open polylines in an axis-aligned plane, segment BVH, signed distance, 2-D fast winding number, repair, batched CPU/device queries |

`fossil.f90` also re-exports:

- `SIGN_PSEUDO_NORMAL` (default), `SIGN_RAY_INTERSECTIONS`, `SIGN_SOLID_ANGLE`, `SIGN_FAST_WINDING_NUMBER` — sign algorithms for signed-distance queries.
- `AABB_TREE_SAH_BVH` (default), `AABB_TREE_OCTREE` — AABB tree kinds.
- `AABB_AUTO_REFINEMENT` — sentinel for octree depth auto-tuning.
- `BOX_CUT`, `BOX_NEAR`, `BOX_INSIDE`, `BOX_OUTSIDE` — box classification codes (issue #38).
- `STATUS_OK`, `STATUS_ALLOC_FAIL`, `STATUS_AMBIGUOUS_ARGS`, `STATUS_FILE_NOT_FOUND`, `STATUS_FILE_OPEN_FAIL`, `STATUS_INVALID_INPUT`.

`file_stl_object` was deleted in commit `7ffe7f8` — all I/O is now type-bound on `surface_stl_object` (`load_from_file`, `save_into_file`).

### Typical usage flow

```fortran
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)  ! runs analyze internally
call surface%sanitize                                                    ! full repair pipeline
d = surface%distance(point=p, is_signed=.true., is_square_root=.true.)   ! SAH BVH + pseudo-normal
call surface%save_into_file(file_name='out.stl')
```

1. `surface%load_from_file(file_name=, guess_format=, ..., aabb_tree_kind=, status=)` — reads STL, runs `analyze` (bounding box, volume, centroid, connectivity, AABB tree, pseudo-normals). Returns `STATUS_INVALID_INPUT` for NaN/Inf coords.
2. `surface%sanitize` — orchestrator: `remove_degenerate_facets` → (`connect_nearby_vertices` if needed) → `analyze` → `remove_duplicate_facets` → `analyze` → `sanitize_normals` → (optional `fill_holes`, with `fill_holes=.true.`; `refine_holes` / `fair_holes` refine and thin-plate-fair the patches, fairing solved by `fossil_linear_solver%cg_solve`). Emits one stderr warning per defect class.
3. `surface%analyze(aabb_refinement_levels=, aabb_tree_kind=)` — invoked automatically by load/sanitize/clip. Rebuilds connectivity and the AABB tree. The optional args are passed through to `aabb%initialize`.
4. Manipulation: `clip`, `rotate`, `translate`, `mirror`, `resize`, `merge_solids`.
5. Queries: `distance`, `is_point_inside`, plus `compute_distance` (richer form with closest facet/edge/vertex, closest point and normal outputs). IB queries (issue #37): `query_point` / `query_many` / `query_many_device` (closest point, facet, `FEATURE_FACE/EDGE/VERTEX`, pseudo-normal, body; optional `ambiguous` / `second_foot_points` on the CPU), `probe_many`, `image_points`; shells and bodies from `analyze` (`facet_shell`, `facet_body`; `body_payload` is device-mapped). Curves have the same set. `distance_field_object%compute(store_closest=.true.)` keeps `closest`, `closest_id`, `closest_normal` on band cells. AMR tagging (issue #38): `box_classify_many` / `box_classify_many_device` (`BOX_CUT/NEAR/INSIDE/OUTSIDE`), `box_intersects_surface`, `tag_block` (signed Chebyshev tag bands on a block); curve twins. These default to `SIGN_FAST_WINDING_NUMBER` (issue #39); the point queries keep `SIGN_PSEUDO_NORMAL`.
6. Validity predicates: `is_watertight()`, `is_manifold()`, `is_volume()`.
7. `surface%save_into_file(file_name=, is_ascii=, status=)` — writes result.

The components of `surface_stl_object` are **private**. All reads go through `get_facets_number()`, `get_bmin()`, `get_bmax()`, `get_volume()`, `get_centroid()`, `get_non_manifold_edges_number()`, `get_degenerate_facets_removed()`, `get_duplicate_facets_removed()`, and `facet_at(i)` (pointer-returning, may be `null()` for out-of-range).

### AABB acceleration

`fossil_aabb_tree_object.f90` provides **two tree kinds** behind the same traversal code:

- **`AABB_TREE_SAH_BVH` (default)** — binary BVH built top-down by partitioning *triangles* with the bucketed surface-area heuristic (16 buckets per axis). Adapts to triangle density. On dragon-fine (24k facets, 32³ query grid): ~110× faster than the octree.
- **`AABB_TREE_OCTREE`** — legacy 8-way space-partitioning octree, kept for benchmarking and as a fallback. Depth is set by `refinement_levels` (or auto-tuned via `AABB_AUTO_REFINEMENT`).

Selection is per-surface via the optional `aabb_tree_kind` argument on `load_from_file`, `adopt_facets`, or `analyze`. Both kinds produce **bit-exact** distance results — best-first traversal with d² pruning never prunes the true closest facet.

### Signed-distance sign algorithm

Set via the optional `sign_algorithm` argument on `distance` / `compute_distance` / `is_point_inside`:

- **`SIGN_PSEUDO_NORMAL` (default)** — Bærentzen–Aanæs angle-weighted pseudo-normal at the closest point. Fused with the distance traversal (no second pass). Requires outward-oriented normals (produced by `sanitize`).
- `SIGN_RAY_INTERSECTIONS` — ray casts in three fixed generic directions (majority vote); odd intersections = inside.
- `SIGN_SOLID_ANGLE` — sum of projected solid angles ≈ ±4π (exact, O(N) per query).
- `SIGN_FAST_WINDING_NUMBER` — fast winding number > 0.5 over the per-node cache `wn_payload` built by `analyze`; robust on surfaces with holes.

Measured by the issue #39 benchmark (`fossil_test_sign_bench`, `_2d`, `_device`; bodies in `src/tests/fossil_sign_bench_bodies.f90`, a test-side module): all four exact on watertight meshes down to 1e-12 L, needles of aspect 1e3 included; FWN the most robust on dirty meshes. The point-triangle squared distance is the norm of the residual (`triangle_closest_st`), not the Eberly quadratic form: the latter is noise below ~1e-8 L and flipped pseudo-normal signs. Needle triangles (det < 1e-4 max(a, c)^2) take `robust_st` (triple products with n, error ~eps/|n| instead of Eberly's ~eps/|n|^2); the others keep the Eberly solve (20 % cheaper). Keep both (`fossil_test_distance_tiny` guards them, the plate included).

### Third-party dependencies (`src/third_party/`)

| Dependency | Key exports used in FOSSIL |
|-----------|---------------------------|
| **PENF** | Portable numeric types: `I4P` (int32), `R8P` (real64), `str()` |
| **VecFor** | 3D vectors: `vector_R8P`, unit vectors `ex_R8P/ey_R8P/ez_R8P`, `rotation_matrix_R8P`, `mirror_matrix_R8P` |
| **FLAP** | CLI argument parsing (`command_line_interface`) used by `fossilizer` |
| **FACE** | ANSI color/style output |
| **BeFoR64** | Base64 encoding |
| **FoXy** | XML/XDMF support |
| **StringiFor** | String utilities |
| **VTKFortran** | VTK file I/O |

They are not git submodules: they are declared in `fobos` `[dependencies]` and fetched
(pinned by `src/third_party/fobos.lock`) with:
```bash
fobis fetch            # --update to re-fetch
```

### Key internal modules

- `fossil_utils.f90` — constants (`EPS`, `PI`, `FRLEN`) and helpers (`is_inside_bb`)
- `fossil_list_id_object.f90` — dynamic integer ID list used for connectivity tracking
- `fossil_aabb_object.f90` / `fossil_aabb_node_object.f90` — AABB primitives. Ray queries (issue #82): one hit per
  crossing by exact tests (the ray as a segment to a far point past the box, rounded once; three orient3d side tests
  with one symbolic perturbation of the endpoints; axis-aligned rays through the `fossil_line_kernel` kernel shared with
  `line_crossings`, `tri_line_hit_oac` included); octree node boxes padded by the tree's `reach` (facets are assigned by
  centroid and stick out of their node). `intersect_ray_all(exclude=)`; t = 0 hits kept. `surface%tree_boxes` /
  `curve%tree_boxes` (issue #89): node boxes breadth first down to `max_depth`, both tree kinds; `save_aabb_into_file`
  uses it.
- `fossil_block_object.f90` / `fossil_block_aabb_object.f90` — structured mesh block with AABB (`extrapolate_distances` deprecated, issue #32)
- `fossil_distance_field.f90` — `distance_field_object` (re-exported by `fossil`, issue #32): signed-distance field on a rectilinear grid; narrow band via `distance_many(_device)`, flood-fill signs, fast sweeping in hyperplane order (OpenMP bit-exact). It imports the surface/facet modules, never `fossil` (circular). `fossil_ib_generator` is built on it. `initialize` without `z` makes a 2-D field for a curve source; `compute_band` / `compute_signs` / `compute` / `default_band` are generics over `surface=` / `curve=`.
- `fossil_curve_object.f90`, `fossil_segment_bvh.f90`, `fossil_curve_slice.f90` — 2-D geometry (issue #36): `curve_object` stores (u, v) as SoA arrays in the plane normal to `normal_axis` (u, v = the other axes in cyclic order); the solid is on the left of every segment (outer loops CCW, holes CW), signed distance negative inside; open loops allowed (`# open` in the ASCII format). The segment BVH's traversals are `_oac` kernels shared by the CPU and the device. `slice_surface` cuts an STL by an axis-aligned plane (canonical crossings, no tolerances). `fossil_predicates%orient2d` and `orient3d` (`orient3d_xyz`) are exact (filter + exact expansion with a bit-mask split, safe under nvfortran FMA contraction; orient3d, host only: exact zero on a shared coordinate plane, 6 triple products of exact differences or 24 of input coordinates, each exact as four doubles by error-free products (Dekker, bit-mask split; no `ieee_fma`, which nvfortran 26.9 lacks) since issue #65; no `predicates_initialize` needed).
- Issue #42 (WP6): `distance_field_object%gradient_diagnostic` (eikonal q per node), `wall_distance_many(_device)` (BVH started at cutoff^2; curve twin). `query_many(curvatures=)` (foot curvature: (H, K) from `jet_curvature` — cubic fit over the 2-ring, O(h^2) pointwise; the #18 cotangent / angle-defect operators stall near 15 % at valence-5 vertices — by barycentric weights; 2-D turning angle; CPU only), `surface_markers` (remesh to sqrt(3) dx + one calibration pass) / `curve_markers` (exact arc-length resampling). `isotropic_remesh(preserve_features)` re-detects features after each split (else crease midpoints were relaxed off and edges chamfered). Since issue #79 the collapse pass coarsens (link condition, one collapse per vertex per pass, no edge above 4L/3, locked endpoint kept); since issue #91 the relaxation is tangential and `REMESH_ERROR_TOL` = 0.07 L refuses any collapse, flip or move whose new facets leave the input farther (cube-01 at L = 1 without features keeps 0.949 of its volume). `decimate` (issue #81): link condition and zero-area refusal; `target_facets <= 0` is BAD_INPUT; patch rims coarsen along themselves onto input rim vertices (`rim_collapses`, default on, corners fixed; issue #91); optional `max_error` ceiling (a heuristic, not a certificate). `extract_isosurface` (issue #87) snaps a node to the level within 1e-12 of its edge crossings (round-off slivers). `jet_curvature` merges near-coincident ring vertices and returns NaN + `CURV_STATUS_FIT_FAILED` when the fit stays rank-deficient (issue #91). `fill_holes(manifold_only=)` / `sanitize(manifold_holes=)` chain loops through true boundary edges only (opt-in).
- `fossil_rigid_pose.f90` — moving bodies (issue #41): `rigid_pose_t` / `rigid_pose_2d_t`; optional `pose=` on the point queries (body-frame, CPU and device), `move(pose)` with `aabb_tree%refit` / `segment_bvh%refit_vertices` (O(N); `aabb_node%set_box` writes the box only — `node%initialize` wipes children), `state_changes`; `fossil_swept_cells.f90`: `swept_cells_2d(curve, ...)` / `swept_cells_3d(surface, ...)` (free routines, not TBPs: the module uses curve/surface/cut_cells_3d) — swept volume = cut-cell volume change (exact), face apertures time-averaged (translation exact via breakpoints: 2-D per face, 3-D per grid plane — mesh vertex on the moving plane, mesh edge on a moving grid line, grid node on a triangle; area quadratic per piece, 2-point Gauss; `exact=.false.` forces Gauss. Rotation: 2-D exact since issue #69 (`faces_rotation_2d`: certified breakpoints via `isolate_zeros` on g(s) = (p0+p1 s) cos phi + (q0+q1 s) sin phi + k0 + k1 s with the |g''| bound, Gauss doubled to round-off per piece, faces outside the swept circle box closed form); 3-D exact too (`faces_rotation_3d`, per grid plane: vertex on plane, edge on grid line, node on triangle, all of the same g form by Rodrigues; zero zones and 1e-10 merge in `isolate_zeros`); Gauss-Legendre `nq` with `exact=.false.` or explicit poses, nq-vs-2nq estimate), `wall` = -sum A_f n_f, `residual` = rigid-field volume balance; `pose_interpolate` (screw motion), `gauss_legendre_01`.
- `fossil_mass_properties.f90` — mass properties (issue #59, roadmap #58 WP1): pure free routines (compensated moment accumulator over signed tetrahedra, cyclic Jacobi `symmetric3_eigen`); `surface%mass_properties(mass, center, inertia, principal_moments, principal_axes, density, about, volume, status)` sums about the bounding-box centre (a body at 1e6 keeps its tensor), vertices only, rejects non-watertight or V <= 0 (`STATUS_INVALID_INPUT`). `compute_centroid` uses the same vertex-only tetrahedra (it scaled with the facet size before #59). `inertia_to_world` (R I R^T), `principal_frame` (rigid_pose_t); 2-D twin `curve%mass_properties` (Green sums per unit depth: `second_moments` = [int u^2, int v^2, int uv] about the centroid, `polar_moment`, `principal_angle` of the long axis).
- `fossil_deviation.f90` — surface deviation (issue #60, roadmap #58 WP2): one-sided Hausdorff h(A, B) through a batched distance callback (no surface type in the module); `surface%deviation(other, lower, mean, rms, per_vertex, ...)` samples vertices, edge midpoints and sub-triangle centroids (`samples_per_area`, default 4 per mean facet area; centroid quadrature for mean/rms), `surface%hausdorff`; `DEV_STATUS_*`. Certified `upper` by branch and bound (`deviation_certified`): per triangle min(max d_k + rho_T, min_k max_j |c_j - q_k|) (Lipschitz; convexity of |x - foot|, quadratic convergence), or max d_k when the three feet share a facet of B (coincident surfaces settle at once), split in four until upper - lower <= tolerance (default 1e-3 bbox diagonal) or `max_triangles` (DEV_STATUS_NOT_CONVERGED). B enters as an extension of `deviation_evaluator_t` (an internal procedure as actual argument made gfortran put trampolines on an executable stack). `device=.true.`: batches through `query_many_device` (other entered), bit-identical; bench line `deviation device vs serial mismatches`. 2-D twin `curve%deviation`/`hausdorff` (`deviation_sampled_2d`, `deviation_certified_2d`: segments halved; `deviation_evaluator_2d_t`, `curve_evaluator_t`), device line in `fossil_test_curve_bench` too. Adopted (issue #60 step 5) as certified checks in the remesh, decimate, smoothing, hole_fair, hole_refine and alpha_wrap tests (bounds measured, stated with margins). Coplanar surfaces with crossing facets refine along the crossings (cost ~1/tolerance): use the default tolerance there.
- `fossil_feature_edges.f90` — feature edges (issue #61, roadmap #58 WP3): `find_feature_edges` (stable merge sort of half-edges, sharp by `is_sharp_dihedral` with unit normals, boundary, non-manifold; corners; polylines CSR), `surface%feature_edges(edges, kinds, corners, polyline_start, polyline_vertices, angle, ...)`; constants `FEATURE_EDGE_SHARP/_BOUNDARY/_NONMANIFOLD` (not the query `FEATURE_EDGE`), default 30 degrees. `isotropic_remesh` locks with the same predicate and `feature_cos` (bit-identical). `build_features` (explicit: queries are intent(in); dropped by `analyze`/`move`) stores segments + BVH (`fossil_feature_bvh`, plain arrays, `feature_box_touch_oac` shared CPU/device); `box_feature_many(_device)`; `tag_block(feature_band=, feature_tags=)`; enter_device maps `feature_seg/box/link` if built; `query_many(_device)(feature_distances=)` (`feature_point_distance_oac`, exact clamped projection, bit-identical to brute force); bench line `feature device vs serial mismatches` (boxes and distances). 2-D: `curve%feature_corners` (turning angle, open-loop ends; computed per call, no build), curve `tag_block(feature_band, feature_tags, feature_angle)`, `query_many(_device)(feature_distances, feature_angle)`; curve bench line too. `save_features_into_file` (legacy ASCII VTK POLYDATA: POINTS, VERTICES = corners, LINES = edges, cell scalar `kind`), `fossilizer --features file.vtk [--feature_angle a]`.
- `fossil_thickness.f90` — local thickness and gaps (issue #62, roadmap #58 WP4): pure helpers (fixed golden-angle cone, lower median, `ball_radius`); `surface%thickness(thickness, gap, method, cone_angle, rays, vertex_thickness, vertex_gap, invalid, status, separation_angle)` per facet: THICKNESS_RAY (shape diameter, 60-degree cone of 31 rays from the centroid, first hit by `intersect_ray_first(exclude=)`, leaving/entering test, lower median; cone bias in closed form) or THICKNESS_MEDIAL (shrinking ball started at half the axial chord, so never above it; batched `query_many`; Ma's separation-angle stop, default 30 = local feature size, small next to sharp convex creases: the default is an open owner decision). `thin_regions` (union-find over edges), `build_thin` (local size per BVH payload triangle + per-node minimum; `analyze` drops, `move` keeps), `box_thin_many(_device)` through `bvh_box_thin_oac`, `tag_block(thin_factor, thin_band, thin_tags)`, enter_device maps `thin_pay/thin_node` if built; bench line `thin device vs serial mismatches` (3-D and curve). 2-D: `curve%thickness` (ray fan from the segment midpoint, `ray_first_2d`; medial disc), `thin_regions` (runs along loops), `build_thin`, `box_thin_many(_device)` (`curve_box_thin_oac`, exact), curve `tag_block(thin_*)`. `save_thickness_into_file` (legacy VTK), `fossilizer --thickness file.vtk [--thickness_method] [--h --thin_factor]`. `segment_sdf` keeps its own SDF (unchanged).
- `fossil_line_query.f90`, `fossil_cut_cells_2d.f90`, `fossil_cut_cells_3d.f90` — cut geometry (issue #40): exact axis-aligned line crossings (orient2d on the projected triangle, symbolic tie-break; exact facet vertices in `surface%tri_vertices`, device-mapped), 2-D cut cells (exact grid walk, Green, union-find pieces), 3-D faces from `slice_surface(on_plane_low=.true.)` sections and cell volumes by the divergence theorem. Rule everywhere: a point on a grid line/plane belongs to the low side (cells (x_i, x_i+1]). 3-D pieces (`cut%pieces`, -1 = untrusted cell): per-cell union-find over flagged clipped wall polygons + cell-edge fluid intervals (corner states = exact perturbed node parities, piercings = exact edge crossings under the global perturbation, issue #70), joined by shared vertices, piercings, fluid corners, face islands (ray); touching sets count, closed sets only if a void. Per-piece measures (`np`, `first_piece`, `piece_cell`, `piece_volume`, `piece_centroid`, `piece_wet_area/_centroid/_normal`, `piece_face(6, np)`): wall terms per set, face areas by Green over edge intervals + traces oriented by the polygon winding (high face as wound, low reversed), inclusions attached by a +x ray; single-piece cells take the totals. One global symbolic perturbation (issue #69): grid points at P + (eps1, eps2, eps3), the line kernels order ties x, y, z (`psign(bfirst)` on y lines), `line_crossing_before(_2d)` decides a crossing against a point exactly (one orient3d / orient2d, ties by the perturbed sign): node parities agree across axes. Since issue #70 traces are read on the perturbed face (in-plane polygons by their solid side, facets ending on the plane only if the facet across crosses it) and a multi-piece cell is trusted only if `split_closes`; lattice bodies: 19 of 1.57 million cells still flagged (947 before; a facet touching the cell along a cell-edge segment only), none wrong (`fossil_test_cut_lattice_pieces`). Solver support (issue #70, roadmap #68 WP1; design in `docs/guide/advanced/cut-solver.md`): `fossil_cut_graph` — `cut_links_t` (CSR of piece links: piece across, face 1-6 / 1-4, shared area, `unmatched`), `cut_links_2d(cut, links)` (equal fluid intervals across an edge; 2-D cut cells now give `piece_aperture`, `piece_iv_*`), `cut_cells_3d(..., links=)` (regions = the plane section's 2-D pieces, matched to the 3-D pieces through the shared fluid edge intervals; one-piece sides take all; islands by the unexplained face area). Edge fractions (#70): `cut_cells_3d(..., edges=.true.)` gives `cut%edge_x/_y/_z` (one `line_crossings` per grid line, crossings placed between the perturbed nodes by `line_crossing_before`; off by default, the swept geometry calls the cut cells many times); `cut%node_in` always (z lines; the 3-D corner states); 2-D `cut%node_in`; `curve%cut_cells(intervals=.false.)` skips `piece_iv_*` (the 3-D sections do unless `links=`). SRD (#70 step 5): `srd_neighbourhoods(links, piece_cell, piece_volume, piece_centroid, wet_normal, x, y, srd, z=, tau=)` → `srd_t` (CSR members with the piece first, `nov`, `vhat`, `xhat`, `eta_nb`, `eta_self`, `status` SRD_REGULAR/COMPLETE/INCOMPLETE/NO_VOLUME): growth through links inside the 3^d then 5^d tile, ordered by n.o/|o| (then |o|_1, lexicographic o, piece); zero-volume pieces never members. FRD (#70 step 6): `frd_neighbourhoods(links, piece_cell, piece_volume, x, y, frd, z=)` → `frd_t` (members reached from a cut piece by monotone index paths in the 3^d block, breadth first; weights kappa_j / sum kappa). Merge targets (#70 step 7): `merge_targets(links, piece_cell, piece_volume, wet_normal, x, y, mt, z=, tau=, mask=, exclude=)` → `merge_t` (single `target` by n.o, chains through unsuitable pieces to `final_target`, MERGE_OK/CHAIN/CYCLE/ALONE; normal-directed `mix_piece`/`mix_weight` at offsets o_a in {0, sign n_a}, beta weights after Hu et al. 2006, general knowledge). Moments (#70 step 8, `fossil_cut_moments`: `mono_count/_exponents`, `segment_green_moments`, `triangle_moments`, `cell_moments`, MOMENTS_MAX = 4): `cut_cells_3d(..., moments_order=)` → `piece_moments`, `piece_wall_moments`, `piece_face_moments(:, 6, np)`, `face_x/y/z_moments`, int X^a Y^b Z^c in scaled cell coordinates about the cell centre; volumes by divergence along x (walls + the two x faces; one-piece cells' faces from the plane sections); 2-D `curve%cut_cells(moments_order=)` → `piece_moments`, `piece_wall_moments`. Quadrature (#70 step 9): `cut_cells_3d(..., quadrature=.true.)` → `quad_start/_point/_weight` (signed cones from the cell centre over the oriented fluid boundary, 36 points per cone via `tet_rule`; one-piece cells' faces from the plane sections' `piece_seg`, multi-piece from the face arcs; fluid cells 27-point Gauss; `quad_status` 1 if the weights miss the volume), `wall_quad_*` (7-point rule on the wall fans); 2-D `curve%cut_cells(boundaries=, quadrature=)` → `piece_seg*`, `quad_*`, `wall_quad_*`. Coarsening (#70 step 10, `fossil_cut_coarsen`): `cut_coarsen_3d(fine, links, x, y, z, ratio, coarse, coarse_links, fine_to_coarse, status)` / `cut_coarsen_2d` — coarse pieces = components of fine pieces linked inside the coarse cell, everything else summed (moments re-centred), equal to the direct coarse cut cells. Performance: quiet 2-D cells (no wall piece, no crossing) take `quiet_cell` (one_cell's large work arrays are NaN-initialised in debug builds).
- `fossil_curvilinear_block.f90`, `fossil_block_geometry.f90`, `fossil_block_locate.f90`, `fossil_block_stencil.f90` —
  curvilinear structured blocks (issue #71, roadmap #68 WP2; design in `docs/guide/advanced/curvilinear-blocks.md`):
  `curvilinear_block_t` (vertices SoA `x, y, z(ni, nj, nk)`; cell (i, j, k) = vertices i..i+1, ...; faces 1-6 =
  -i, +i, -j, +j, -k, +k) and `curvilinear_block_2d_t` (`u, v`). Topology records `topo_record_t` (TOPO_BOUNDARY,
  TOPO_ABUT, TOPO_PERIODIC, TOPO_SINGULAR; CGNS transform, `record_map` valid off the face) proposed by
  `detect_topology(_2d)(blocks, tolerance, translations)` (exact coincidence by default; sets `blocks(b)%id`).
  Geometry: exact trilinear volumes (2-point Gauss, degree 2 per variable), face vectors 1/2 d1 x d2, centroids.
  Validity: `hex_validity` certifies sign(det J) by 27 Bernstein coefficients with error bounds + subdivision
  (CELL_VALID / SINGULAR (exact zeros from coincident vertices: poles) / INVALID / UNDETERMINED); 2-D exact by
  orient2d at the corners. Location: `build_search` then `locate_many` (cell-box BVH, 2-D exact orient2d containment +
  closed-form inverse, 3-D Newton within 1e-12; lowest linear index wins = low-side rule; LOC_INSIDE / OUTSIDE /
  WARPED / DEGENERATE / NOT_CONVERGED), `walk=.true.` (stencil jumping, equal to the robust path bit for bit on blocks
  that do not overlap themselves), `locate_many_device` (the `_oac` kernels on copies of the arrays; bench line
  `block location device vs serial mismatches` in both device benches), `locate_in_blocks(_2d)`. `stencils` (Lagrange
  width q <= 8 on the vertices, across records to the same block, shifted or reduced at physical boundaries,
  order-matched xi by Newton on the stencil's isoparametric map: q = 4 converges at 4, the trilinear xi caps it at 2).
  `quality` (equiangle skewness, orthogonality, aspect, warp, certified min J / mean J). Validation:
  `fossil_val_curvilinear_block` (spherical shell 128 x 64 x 32 with poles, C-grid 800 x 200).
- `fossil_overset.f90` — overset assembly kernels (issue #72, roadmap #68 WP3; design in `docs/guide/advanced/overset.md`),
  on sets of `curvilinear_block_t` / `_2d_t` (points = vertices; cell-centred data: the cell-centre blocks with a ghost
  or extrapolated layer). `wall_cutter(_2d)`: wall faces (TOPO_BOUNDARY named `wall`) -> surface/curve, shorter
  diagonal, solid inside, overlapping patches trimmed (a coarser patch's element dropped only if its vertices and
  centroid lie within the tolerance of the finer patch: no gap, a double-covered strip; `strict=.false.` = the
  centroid rule of #72, issue #91; `shared` walls). `hole_flags(_2d)`: FWN signed distance < offset (on the wall:
  field), CPU/device bit-exact, `cuts=` (exclude the posed owner blocks), `poses=`/`cutter_pose=`. `hole_map(_2d)`:
  Cartesian map, exact BOX_CUT, paint from corners, line of sight along a grid axis with exact crossing parity (exact;
  `pegasus=.true.` shows PEGASUS 5's thin-body failure); `xray_flags(_2d)`. `donor_candidates(_2d)`: all candidates per
  vertex (WP2 location per block pair; `device=` with block `enter_device` persistent state; `shifts=` wall
  projection; `poses=`; `previous=` walks), quality (donor/receptor volume, ratio, certified min J), wall distance;
  `overlap_matrix` (PAIR_NONE/OVERLAP/SHARED). `classify(_2d)`: SELECT_EXPLICIT (hole rings) / SELECT_IMPLICIT
  (ratio < 1), mandatory = faces named `overlap`, immune wall layers, stencil validity (no hole; explicit: no receptor,
  mandatory vertices count as receptors), Jacobi fixed point (block-order independent), circular pairs; iota 1 / 0 /
  -k, orphans ORPHAN_* by cause. `fringe_levels` (SUGGAR node rules), `overlap_width` (Chesshire's d),
  `wall_projection(_2d)` (shift to the donor wall, nudged 1e-9 into the fluid, linear decay), `exposed_points`
  (Meakin & Suhs violations), `connectivity_t` / `build_connectivity` / `save_connectivity`. Bench lines `hole flags`
  and `donor candidates device vs serial mismatches`. Validation `fossil_val_overset` (twilight-zone orders q through
  a cylinder assembly; 3-D sphere assembly timed). nvfortran 26.9: no self-referencing array growth
  (`v = reshape(v, ..., pad=v)` corrupted the heap).
- `fossil_implicit_points.f90` — implicit points and exact predicates on them (issue #65, the boolean's exact
  arrangement): `implicit_point_t` of kind `IPK_EXPLICIT` (input vertex), `IPK_LPI` (line p,q meets plane r,s,t),
  `IPK_SSI` (coplanar segments meeting, solved in the projection dropping `drop`), canonical input order (one
  construction = one identity, `ip_same`; `ip_approx` rounds by one formula: equal identities, equal bits);
  `ip_orient3d` (implicit vs explicit plane), `ip_orient2d` (three implicit points in a projection), `ip_compare`
  (coordinate order): homogeneous (X, W) polynomials, a running-error filter in doubles, then exact expansion
  arithmetic (18-bit bit-mask splits, FMA-safe). Golden signs in `src/tests/implicit_points_golden.txt` (exact
  rational reference). Also `IPK_EEI` (two exactly coplanar lines), `ip_lex_compare` (exact lexicographic order),
  `ip_identity_compare`, `ip_rank`. Also `IPK_TPI` (issue #91: three planes through three input points each, Cramer about a1,
  W degree 6, X degree 7; identity = sorted plane triple; rank 5); stage-2 intermediates are compressed (Shewchuk
  Compress: signs only, the older kinds round to the same bits). `fossil_arrangement%add_triple_points`: two cuts of
  one facet made by distinct transversal facets that cross inside it get their TPI in all three facets. `fossil_cdt_implicit%cdt_facet`: constrained triangulation of a facet on implicit
  points (insertion by split, Sloan recovery: exact orient2d for validity; Lawson incircle on rounded coordinates for
  quality only).
- Named patches (issue #66, roadmap #58 WP6): `facet_object%patch` (default 1) indexes the surface's private
  `patch_names(:)` (default [header]); `sync_patch_names` (called by `adopt_facets`) extends the table to every id.
  The ASCII reader opens a patch per `solid <name>` block (`stl_ascii_next_facet`; before #66 a multi-solid file
  aborted at the first `endsolid`), `load_from_files` one per file, the ASCII writer one block per patch (a
  single patch writes the old file). API: `get_patches_number`, `patch_name`, `facet_patch`, `set_patch`.
  The facet id travels through every copy; `merge_solids`, `boolean` merge tables by name; `fill_holes` uses the
  boundary's patch or `filled`; remesh/decimate/smooth lock patch boundaries; `alpha_wrap` takes the nearest input
  facet's. Queries: `query_many(_device)(patches=)`, `probe_many`, `image_points` (device: host lookup from the
  returned facet ids; bench line `patch device vs serial mismatches`).
  2-D: `curve_object` (private `seg_patch`, `patch_names`): `# patch <name>` lines (do not end a loop), `repair`,
  `slice_surface` (facet's patch), `set_patch`, queries `patches=` (curve bench line too). Empty patches write no block.
  `fossilizer --patches prefix_` lists them and writes one STL per patch.
- `fossil_body_distance.f90` — body-body distance (issue #67, roadmap #58 WP7): `tri_tri_distance` (0 when `pair_cut`
  says the triangles meet, coplanar overlap by exact orient2d; else min of 6 vertex-triangle + 9 edge-edge; kernels stable to ~eps L since issue #69: `closest_segment_segment` by projection on the plane normal to the longer segment, `closest_on_triangle` = `triangle_closest_st`),
  `body_distance` (best-first descent of two AABB trees, heap by box distance, tie rule smallest (facet_a, facet_b):
  bit-identical to `body_distance_brute`), `surface%distance_to(other, ..., pose_self, pose_other)` (poses: facets
  moved by `to_world`, node bounds = box of the 8 moved corners; equal to `move` bit for bit), `proximity_pairs(cutoff)`
  over `facet_body` bodies and the free `proximity_pairs_of(surfaces, cutoff)` (box broad phase). Curve twin
  `curve%distance_to` (exact orient2d contact, box-filtered loop over segment pairs, same tie rule).
  `fossilizer --distance other.stl`. CPU only.
- `fossil_ccd.f90` — continuous collision detection (issue #73, roadmap #68 WP4; design in
  `docs/guide/advanced/collision.md`): the motion of `pose_interpolate` (rotation about the axis through the BODY-FRAME
  origin, angle and translation linear in s) has the exact speed bound `motion_bound_3d/_2d` (|dt| + |theta| max_v rho_v,
  rho = distance to that axis); `ccd_advance` (Mirtich's conservative advancement over an abstract `ccd_evaluator_t`:
  steps (d - tol/2) / mu, d >= tol/2 on [0, toi], d(toi) <= tol, no tunnelling; CCD_CONTACT / SEPARATE /
  INITIAL_CONTACT / NOT_CONVERGED (cap 1000, toi a certified lower bound) / INVALID_INPUT). `surface%time_of_impact`,
  `curve%time_of_impact` (over `distance_to` with poses; default tol 1e-6 of the larger bbox diagonal);
  `directional=.true.` adds C2A's per-pair separating-plane step (`body_safe_step` in `fossil_body_distance`; valid per
  convex triangle / segment pair only, never for a whole non-convex body) and takes the larger certified step;
  `.false.` plain only; absent = adaptive (issue #75): plain until 8 consecutive steps stall (< 1 % of 1 - s AND < 1 %
  of the margin d - tol/2 closed), then the larger step; never-stalling loops are the plain loop bit for bit. `first_impact` / `first_impact_2d`: swept-box broad phase, candidates by lower bound, equal to the loop over all
  pairs. CPU only. Validation `fossil_val_ccd`.
- `fossil_projection.f90`, `fossil_block_report.f90`, `fossil_hyperbolic.f90` — body-fitted grid support (issue #74,
  roadmap #68 WP5; design in `docs/guide/advanced/body-fitted.md`): `projection_target_t` (`surface%projection_target`,
  `curve%projection_target`: elements by patch with one tree per patch; constraint curves = feature edges of #61 +
  patch boundaries, chained between corners; corners); `project_points(_2d)` with PROJ_FREE / PATCH / CURVE / CORNER
  (ties to the lowest index; corners bit for bit), `classify_points(_2d)` (corner, else curve, else the patch of the
  closest facet, within a tolerance), `project_block_face` / `project_block_side` (in place). `block_report(_2d)`
  (own module: fossil_projection cannot use the surface type): y1 along the grid line and to the wall cutter,
  y1+ (u_tau given, or `utau_flat_plate` flagged as an estimate), growth ratio, WP2 quality, `cell_centred` (wall
  cells), `save_block_report` (text). `hyperbolic_grid_2d`: normal steps of mitre length h/cos(beta/2) (exact polygon
  offsets), implicit smoothing of the directions; folds at right-angled concave notches are reported (status 2), not
  prevented; `geometric_spacing`. CPU only. Validation `fossil_val_body_fitted`.
- `fossil_box_predicates.f90`, `fossil_box_query.f90` — AMR tagging (issue #38): box-triangle separating-axis test (conservative: run on the box inflated by `box_tolerance`, 1e-12 of the problem size, never misses a contact), exact 2-D box-segment test; `bvh_box_contact_oac` / `curve_box_contact_oac` (shared CPU/device traversals), `cell_bands` / `cell_components` (tag-band BFS, 26-connected components). The device box path follows the #20 rules and needs no new device arrays.

## Coding conventions

- Fortran 2003+ OOP: all types use `type :: foo_object` with `contains` bound procedures
- All bound procedures use `pass(self)` explicitly
- `destroy` / `initialize` pattern: objects are reset with a fresh local instance (`self = fresh`)
- Numeric types from PENF (`I4P`, `R8P`, etc.) — never use bare `integer` or `real` for library code
- 3D geometry via VecFor `vector_R8P` — not plain arrays

## Known limitations

- **`surface%boolean` / `surface%resolve_self_intersections`**: exact
  arrangement since issue #65, no tolerance anywhere: cut points are
  implicit points (`fossil_implicit_points`), pairs are classified by
  `pair_cut` / `coplanar_cut` (`fossil_tritri_exact`), the registry merges
  points exactly and gives distinct points distinct rounded coordinates,
  facets are triangulated by `cdt_facet` (`fossil_cdt_implicit`, exact
  orient2d), and the selection goes by patches (exact decision at a cut
  edge, `local_inside`; winding number only for a patch without cuts).
  Faces coplanar by construction but off by round-off (rotated) cross
  exactly: the result keeps the slivers, watertight with the right volume,
  and its rounding to doubles can fold a sliver by an ulp (self-intersecting
  pairs; the harness allows them only up to a depth of 1e-12 of the diagonal,
  measured <= 1e-15). Opt-in `boolean(..., weld=tol)`: B's vertices within
  tol of A's take A's coordinates bit for bit (unique nearest, never two
  onto one, B kept non-degenerate; `weld_vertices` in `fossil_boolean`),
  which merges inputs that share vertices up to round-off; default off.
  Contact-only and touching-symdiff results are non-manifold by definition.
  `resolve_self_intersections` (issue #78): the arrangement of the surface
  with itself from `find_self_intersections`' pairs, then Zhou et al.'s
  selection (patches with winding number 0 in front, from a rounded
  solid-angle sum, O(patches x facets)); `BOOL_STATUS_NOT_CLOSED` for an
  open selection. Three facets meeting inside a facet get a three-plane
  point (`IPK_TPI`, issue #91): the dragon (open) gives `BOOL_STATUS_NOT_CLOSED`;
  still `CDT_FAILED`: a cut along a coplanar partner's edge crossing another cut
  inside a facet (inferred, not triggered).
  Validated by identities that need no analytic answer (test-side harness
  `src/tests/fossil_boolean_harness.f90`: volume, second moments, area,
  deviation, membership, validity; self-intersections by the exact
  `find_self_intersections` since issue #69): curved bodies, the dense degenerate sweep
  (`fossil_val_boolean_sweep`, 360 cases, 1e-4..1e-14), CSG chains.
- **`find_self_intersections`**: every decision exact since issue #69 (`pair_cut` narrow phase, coplanar overlaps
  in the projection, exact coplanar pairs at a shared vertex, folds over a shared edge, sheets crossing along a
  non-manifold edge by the exact fan order); duplicate facets (three shared vertices) are left to
  `remove_duplicate_facets`.

<!-- rtk-instructions v2 -->
# RTK (Rust Token Killer) - Token-Optimized Commands

## Golden Rule

**Always prefix commands with `rtk`**. If RTK has a dedicated filter, it uses it. If not, it passes through unchanged. This means RTK is always safe to use.

**Important**: Even in command chains with `&&`, use `rtk`:
```bash
# ❌ Wrong
git add . && git commit -m "msg" && git push

# ✅ Correct
rtk git add . && rtk git commit -m "msg" && rtk git push
```

## RTK Commands by Workflow

### Build & Compile (80-90% savings)
```bash
rtk cargo build         # Cargo build output
rtk cargo check         # Cargo check output
rtk cargo clippy        # Clippy warnings grouped by file (80%)
rtk tsc                 # TypeScript errors grouped by file/code (83%)
rtk lint                # ESLint/Biome violations grouped (84%)
rtk prettier --check    # Files needing format only (70%)
rtk next build          # Next.js build with route metrics (87%)
```

### Test (90-99% savings)
```bash
rtk cargo test          # Cargo test failures only (90%)
rtk vitest run          # Vitest failures only (99.5%)
rtk playwright test     # Playwright failures only (94%)
rtk test <cmd>          # Generic test wrapper - failures only
```

### Git (59-80% savings)
```bash
rtk git status          # Compact status
rtk git log             # Compact log (works with all git flags)
rtk git diff            # Compact diff (80%)
rtk git show            # Compact show (80%)
rtk git add             # Ultra-compact confirmations (59%)
rtk git commit          # Ultra-compact confirmations (59%)
rtk git push            # Ultra-compact confirmations
rtk git pull            # Ultra-compact confirmations
rtk git branch          # Compact branch list
rtk git fetch           # Compact fetch
rtk git stash           # Compact stash
rtk git worktree        # Compact worktree
```

Note: Git passthrough works for ALL subcommands, even those not explicitly listed.

### GitHub (26-87% savings)
```bash
rtk gh pr view <num>    # Compact PR view (87%)
rtk gh pr checks        # Compact PR checks (79%)
rtk gh run list         # Compact workflow runs (82%)
rtk gh issue list       # Compact issue list (80%)
rtk gh api              # Compact API responses (26%)
```

### JavaScript/TypeScript Tooling (70-90% savings)
```bash
rtk pnpm list           # Compact dependency tree (70%)
rtk pnpm outdated       # Compact outdated packages (80%)
rtk pnpm install        # Compact install output (90%)
rtk npm run <script>    # Compact npm script output
rtk npx <cmd>           # Compact npx command output
rtk prisma              # Prisma without ASCII art (88%)
```

### Files & Search (60-75% savings)
```bash
rtk ls <path>           # Tree format, compact (65%)
rtk read <file>         # Code reading with filtering (60%)
rtk grep <pattern>      # Search grouped by file (75%)
rtk find <pattern>      # Find grouped by directory (70%)
```

### Analysis & Debug (70-90% savings)
```bash
rtk err <cmd>           # Filter errors only from any command
rtk log <file>          # Deduplicated logs with counts
rtk json <file>         # JSON structure without values
rtk deps                # Dependency overview
rtk env                 # Environment variables compact
rtk summary <cmd>       # Smart summary of command output
rtk diff                # Ultra-compact diffs
```

### Infrastructure (85% savings)
```bash
rtk docker ps           # Compact container list
rtk docker images       # Compact image list
rtk docker logs <c>     # Deduplicated logs
rtk kubectl get         # Compact resource list
rtk kubectl logs        # Deduplicated pod logs
```

### Network (65-70% savings)
```bash
rtk curl <url>          # Compact HTTP responses (70%)
rtk wget <url>          # Compact download output (65%)
```

### Meta Commands
```bash
rtk gain                # View token savings statistics
rtk gain --history      # View command history with savings
rtk discover            # Analyze Claude Code sessions for missed RTK usage
rtk proxy <cmd>         # Run command without filtering (for debugging)
rtk init                # Add RTK instructions to CLAUDE.md
rtk init --global       # Add RTK to ~/.claude/CLAUDE.md
```

## Token Savings Overview

| Category | Commands | Typical Savings |
|----------|----------|-----------------|
| Tests | vitest, playwright, cargo test | 90-99% |
| Build | next, tsc, lint, prettier | 70-87% |
| Git | status, log, diff, add, commit | 59-80% |
| GitHub | gh pr, gh run, gh issue | 26-87% |
| Package Managers | pnpm, npm, npx | 70-90% |
| Files | ls, read, grep, find | 60-75% |
| Infrastructure | docker, kubectl | 85% |
| Network | curl, wget | 65-70% |

Overall average: **60-90% token reduction** on common development operations.
<!-- /rtk-instructions -->