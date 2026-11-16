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
- No `norm2` intrinsic in device code, and above all none of an array constructor: nvfortran 26.9's device
  `norm2([x, y])` returns 0 (`-O0` and `-O2`; `sum`/`maxval` of the same constructor are right), and the device
  `norm2` of an array is 1 ulp off the host's. Use an explicit replica (`norm2_xy` in `fossil_cut_cells_2d`).
  Reproducers: issue #101, items 12 and 9a.

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
`device vs serial max |d| err 0.000E+00`; the `rich`, `box`, `posed`, `sign bench`, `line crossings`, `deviation`, `feature`, `thin`, `cut cells 2-D` and `cut cells 2-D chain` mismatch lines must be 0 too (the 3-D bench prints the same lines, the 2-D cut cells excepted, plus `orient3d`, `cut cells 3-D sections`, `cut cells 3-D nodes`, `cut cells 3-D` and `swept cells 3-D`, `cut cells 3-D link inputs`, `cut links 3-D`, `cut walls 3-D` and `piece map 3-D`; the curve bench also `cut cells 2-D blocks grids`, `swept cells 2-D`, `cut links 2-D` and `piece map 2-D`; both print `probe`, WP1 step 01, `curvature`, WP1 step 03, and `closure`, WP1 step 04). Validated on the GPU 2026-11-05 (26.9, RTX 4070, 512²).

The `tests-nvidia-chain` mode (WP2 #100 step 6b) builds `src/tests/validation/fossil_val_device_chain.f90` only: the
per-step chain of a moving-body cut-cell solver (move, tag, cut + links, swept, and on the Gauss path the piece map
since #105 step 4; 3-D sphere and cone, 2-D ellipse and NACA 0012), every stage's device output compared bitwise with
the CPU (pass criterion `device chain vs serial mismatches: 0`; the map by `fossil_piece_map_compare`), then timed on
the CPU, the device copied back and the device resident (the map column = the median of the paired differences of the swept
call with `map=` and without, with their spread, `unresolved` when not above it; the map has no resident form). The default run is the gate (3-D 32³/64³, 2-D 256²; the map checked at every
size, timed at 32³/256²; 7-10 min on GPU 1); `fossil_val_device_chain full` adds the step-0 sizes 128³ and 1024² and
times the map everywhere (about 35 min before the map stage). Built with gfortran it is a validation program: the device routines run as host twins at 8³/48².

The NVIDIA modes compile with `-gpu=nofma -Mnofma` (`$CSTATIC_NV`): with FMA contraction the
device lost host bit-exactness (1-ulp distances, sign and facet-id flips at near-ties); the device
kernels are 15-30 % slower for it. Do not remove it unless the gates stay at 0. The GNU modes have the twin guard
`-ffp-contract=off` (`$CSTATIC_GNU`, `$CSHARED_GNU`): a gfortran targeting FMA hardware (`-march=x86-64-v3`)
contracted at -O2 and broke the bitwise twin checks (296 chain mismatches without it, 0 with it; outputs and perf
hashes on baseline x86-64 unchanged). VecFor's explicit `fma` calls stay FMAs (correctly rounded either way).

Subnormals: the nvfortran host code flushes subnormal results to zero (FTZ) and treats subnormal operands as zero
(DAZ); the device does neither, so a result below ~2.2e-308 can differ host vs device. Host/device (and gfortran)
agreement there needs both `-Mnoflushz -Mnodaz` on the host: `-Mnoflushz` alone is not enough. Neither is in
`$CSTATIC_NV` today (an open owner decision; issue #101, item 9b).

The device build scope is **BVH-only** (the octree stays CPU-only). The 2-D and 3-D cut cells also run on the device (`curve%cut_cells_device`, `cut_cells_3d_device`, WP2 #100), with their piece links (`links=`, `cut_links_2d_device`) and the swept geometry (`swept_cells_3d_device`, `swept_cells_2d_device`): the whole per-step chain of a moving body (move, tag, cut + links, swept), bit for bit the CPU's.

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
./scripts/run_validation.sh -j 4     # 4 programs at a time (single-threaded, no files written); results in name
                                     # order, a failing test prints its log tail (whole log: exe-validation/NAME.log)
```

The memory-heavy programs (`SOLO` in the script: `fossil_test_distance_bench` 11.9 GB peak resident,
`fossil_val_cut_solver` 5.0 GB; every other one under 1.6 GB) run first, alone: a GitHub runner has 16 GB.

`fossil_test_distance_bench` prints its five timing tables of the 3-D device kernels only when built with OpenACC (or
with `FOSSIL_BENCH_TIMING=1`): without it they time the host twins against the CPU (gfortran, dragon-fine 40^3: 1375 s
with them, 384 s without). The nightly `validation.yml` (project-owned) calls fortran-toolchains' `matrix.yml@v1` on
`linux-gnu`: `scripts/ci/validation.sh` builds `validation-gnu` and runs `-j $(nproc)`; `exe-validation/*.log` is
uploaded as the artifact `validation-logs-<toolchain>` on every run.

CI toolchain: the scaffold workflows are thin callers of the reusable workflows of `szaghi/fortran-toolchains` pinned
at `@v1` (repo `~/fortran/fortran-toolchains`), rendered from `[scaffold] toolchain = linux-gnu` in `fobos`: `ci.yml`
(tests, coverage) and `docs.yml` run in its `gnu-16.2` image (gfortran 16.2 on the official `gcc` image, Open MPI
5.0.11, HDF5 2.3.0, FoBiS, formal, gcovr; nothing is installed per job), as the project-owned `docs-snippets.yml`
and `validation.yml` (through the reusable `matrix.yml`); `release.yml` publishes a tag's release, then `install.yml` installs it on a fresh runner (the
linux-bare toolchain). The ubuntu-26.04 runners install amd64v3 packages: their gfortran emitted FMAs at
-O2 and broke the bitwise device-twin checks of `fossil_val_device_chain` / `fossil_test_distance_bench` (run
43194838880); the image's gfortran targets baseline x86-64 (no FMA), and `-ffp-contract=off` guards any other. Every `tests-*` mode excludes `src/tests/validation/`. `fobis clean` leaves stale executables of
moved sources in `exe/`: delete them by hand, or `run_tests.sh` runs them.

Performance (issue #94): `fossil_val_perf` (validation) times 24 kernels (query batches, body distance, cut cells,
swept cells, CCD, remesh/decimate, overset), each line with t_min / t_median of N repetitions after a warm-up and a
64-bit hash of the outputs' bits. `scripts/run_perf.sh NAME` writes the named baseline
`src/tests/validation/perf/NAME.txt` (one thread; an existing name is refused: baselines are never overwritten);
`scripts/run_perf.sh --compare A B` prints the speed-ups and fails on any hash difference. A performance change must
keep every hash (bit-identical outputs). Baselines are machine-specific (the header names host, compiler, flags,
commit): compare runs from one machine. `run_validation.sh` runs the bench with one repetition.

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
`check_doc_snippets.sh`). Brand: `docs/brand/marks.py` generates `docs/public/{logo,logo-dark,favicon,favicon-dark}.svg` (the Sumi ammonite; PNG favicons rasterised from `favicon.svg`) and the 18 landing-card marks `docs/public/icons/<slug>{,-dark}.svg` (same brush, drawn for 48 px; the `icon:` objects of `docs/index.md`); never hand-edit the SVGs;
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
5. Queries: `distance`, `is_point_inside`, plus `compute_distance` (richer form with closest facet/edge/vertex, closest point and normal outputs). IB queries (issue #37): `query_point` / `query_many` / `query_many_device` (closest point, facet, `FEATURE_FACE/EDGE/VERTEX`, pseudo-normal, body; optional `ambiguous` / `second_foot_points` on the CPU), `probe_many` / `probe_many_device`, `image_points` / `image_points_device` (WP1 step 01: `lengths(m, n)` -> `probe_points(3, m, n)`, `visible(m, n)`, a generic over the one-length form; `pose=` on host and device; visibility by `fossil_segment_visibility`, below); `query_many(_device)(principal_curvatures=, principal_directions=)` after `build_curvature` (WP1 step 03, `fossil_shape_operator` below); shells and bodies from `analyze` (`facet_shell`, `facet_body`; `body_payload` is device-mapped). Curves have the same set. `distance_field_object%compute(store_closest=.true.)` keeps `closest`, `closest_id`, `closest_normal` on band cells. AMR tagging (issue #38): `box_classify_many` / `box_classify_many_device` (`BOX_CUT/NEAR/INSIDE/OUTSIDE`), `box_intersects_surface`, `tag_block` (signed Chebyshev tag bands on a block); curve twins. These default to `SIGN_FAST_WINDING_NUMBER` (issue #39); the point queries keep `SIGN_PSEUDO_NORMAL`.
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
  uses it. The per-ray kernels (`ray_prepare_oac`, `ray_crossing_oac`, side tests) live in `fossil_ray_kernel.f90` on
  plain arrays, `!$acc routine seq` (WP1 step 01): `ray_prepare` / `ray_facet_crossing` call them, bit for bit the
  pre-split code (sums written as gfortran evaluates `sum` / `dot_product`: from 0, left to right).
- `fossil_segment_visibility.f90` — probe/image visibility (WP1 step 01), one any-hit kernel per dimension shared by the
  host `probe_many` / `image_points` and their `_device` twins: `bvh_segment_blocked_oac` (SAH BVH, the exact ray-facet
  crossing of `fossil_ray_kernel` on the ray `intersect_ray_all` prepares against the root box, hit at tol < t <= |b - a|,
  tol = 1e-9 bbox diagonal; nodes pruned by the padded slab test and by the least projection of the padded box beyond
  |b - a| + pad; equal to the former `intersect_ray_all` + filter: 0 differences on 1 200 000 probes/images of nine surfaces and three curves (validation `fossil_val_probe_visibility`),
  and in `fossil_test_probe_device_twin`); `curve_segment_blocked_oac` (segment BVH, the `segments_in_box` payload filter
  then exact `segments_intersect_2d`, the foot's segments skipped). An octree surface keeps `segment_clear` (host ray
  query). The device probes take feet and normals from `query_many_device`, then one loop over (probe, point) pairs;
  CI twins `fossil_test_probe_device_twin`, `_2d`.
- `fossil_block_object.f90` / `fossil_block_aabb_object.f90` — structured mesh block with AABB (`extrapolate_distances` deprecated, issue #32)
- `fossil_distance_field.f90` — `distance_field_object` (re-exported by `fossil`, issue #32): signed-distance field on a rectilinear grid; narrow band via `distance_many(_device)`, flood-fill signs, fast sweeping in hyperplane order (OpenMP bit-exact). It imports the surface/facet modules, never `fossil` (circular). `fossil_ib_generator` is built on it. `initialize` without `z` makes a 2-D field for a curve source; `compute_band` / `compute_signs` / `compute` / `default_band` are generics over `surface=` / `curve=`.
- `fossil_curve_object.f90`, `fossil_segment_bvh.f90`, `fossil_curve_slice.f90` — 2-D geometry (issue #36): `curve_object` stores (u, v) as SoA arrays in the plane normal to `normal_axis` (u, v = the other axes in cyclic order); the solid is on the left of every segment (outer loops CCW, holes CW), signed distance negative inside; open loops allowed (`# open` in the ASCII format). The segment BVH's traversals are `_oac` kernels shared by the CPU and the device. `slice_surface` cuts an STL by an axis-aligned plane (canonical crossings, no tolerances). `fossil_predicates%orient2d` and `orient3d` (`orient3d_xyz`) are exact (filter + exact expansion with a bit-mask split, safe under nvfortran FMA contraction; `orient2d_xy` returns +0 right after the filter when two inputs coincide, the exact stage's value (#100 step 3a); `grow_expansion` works in place, no automatic array; orient3d, a device routine since #100 step 3c (`o3d_*` module helpers, expansion `e(96)`: at most 24 four-double products; ptxas 1064 B): exact zero on a shared coordinate plane, 6 triple products of exact differences or 24 of input coordinates, each exact as four doubles by error-free products (Dekker, bit-mask split; no `ieee_fma`, which nvfortran 26.9 lacks) since issue #65; no `predicates_initialize` needed).
- Issue #42 (WP6): `distance_field_object%gradient_diagnostic` (eikonal q per node), `wall_distance_many(_device)` (BVH started at cutoff^2; curve twin). `query_many(curvatures=)` (foot curvature: (H, K) from `jet_curvature` — cubic fit over the 2-ring, O(h^2) pointwise; the #18 cotangent / angle-defect operators stall near 15 % at valence-5 vertices — by barycentric weights; 2-D turning angle; CPU only), `surface_markers` (remesh to sqrt(3) dx + one calibration pass) / `curve_markers` (exact arc-length resampling). `isotropic_remesh(preserve_features)` re-detects features after each split (else crease midpoints were relaxed off and edges chamfered); since roadmap #92 WP5 (c) (only with `preserve_features`; the default path bit for bit): `crease_runs` cuts the feature lines (`find_feature_edges`) into runs between fixed vertices (corners, non-manifold ends, patch-boundary vertices when patched); a straight run (interior within `STRAIGHT_TOL` = 1e-9 of its length from its chord, no non-manifold edge) collapses along its line in the 4L/5..4L/3 window through `try_collapse_edge` (onto the fixed end, or the midpoint first) and its vertices relax along the chord (`apply_crease_relax`, clamped between the run neighbours: axis-aligned coordinates bit for bit); curved runs and patch boundaries stay locked; each iteration splits every edge above 4L/3 in up to `SPLIT_ROUNDS` = 8 rounds after the high-water slot (one pass starved: 144 of 149 crease collapses refused by `high`); no collapse or slide may leave a facet below `SLIVER_FRACTION` = 0.1 of its area vector (`facet_shrinks`; the flip test accepted folds flat onto three crease vertices: tutorial 12 bracket min area 5e-12 -> 1.96e-4). cube-01 at L = 0.2: 338 facets, edges 0.146..0.264 (406, 0.095..0.335); at L = 1: 20 facets (36). Markers across a 90-degree crease are 0.707 dx apart by construction (ideal cube fraction ~1 - 3 dx = 0.70 at dx 0.1): cube 2 x 2 at dx 0.1 0.528 -> 0.726 (`fossil_test_surface_markers`, 7-8 vertices per cube edge, spacings 0.807..1.019 L); sphere bit for bit; SLIVER_FRACTION and SPLIT_ROUNDS not swept. 12 docs outputs moved (tut12-*, ref_isotropic_remesh, ref_patches, ref_ib_queries, ck_markers, ck_remesh). Since issue #79 the collapse pass coarsens (link condition, one collapse per vertex per pass, no edge above 4L/3, locked endpoint kept); since issue #91 the relaxation is tangential and `REMESH_ERROR_TOL` = 0.07 L refuses any collapse, flip or move whose new facets leave the input farther (cube-01 at L = 1 without features keeps 0.949 of its volume). `decimate` (issue #81): link condition and zero-area refusal; `target_facets <= 0` is BAD_INPUT; patch rims coarsen along themselves onto input rim vertices (`rim_collapses`, default on, corners fixed; issue #91); optional `max_error`: h(result, input) <= max_error certified to round-off (roadmap #92 WP5 (b): `cert_collapse` checks the target and every new one-ring facet by `cert_triangle`, a branch and bound against the input's AABB tree (wrapper passes `reference_facet` / `reference_tree`, required with max_error): least of Lipschitz max d_k + rho_T, max corner distance to one input facet, and the same over a split along the edge between two input facets (`edge_bound`) or the wedges of an input vertex's fan (`wedge_bound`, `fan_build`), else split in four; caps 40 levels / 4096 triangles refuse as undecided, never reached); a refused collapse is skipped, its edge re-queued when a neighbour moves an end; the plane-cost stop (sqrt cost > max_error -> NO_PROGRESS) is kept as the only guard of the input -> result side, which is not certified (certificate alone: dragon 1e-2 diag 2960 -> 1020 facets, that side 3-5 max_error). Dragon 3e-4 diag: [1.70, 1.74] -> [0.83, 0.87] max_error (6534 -> 6538 facets, 2 refusals), bunny 3e-3: [1.11, 1.13] -> [0.53, 0.58]; smooth bodies 0 refusals (every docs output and the bracket validation unchanged); +15-60 % on the dragon, +6 % on the bunny (`find_shared_facets` O(n) per collapse dominates there); `fossil_test_decimate_certified`. `extract_isosurface` (issue #87) snaps a node to the level within 1e-12 of its edge crossings (round-off slivers). `jet_curvature` merges near-coincident ring vertices and returns NaN + `CURV_STATUS_FIT_FAILED` when the fit stays rank-deficient (issue #91). `fill_holes(manifold_only=)` / `sanitize(manifold_holes=)` chain loops through true boundary edges only (opt-in).
- `fossil_rigid_pose.f90` — moving bodies (issue #41): `rigid_pose_t` / `rigid_pose_2d_t`; optional `pose=` on the point queries (body-frame, CPU and device), `move(pose)` with `aabb_tree%refit` / `segment_bvh%refit_vertices` (O(N); `aabb_node%set_box` writes the box only — `node%initialize` wipes children), `state_changes`; `fossil_swept_cells.f90`: `swept_cells_2d(curve, ...)` / `swept_cells_3d(surface, ...)` (free routines, not TBPs: the module uses curve/surface/cut_cells_3d) — swept volume = cut-cell volume change (exact), face apertures time-averaged (translation exact via breakpoints: 2-D per face, 3-D per grid plane — mesh vertex on the moving plane, mesh edge on a moving grid line, grid node on a triangle; area quadratic per piece, 2-point Gauss; `exact=.false.` forces Gauss. Rotation: 2-D exact since issue #69, by local integration since #104 (`faces_rotation_local_2d`: the piece map's machinery on every candidate cell, `all_cand` + `all_meas`, each integrated per event interval by the doubling rounds to round-off on its local wall; faces of no candidate cell closed form, `wall_free_state`; untrusted / failed cells, anomalies and open curves through `faces_rotation_2d` on a face mask: certified breakpoints via `isolate_zeros` on g(s) = (p0+p1 s) cos phi + (q0+q1 s) sin phi + k0 + k1 s with the |g''| bound, Gauss doubled to round-off per piece; equal to the old path to round-off, 36-286x faster up to 2 cells per step, 0.07-0.65 of a Gauss call, 7.5-9x at 20 degrees; `swept_cells_2d_reference` = the old path, for tests; `sw%evaluations` there = the largest node count of one cell; `fallback_faces`; optional `outflow_u/_v`; the map's device passes on `swept_cells_2d_device`); 3-D exact too (`faces_rotation_3d`, per grid plane: vertex on plane, edge on grid line, node on triangle, all of the same g form by Rodrigues; zero zones and 1e-10 merge in `isolate_zeros`); Gauss-Legendre `nq` with `exact=.false.` or explicit poses; `error_estimate` = max(nq-vs-2nq face change, max |residual| / largest rigid-field speed) on the screw Gauss path since issue #104; `refine=.true.` (screw Gauss path, opt-in, #104): `piece_map_*_build(all_meas=, measured=)` integrates every event cell per event interval on its local wall, `refine_faces_2d/3d` writes those faces (the low cell's value on a shared face), residual to round-off, 1.1-2.2x Gauss in 2-D and 2.7-22x in 3-D; the device path then cuts the end poses in full (host in 2-D, `cut_cells_3d_device` without `alpha_only` in 3-D) and refines on the host), `wall` = -sum A_f n_f, `residual` = rigid-field volume balance; `pose_interpolate` (screw motion), `gauss_legendre_01`. Performance (issue #94): the exact 3-D rotation slices the body-frame surface with `slice_surface(pose=)` (each vertex placed by `to_world` as read, equal to slicing a moved copy bit for bit) instead of copy + `move` per Gauss node; the 2-D rotation places one work copy per node with `place_for_crossings` (vertices, bbox, node boxes via `segment_bvh%refit_boxes`, the rest stale by design) and memoises the crossing events (wall vertex on a grid line once per line, node on a segment once per node). Device (WP2 #100 step 4): `swept_cells_3d_device` (free; shares `swept_3d_run` with the host path) = `swept_cells_3d` bit for bit on every path: end poses by `cut_cells_3d_device(alpha_only)`; faces batched: the host builds the (plane, node) list as `plane_average` / `plane_piece` (breakpoints, events, `isolate_zeros` on the host: <= 1 % of a call), slices through the plane sweep (posed: `slice_planes_prepare(spans=)` from the step's vertex paths with a margin, `slice_surface_plane(pose=)`), one `cut_cells_device_many` per chunk on per-block grids, then one thread per face adds over its nodes in node order (`translation_reduce_dev`, `rotation_reduce_dev`; no atomics; the rotation's Gauss doubling in rounds with the host's decisions); Gauss path via `cut_cells_3d_sections(pose=)` (no copy + move per node); public optional `chunk_cells` (default 2^22) caps the launch memory without changing a bit. 1.7-4.7x on an RTX 4070 at 32^3-64^3; the exact rotation at 64^3 is 250-600 s per host call (650-15 400 sections per plane). `swept_cells_2d_device`: only the end poses on the device (the 2-D faces are 94-99.9 % of a call: `collect`, the node breakpoints, `face_measure`, `place_for_crossings`; gains left are host algorithms): 1.00-1.04x; the exact rotation (local integration, #104) runs the map's per-item passes on the device (`faces_rotation_local_2d(device=)`), 1.5x / 1.9x at 256^2 / 1024^2. Bench lines `swept cells 3-D` / `2-D device vs serial mismatches`; CI twins `fossil_test_swept_3d_device_twin`, `fossil_test_swept_2d_device_twin`.
- `fossil_piece_map.f90` — the piece map of a time step (WP3 #102 of roadmap #92): `swept_cells_2d(..., map=)` / `swept_cells_3d(..., map=)` (refused with `poses=`) return `piece_map_2d_t` / `piece_map_3d_t`; `swept_cells_2d_device(..., map=)` / `swept_cells_3d_device(..., map=)` the same maps bit for bit (#105, below).
  - **Tracks:** the components of a cell's space-time fluid region over the step. Each pose0 piece and each pose1 piece lies in one track: TRACK_CONTINUING / FRESH / DEAD / SPLIT / MERGE / SPLIT_MERGE / TRANSIENT, or UNTRUSTED (-1: one track with all the cell's pieces and its values).
  - **Layout:** CSR cell -> tracks (`tracks`, `first_track`), old and new pieces (`old_start/old_piece`, `new_start/new_piece`, numbered as the end-pose cut cells number them), inverse maps `old_track` / `new_track`.
  - **Per track:** `volume0/1`, `faces` (4 or 6: time-averaged fluid length or area), `wall`, `outflow` (the step's rigid-field outflow per face part; under rotation the velocity varies along a face, so a conservative update needs it) and `residual` = V1 - V0 + sum outflow.
  - **Construction:** a cell's topology changes only at the events of its boundary, the swept breakpoints: 2-D wall vertex on a face line, corner on a segment; 3-D mesh vertex on a face plane, mesh edge on a cell edge, corner on a facet. A rigid wall cannot meet itself inside the cell.
    - One-cell pieces are evaluated at the end poses (bit for bit the end cuts, else untrusted) and inside each interval between events: `cut_cell_2d_pieces`, `cut_cell_3d_pieces`.
    - Pieces sharing a feature (a wall segment or facet id, a fluid corner) are one track. There is no "volatile label" rule: a feature on both sides of one event persists through it.
    - Only event-free cells are trivial. "One piece at each end" is wrong (step 0: 628 bench cells, 148 in 2-D and 480 in 3-D).
    - Literature (llm-wiki `immersed-boundary-methods`, #104): Bayyuk 1993 merges to keep a cell's topology over a step (thin bodies and steps over one cell excluded), Meinke 2013 §3.2 split cells, Murman 2003 §4 / Badia 2023 §2 space-time cells; none recorded as linking pieces through a step. Checked against them in #104: `docs/guide/advanced/moving-bodies.md`, `cut-solver.md`, tutorial 11.
  - **Integration:** multi-track cells integrate per event interval on every path, Gauss included (global nodes miss short-lived tracks: 0.39 of a cell). Single-track cells copy `sw`.
  - **Explicit tolerances:**
    - `EVENT_GAP` 1e-10, with per-event round-off widths (issue #104: `EVENT_RES` on the translation formulas, `isolate_zeros(width=)` on the screw ones): an interval narrower than `EVENT_GAP` is evaluated if wider than the sum of its two events' widths (`interval_wide`; counted in `narrow_cells`), else its events act as one; its quadrature nodes may be skipped (the `narrow` flag);
    - pieces below 1e-9 of the cell side with no known feature are dropped;
    - `isolate_zeros(uncertain=)`: a double root within round-off (exact tangential contact under rotation) marks the cell untrusted and counts it in `contact_cells`.
  - **Gates:**
    - `fossil_test_piece_map_2d`, `_3d`, `_3d_bodies` (shared checks `fossil_piece_map_3d_checks`): thin plate dead + fresh, narrow gap continuing, and the other hand cases;
    - validation `fossil_val_piece_map_3d` and `fossil_val_moving_conservation`: a passive scalar through the map in 2-D and 3-D, mass to about 1e-14 and a uniform flow moving with the body exact on the exact paths. Its first-order upwind smears the leakage control, so tracks beat the end-pose overlap there by only 0.1-1 %.
  - **Host passes (#105 step 1, every map bit for bit the f855194 one):** the plan (`map_plan_2d` / `map_plan_3d` -> `map_plan_t`: candidates, per cell status, local elements, events, interval layout, evaluation list as flat CSR, one pose table of `pose_interpolate` poses, the Gauss node sets; count, scan, fill); the evaluations (`map_evals_t`: one item per (cell, evaluation), sized by a count pass, evaluated in host chunks of `MAP_CHUNK_BYTES` = 2^19 bytes with one workspace per item: `map_item_2d/3d_count`, `_eval`, then scans and `map_evals_fill`); the tracks (`map_tracks`); the measures by doubling rounds over flat (cell, interval, node) items (`map_pairs_t`, `map_round_*`, sums in node order: the old second pass is gone); the assembly by count / scan / fill (`map_asm_*`). `isolate_zeros` has no internal procedures (`iz_g`, `iz_add`, stack `IZ_STACK`). OpenMP (`!$omp` over cells and items) checked bit for bit apart (4 threads, the 32 quick reference cases).
  - **Device (#105 steps 2-3):** `piece_map_2d_build(device=)` / `piece_map_3d_build(device=)` run the per-item passes as routine-seq kernels (`map_count_*_dev`, `map_eval_*_dev`, `map_fill_dev`, `map_cells_dev`, `map_node_*_dev`, `map_round_test_dev`, `map_reduce_3d_dev`, `map_asm_*_dev`); the placed geometry, the plan's local elements and the evaluations stay resident for one build. The host keeps the plan (candidates, events, poses, Gauss nodes), the 2-D corner flags / 3-D corner states, the scans and the round decisions: device `cos` / `sin` / `atan2` are 1 ulp off the host's in 17 / 14 / 1.8 % of random arguments (RTX 4070, 26.9, nofma; #105 step 0). With `map` the device end poses are cut in full (3-D, not `alpha_only`) or on the host (2-D). `chunk_bytes` (default `MAP_DEV_CHUNK_BYTES` = 2^28) caps one launch set; no output depends on it. CI twins `fossil_test_piece_map_2d_device_twin`, `_3d_device_twin` (4 KB launch sets, ulp negative control; the 3-D retry paths are reached by no test); bench lines `piece map 2-D` / `3-D device vs serial mismatches`. RTX 4070 (shared, indicative): 3-D refine machinery 2.1-6.1x (3.7-10.7x in an earlier unpaired run), 3-D Gauss map 0.40-0.89x in the bench (3.7x sphere 64^3), 1.2-1.9x in the chain (the host plan dominates: 85-91 % of the cone's device step is the map), 2-D no gain (0.64-0.77x in the chain, 0.34-0.89x in the curve bench). Map timing (all three benches): no timer of its own, so the share is the median of differences (call with `map=` - call without) paired within a repetition (call order reversed every other repetition in the two benches), spread = 2 x their median absolute deviation (`pair_spread`), `not resolved` / `unresolved` when the median is not above it (a difference of bests gave negative maps and a -2.85 s order bias before).
  - **Map reference:** `fossil_val_piece_map_ref` (validation; 64 cases: ellipse, NACA 0012 at 256^2/1024^2, sphere, cone at 32^3/64^3; translation, exact rotation, Gauss nq = 4, Gauss + refine; 0.5 and 2 cells per step) writes or checks a 64-bit hash per map field against `src/tests/validation/perf/map-ref-f855194.txt` (`check FILE [quick] [case=N]`; never overwritten: a deliberate change writes a new named file). Comparator `fossil_piece_map_compare` (`piece_map_2d_mismatches` / `_3d_`: every field, reals by bits; CI `fossil_test_piece_map_compare`).
  - **Fresh and dead are keyed on piece counts.** Since issue #103, `cell_totals_3d` gives a cell without facets whose three faces agree an alpha of exactly 0 or 1. Before, its face-sum volume left ±1.3e-16 on solid cells, and swept `cleared` / `covered` undercounted. A solid cell keeps the cell centre as its centroid. Cells whose faces disagree (open surfaces) keep their fraction and stay untrusted.
- `fossil_mass_properties.f90` — mass properties (issue #59, roadmap #58 WP1): pure free routines (compensated moment accumulator over signed tetrahedra, cyclic Jacobi `symmetric3_eigen`); `surface%mass_properties(mass, center, inertia, principal_moments, principal_axes, density, about, volume, status)` sums about the bounding-box centre (a body at 1e6 keeps its tensor), vertices only, rejects non-watertight or V <= 0 (`STATUS_INVALID_INPUT`). `compute_centroid` uses the same vertex-only tetrahedra (it scaled with the facet size before #59). `inertia_to_world` (R I R^T), `principal_frame` (rigid_pose_t); 2-D twin `curve%mass_properties` (Green sums per unit depth: `second_moments` = [int u^2, int v^2, int uv] about the centroid, `polar_moment`, `principal_angle` of the long axis).
- `fossil_deviation.f90` — surface deviation (issue #60, roadmap #58 WP2): one-sided Hausdorff h(A, B) through a batched distance callback (no surface type in the module); `surface%deviation(other, lower, mean, rms, per_vertex, ...)` samples vertices, edge midpoints and sub-triangle centroids (`samples_per_area`, default 4 per mean facet area; centroid quadrature for mean/rms), `surface%hausdorff`; `DEV_STATUS_*`. Certified `upper` by branch and bound (`deviation_certified`): per triangle min(max d_k + rho_T, min_k max_j |c_j - q_k|) (Lipschitz; convexity of |x - foot|, quadratic convergence), or max d_k when the three feet share a facet of B (coincident surfaces settle at once), split in four until upper - lower <= tolerance (default 1e-3 bbox diagonal) or `max_triangles` (DEV_STATUS_NOT_CONVERGED). B enters as an extension of `deviation_evaluator_t` (an internal procedure as actual argument made gfortran put trampolines on an executable stack). `device=.true.`: batches through `query_many_device` (other entered), bit-identical; bench line `deviation device vs serial mismatches`. 2-D twin `curve%deviation`/`hausdorff` (`deviation_sampled_2d`, `deviation_certified_2d`: segments halved; `deviation_evaluator_2d_t`, `curve_evaluator_t`), device line in `fossil_test_curve_bench` too. Adopted (issue #60 step 5) as certified checks in the remesh, decimate, smoothing, hole_fair, hole_refine and alpha_wrap tests (bounds measured, stated with margins). Coplanar surfaces with crossing facets refine along the crossings (cost ~1/tolerance): use the default tolerance there.
- `fossil_feature_edges.f90` — feature edges (issue #61, roadmap #58 WP3): `find_feature_edges` (stable merge sort of half-edges, sharp by `is_sharp_dihedral` with unit normals, boundary, non-manifold; corners; polylines CSR), `surface%feature_edges(edges, kinds, corners, polyline_start, polyline_vertices, angle, ...)`; constants `FEATURE_EDGE_SHARP/_BOUNDARY/_NONMANIFOLD` (not the query `FEATURE_EDGE`), default 30 degrees. `isotropic_remesh` locks with the same predicate and `feature_cos` (bit-identical). `build_features` (explicit: queries are intent(in); dropped by `analyze`/`move`) stores segments + BVH (`fossil_feature_bvh`, plain arrays, `feature_box_touch_oac` shared CPU/device); `box_feature_many(_device)`; `tag_block(feature_band=, feature_tags=)`; enter_device maps `feature_seg/box/link` if built; `query_many(_device)(feature_distances=)` (`feature_point_distance_oac`, exact clamped projection, bit-identical to brute force); bench line `feature device vs serial mismatches` (boxes and distances). 2-D: `curve%feature_corners` (turning angle, open-loop ends; computed per call, no build), curve `tag_block(feature_band, feature_tags, feature_angle)`, `query_many(_device)(feature_distances, feature_angle)`; curve bench line too. `save_features_into_file` (legacy ASCII VTK POLYDATA: POINTS, VERTICES = corners, LINES = edges, cell scalar `kind`), `fossilizer --features file.vtk [--feature_angle a]`.
- `fossil_thickness.f90` — local thickness and gaps (issue #62, roadmap #58 WP4): pure helpers (fixed golden-angle cone, lower median, `ball_radius`); `surface%thickness(thickness, gap, method, cone_angle, rays, vertex_thickness, vertex_gap, invalid, status, separation_angle)` per facet: THICKNESS_RAY (shape diameter, 60-degree cone of 31 rays from the centroid, first hit by `intersect_ray_first(exclude=)`, leaving/entering test, lower median; cone bias in closed form) or THICKNESS_MEDIAL (shrinking ball started at half the axial chord, so never above it; batched `query_many`; Ma's separation-angle stop, default 30 = local feature size, small next to sharp convex creases: the default is an open owner decision). `thin_regions` (union-find over edges), `build_thin` (local size per BVH payload triangle + per-node minimum; `analyze` drops, `move` keeps), `box_thin_many(_device)` through `bvh_box_thin_oac`, `tag_block(thin_factor, thin_band, thin_tags)`, enter_device maps `thin_pay/thin_node` if built; bench line `thin device vs serial mismatches` (3-D and curve). 2-D: `curve%thickness` (ray fan from the segment midpoint, `ray_first_2d`; medial disc), `thin_regions` (runs along loops), `build_thin`, `box_thin_many(_device)` (`curve_box_thin_oac`, exact), curve `tag_block(thin_*)`. `save_thickness_into_file` (legacy VTK), `fossilizer --thickness file.vtk [--thickness_method] [--h --thin_factor]`. `segment_sdf` keeps its own SDF (unchanged).
- `fossil_line_query.f90`, `fossil_cut_cells_2d.f90`, `fossil_cut_cells_3d.f90` — cut geometry (issue #40): exact axis-aligned line crossings (orient2d on the projected triangle, symbolic tie-break; exact facet vertices in `surface%tri_vertices`, device-mapped), 2-D cut cells (exact grid walk, Green, union-find pieces), 3-D faces from `slice_surface(on_plane_low=.true.)` sections and cell volumes by the divergence theorem. Rule everywhere: a point on a grid line/plane belongs to the low side (cells (x_i, x_i+1]). 3-D pieces (`cut%pieces`, -1 = untrusted cell): per-cell union-find over flagged clipped wall polygons + cell-edge fluid intervals (corner states = exact perturbed node parities, piercings = exact edge crossings under the global perturbation, issue #70), joined by shared vertices, piercings, fluid corners, face islands (ray); touching sets count, closed sets only if a void. Per-piece measures (`np`, `first_piece`, `piece_cell`, `piece_volume`, `piece_centroid`, `piece_wet_area/_centroid/_normal`, `piece_face(6, np)`): wall terms per set, face areas by Green over edge intervals + traces oriented by the polygon winding (high face as wound, low reversed), inclusions attached by a +x ray; single-piece cells take the totals. One global symbolic perturbation (issue #69): grid points at P + (eps1, eps2, eps3), the line kernels order ties x, y, z (`psign(bfirst)` on y lines), `line_crossing_before(_2d)` decides a crossing against a point exactly (one orient3d / orient2d, ties by the perturbed sign): node parities agree across axes. Since issue #70 traces are read on the perturbed face (in-plane polygons by their solid side, facets ending on the plane only if the facet across crosses it) and a multi-piece cell is trusted only if `split_closes`; lattice bodies: 19 of 1.57 million cells still flagged (947 before; a facet touching the cell along a cell-edge segment only), none wrong (`fossil_test_cut_lattice_pieces`). Solver support (issue #70, roadmap #68 WP1; design in `docs/guide/advanced/cut-solver.md`): `fossil_cut_graph` — `cut_links_t` (CSR of piece links: piece across, face 1-6 / 1-4, shared area, `unmatched`), `cut_links_2d(cut, links)` (equal fluid intervals across an edge; 2-D cut cells now give `piece_aperture`, `piece_iv_*`), `cut_cells_3d(..., links=)` (regions = the plane section's 2-D pieces, matched to the 3-D pieces through the shared fluid edge intervals; one-piece sides take all; islands by the unexplained face area). Edge fractions (#70): `cut_cells_3d(..., edges=.true.)` gives `cut%edge_x/_y/_z` (one `line_crossings` per grid line, crossings placed between the perturbed nodes by `line_crossing_before`; off by default, the swept geometry calls the cut cells many times); `cut%node_in` always (z lines; the 3-D corner states); 2-D `cut%node_in`; `curve%cut_cells(intervals=.false.)` skips `piece_iv_*` (the 3-D sections do unless `links=`). SRD (#70 step 5): `srd_neighbourhoods(links, piece_cell, piece_volume, piece_centroid, wet_normal, x, y, srd, z=, tau=)` → `srd_t` (CSR members with the piece first, `nov`, `vhat`, `xhat`, `eta_nb`, `eta_self`, `status` SRD_REGULAR/COMPLETE/INCOMPLETE/NO_VOLUME): growth through links inside the 3^d then 5^d tile, ordered by n.o/|o| (then |o|_1, lexicographic o, piece); zero-volume pieces never members. FRD (#70 step 6): `frd_neighbourhoods(links, piece_cell, piece_volume, x, y, frd, z=)` → `frd_t` (members reached from a cut piece by monotone index paths in the 3^d block, breadth first; weights kappa_j / sum kappa). Merge targets (#70 step 7): `merge_targets(links, piece_cell, piece_volume, wet_normal, x, y, mt, z=, tau=, mask=, exclude=)` → `merge_t` (single `target` by n.o, chains through unsuitable pieces to `final_target`, MERGE_OK/CHAIN/CYCLE/ALONE; normal-directed `mix_piece`/`mix_weight` at offsets o_a in {0, sign n_a}, beta weights after Hu et al. 2006, general knowledge). Moments (#70 step 8, `fossil_cut_moments`: `mono_count/_exponents`, `segment_green_moments`, `triangle_moments`, `cell_moments`, MOMENTS_MAX = 4): `cut_cells_3d(..., moments_order=)` → `piece_moments`, `piece_wall_moments`, `piece_face_moments(:, 6, np)`, `face_x/y/z_moments`, int X^a Y^b Z^c in scaled cell coordinates about the cell centre; volumes by divergence along x (walls + the two x faces; one-piece cells' faces from the plane sections); 2-D `curve%cut_cells(moments_order=)` → `piece_moments`, `piece_wall_moments`. Quadrature (#70 step 9): `cut_cells_3d(..., quadrature=.true.)` → `quad_start/_point/_weight` (signed cones from the cell centre over the oriented fluid boundary, 36 points per cone via `tet_rule`; one-piece cells' faces from the plane sections' `piece_seg`, multi-piece from the face arcs; fluid cells 27-point Gauss; `quad_status` 1 if the weights miss the volume), `wall_quad_*` (7-point rule on the wall fans); 2-D `curve%cut_cells(boundaries=, quadrature=)` → `piece_seg*`, `quad_*`, `wall_quad_*`; with `boundaries=` also `piece_seg_src` (issue #106): the curve segment (its first vertex) each wall segment lies on, 0 on cell edges. 3-D `cut_cells_3d(..., walls=)` (issue #106) → `cut_walls_3d_t`: one wall term per (cell, facet) clipped to >= 3 vertices, with its piece (0 in a cell without trusted pieces), vector area (solid to fluid), wet area, centroid, in cell then facet order, the order the cut sums them into `piece_wet_normal` (bit for bit when every polygon of a piece after its first is one fan triangle); the other outputs unchanged bit for bit. On the device since #109 WP2 step 05: `cut_cells_3d_device(..., walls=)` bit for bit (stage 7, after the pieces: count / `scan_dev` / fill per cell through the routine-seq kernels `wall_terms_sets` (a multi-piece cell reads the polygons and sets its fill pass left in its workspace slice: `cell_pieces_layout` offsets) and `wall_terms_clip` (every other cell: the clip of `cell_totals_3d`), shared with the host path, `polygon_terms` written as scalar loops (host bits unchanged)); resident like links (`cut_walls_3d_t%on_device`, `device_key` = host address of cell_start, `enter_device` / `update_host` / `exit_device`, `cut_walls_reset`; `cut_cells_3d` and `cut_cells_3d_device` take walls `intent(inout)`); `alpha_only` resets them; bench line `cut walls 3-D device vs serial mismatches`, CI twin `fossil_test_cut_walls_3d_device_twin`. Coarsening (#70 step 10, `fossil_cut_coarsen`): `cut_coarsen_3d(fine, links, x, y, z, ratio, coarse, coarse_links, fine_to_coarse, status)` / `cut_coarsen_2d` — coarse pieces = components of fine pieces linked inside the coarse cell, everything else summed (moments re-centred), equal to the direct coarse cut cells. Per-cell kernels (WP2 #100 step 1; device-ready, every one `!$acc routine seq`, compiled by nvfortran 26.9): `compute_cut_cells_2d` keeps `walk_all`, the indexes and `node_flags` on the host, then classifies cells; uniform cells take `quiet_cell_2d_kernel`, cut cells `cut_cell_2d_count` (measures + exact counts of pieces, intervals, boundary segments) then, after a prefix sum, `cut_cell_2d_fill` writing the public arrays in place (no growth, no final copy). Each cut cell works in its own slice of one workspace, sized from its crossings X, wall pieces P and moments nm (`cut_cell_2d_sizes`: reals S(12+2nm) + P(13+2nm) + 2nm + 1, integers 7 + 4S + 2X + 14P + 1, S = 4 + X + P), so no cell can overflow (the MAXS = 512 slots and the #95 open-cell fallback are gone). Plain arrays and scalars only: no class, no allocation, no I/O, no internal procedures (keep it so: the device runs the same kernels). `curve%cut_cells_device(x, y, cut, status, intervals, resident)` (WP2 #100 steps 2, 6a): the host walk and indexes are shared with `cut_cells` (`cut_cells_2d_grid`, `_prepare`); the device runs the node flags (`node_parity`, today's exact per-pair tests), the classification and the integer prefix sums (`cut_cells_2d_classify_dev`, `_scan_dev`: chunks of 64 cells, three-level chunk offsets; public, compared with the host `_classify`/`_scan` integer for integer by `fossil_test_cut_cells_device_twin`), the count, the quiet fill and the cut fill; moments order 0, no boundaries or quadrature; bit for bit `cut_cells` (curve bench lines `cut cells 2-D device vs serial mismatches` and `cut cells 2-D chain device vs serial mismatches`; run `fossil_test_curve_bench 256` to cover 1024^2 and print the timing and per-step chain tables). Resident (step 6a): `resident=.true.` leaves the outputs mapped (unstructured `enter data` on locals, then `move_alloc` into the cut: the present table is keyed on the unchanged host addresses); `cut_cells_2d_t%on_device`, `enter_device`, `update_host`, `exit_device` (member arrays); every routine writing a cut (`cut_cells`, `cut_cells_device`, `cut_coarsen_2d`) takes it `intent(inout)` and releases it first (`cut_cells_2d_reset`: `intent(out)` would free mapped arrays); a copy `b = a` of a resident cut is recognised by `device_key` (host address of alpha) and treated as not mapped; a resident cut must get `exit_device` before it goes out of scope (no `final`). The device routine reuses the old cut's host arrays of the same shape (`take`): pageable copies into fresh pages run at 3-4 GB/s against 12-14 GB/s. nvfortran 26.9 workarounds: `norm2_xy` (gfortran's scaled norm2 written out) instead of the `norm2` intrinsic (device and host 1 ulp apart), `!pgi$r opt=1` on `count_core` (miscompiled at -O2 on the device), explicit `present(...)` on every compute construct over `enter data` arrays (an element expression such as `first_piece(i, j) - 1` triggered a bogus implicit lookup, "partially present"), and no `enter data` on an allocatable dummy (it leaks a descriptor mapping): build in locals, then `move_alloc`. RTX 4070, 1024^2, best of 5: CPU 41-68 ms, copied back 27-32 ms (1.4-2.1x), resident 6.4-9.3 ms (6.4-7.2x); per time step (move + tag + cut) 1.8-2.3x / 3.6-4.3x, the resident step bound by the tags (`box_classify_many_device` uploads the cell boxes each call). A ~1 s slow window of pageable copies follows a GPU-idle phase: time after warm-up calls. 3-D per-cell kernels (WP2 #100 step 3a, CPU, every hash unchanged): `cell_candidates_3d` bins the facets by box and keeps a candidate only if `box_triangle_overlap` (conservative, inflated by `box_tolerance` of grid and surface) touches the closed cell box, in list order (needle facets, |n|_1 < 1e-2 e_max^2, are kept untested; a dropped facet's atom was a singleton never counted and set roots are least atoms, so pieces and sums are unchanged; the cone at 128^3 makes 40,616 non-empty calls of 1,048,152); `cell_totals_3d` is the clip + totals as a per-cell gather in increasing facet order (it replaced a facet-to-cell scatter); clip order (issue #107): x, y, z, except a facet nearly parallel to a y or z plane (`flat_axis`: extent > 0, the other normal components <= `FLAT_TAU` = 1e-3 of a nonzero n_a), clipped along that axis first on its own vertices (`clip_cell` / `clip_cell_flags`, used by `cell_totals_3d`, `wall_terms_clip`, `cell_pieces_core`), as the plane section reads it: the x, y, z order rounded later crossings onto the plane and counted the in-plane region in both cells (up to 1/3 of a cell; 1000-ulp slopes 2.7e-4; no existing case hit it, so no baseline moved; `fossil_test_cut_ulp_plane`; `overlap_weights_3d` dropped its `exact_shift` workaround); `cell_piercings` computes the exact edge piercings (routine seq since step 3c, `fill` flag; `line_crossing_before` and `orient3d_xyz` are device routines); `cell_pieces_3d` / `cell_pieces_core` and the former internal procedures are `!$acc routine seq` module procedures on slices of one workspace per cell (`cell_pieces_sizes` / `cell_pieces_layout` from candidates, piercings, moments, quadrature; mean 9-11 KB, max 54 KB at 128^3 without moments or quadrature); `cut_pieces_3d` counts, scans and fills in place (no growth); cells whose candidates all miss and whose corners agree take the quiet path; `norm2_3` replaces `norm2`. `links_3d` keeps its internal procedures and growth (step 5). Sphere 48^3: `cut_cells_3d` 0.162 s (v3.0.0) -> 0.082 s, with links 0.205 -> 0.115 s; cone 128^3 5.1-5.7 -> 2.3 s. Plane sweep (step 3c): `slice_planes_prepare` / `slice_surface_plane` / `slice_surface_planes` (`fossil_curve_slice`) list, once per axis, the facets whose closed span contains each level, in increasing facet index (a superset of what `slice_surface` keeps, in its order: the same curves bit for bit; NaN facets listed everywhere); both section paths use it (bunny 128^3 slicing 2.9-3.2x faster). Device 3-D cut cells (step 3c): `cut_cells_3d_device(surface, x, y, z, cut, status, edges, resident)` (free routine: the module uses the surface type, a TBP would be circular) = `cut_cells_3d` without links bit for bit (moments order 0, no quadrature, no link intervals; the surface needs its SAH BVH, mapped by `enter_device` or copied in). Host: grid check, `cell_candidates_3d`, `open_cells_3d`, slicing. Device: the sections (`sections_axis_device`), the node states (`edge_fractions_dev` / `lines_dev`: one thread per grid line, `bvh_line_count_oac` then `scan_dev` offsets then `line_fractions_oac` with `bvh_line_crossings_oac`'s sorted crossings, the host's order, placed by the device `line_crossing_before`; entry `cut_cells_3d_nodes_device`), `cell_totals_3d`, classification and workspace sizes with `scan3_dev`, `cell_piercings` and the count pass of `cell_pieces_3d`, `pieces_scan_dev`, then `single_piece_3d` or the fill pass. Resident as 2-D: `cut_cells_3d_t%on_device`, `device_key`, `enter_device` / `update_host` / `exit_device`; `cut_cells_3d`, `cut_cells_3d_sections`, `cut_coarsen_3d` take the cut `intent(inout)` and release it (`cut_cells_3d_reset`). Bench lines `orient3d`, `cut cells 3-D nodes`, `cut cells 3-D` device vs serial mismatches (`fossil_orient3d_cases` holds the predicate cases: avoid subnormal coordinates, #101 item 9b); CI twin `fossil_test_cut_cells_3d_device_twin`. ptxas: `cell_pieces_core` 4008 B (716 B spills), `line_fractions_oac` 776 B; the runtime sizes the stack (no setting needed). RTX 4070, 128^3, best of 5: CPU 0.54-0.68 s, copied back 2.1-3.0x, resident 2.7-3.8x; a resident sphere call: kernels ~41 ms, the faces copied back by the sections and mapped again ~87 ms, host candidates/open cells/slicing ~100 ms (what remains to move). Links (WP2 #100 step 5): `links_3d` / `cut_links_2d` build the links without growth or a global sort: per face (per edge) a capacity, a prefix sum and a fill into slots (low piece, high piece, area) through routine-seq kernels (`links_face_capacity`, `links_face_regions` = today's matching + `resolve` in region order, `links_side_piece`, `links_resolve`), then a per-piece CSR pass (`links_face_id`, `links_face_bases`, `links_piece_count_3d`, `links_piece_fill_3d`; 2-D: counting sort by edge + stable insertion sort by start, pairs, per-piece CSR); equal keys come from one face in region order, the order `build_csr`'s stable sort summed them, so `cut_links_t` is bit for bit the same (`build_csr`/`merge_sort_links` stay for `cut_coarsen`). links_3d 3.2-4.7x at 128^3, now N^3. Device: `cut_link_inputs_3d_t` / `cut_cells_3d_link_inputs_device` (piece intervals `piv_*` with `cell_iv_start`, the sections kept resident with intervals; host twin `cut_cells_3d_link_inputs`), `cut_cells_3d_device(..., links=, resident=)` (`cut_links_3d_device`, `links_3d_dev`: areas summed sequentially per piece), `cut_links_2d_device(cut, links, resident)` (owners, gather per edge, the host kernels); `cut_links_t%on_device`, `enter_device` / `update_host` / `exit_device`, `cut_links_reset` (every writer takes links intent(inout)). 128^3 cut + links resident 6.6-9.4x, 2-D 1024^2 10.5-12x. nvfortran 26.9: #101 item 6 hit `cell_pieces_core` (`niv` bumped via intent(inout) by `walk_edge` / `faces_islands`, read stale in the caller's loops at -O2): the interval count/fill lives in `piece_link_ivs`; the `do a = 1, niv` touch loop of `cell_pieces_core` is a latent instance. Step 4 (WP2 #100): `cut_cells_3d_device` keeps the faces on the device between the sections and the cells (`sections_axis_device(keep, full)`, `scatter_fractions`; a resident 128^3 sphere call 0.21 -> 0.15 s) and has `alpha_only=.true.` (alpha only: it depends on the candidates, the cell edges, the +x/+y/+z face fractions and the six closed flags; the rest is skipped; `is_alpha_only`); `cut_cells_device_many` is a generic over one grid (`_grid`) and a grid per block x(:, b), y(:, b) of equal size (`_grids`, `compute_cut_cells_2d_blocks_grids_device`, `cut_cells_2d_blocks_core`; the 2-D AMR batching of equal-size blocks). Multi-block (WP2 #100 step 3b): `cut_cells_device_many(curves, x, y, cut, status, intervals, pieces, resident)` (`fossil_curve_object`; free routine, the curve arrays are private) runs many curves on one grid in one launch set (`compute_cut_cells_2d_blocks_device`: per-block host walk and indexes concatenated, each kernel reads its block's part through the part's start element; cells g = (b-1) nx ny + c; `cut_cells_2d_classify_blocks_dev`; the scan is `_scan_dev` on the stacked grid nx x (ny nb)); `cut_cells_2d_blocks_t` (cell arrays [.., nb], pieces concatenated with `block_piece`/`block_iv` offsets; block b = `curves(b)%cut_cells` bit for bit, `fossil_test_cut_cells_blocks_twin`). `cut_cells_3d(..., device_sections=.true.)` (only without links, moments, quadrature) and `cut_cells_3d_sections(..., device=)` (faces only, for benches) cut the 3(N+1) plane sections that way (cells only, resident) and scatter them into the face arrays on the device (`scatter_sections`); bench line `cut cells 3-D sections device vs serial mismatches`. RTX 4070, 128^3: host 0.31-0.85 s, device 0.11-0.70 s (1.2-2.9x); the host `slice_surface` (every facet per plane) is ~60 % of the device path, the face copy back (fresh pages) ~10 %. Root accumulators, moments included, move to the winner on every union (issue #95).
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
- `fossil_overset_patch.f90` — the geometry of a conservative overset interface (issue #106, roadmap #92 WP4; design in
  `docs/guide/advanced/overset.md`, section "Conservative interface"; literature: llm-wiki
  `methods/conservative-overset-interface.md`). Rectilinear major / donor grid only (edges x, y, z), curvilinear
  minor / receptor blocks; CPU only; every geometric decision is the exact cut cells'. **A, the patched interface**
  (Wang 1995, the only globally conservative family): `overset_patch_2d(blocks, minor, x, y, patch, boundary_name=,
  faces=, status=)` / `overset_patch_3d(..., z, ...)`: the minor faces of the TOPO_BOUNDARY records named `overlap`
  (or `faces=`), oriented with the minor region solid (the inverse of `wall_cutter`), chained into loops by exact end
  coincidence (2-D) / closed by exact directed-edge pairing (3-D faces split along the shorter diagonal, c1 c3 on a
  tie; collapsed triangles dropped), checked (exact self-intersections: overlapping minor regions refused, since the
  cut counts regions by parity; nesting by orientation / signed volume), then cut: `curve%cut_cells(boundaries=)` +
  `cut_links_2d` in 2-D, `cut_cells_3d(links=, walls=)` in 3-D, bit for bit. Smallest faces from the cut's own labels
  (`piece_seg_src`, `walls=`), no geometric matching: CSR per piece (`piece_sf_start`) and per minor face
  (`face_sf_start`, `face_sf`), `sf_piece`, `sf_face`, `sf_minor`, `sf_area` (out of the major piece, into the minor
  cell), `sf_centroid`; `face_vector` (3-D 1/2 (c3 - c1) x (c4 - c2)), `face_covered` = - sum sf_area. Identities:
  per piece sum_f W_f S_f n_f + sum sf_area = 0 and sum sf_area = - piece_wet_normal (bit for bit with one wall chain
  / one fan per extra polygon); face_covered = face_vector inside the major box. Degenerate pieces (a minor boundary
  on a grid line / plane, the low-side rule puts it into the solid-side cells: zero area) flagged, not removed:
  `degenerate`, `merge_into` (the linked non-degenerate piece of largest shared measure; 0 for point-like pieces);
  folding is the solver's. 3-D untrusted cut cells (pieces -1) keep their smallest faces per cell (`uc_cell`,
  `uc_sf_start`, sf_piece 0) with OVERSET_PATCH_UNTRUSTED (0 of 25 688 at 96^3). Statuses OVERSET_PATCH_OK /
  BAD_INPUT / NO_BOUNDARY / OPEN / OVERLAP / NESTING / UNTRUSTED. **B, intersection weights** (Picard et al. 2021:
  states, not fluxes, not globally conservative): `overlap_weights_2d(block, cells, x, y, w, status=)` /
  `overlap_weights_3d`: each receptor (2-D the polygon of its vertices; 3-D 12 triangles, shorter diagonals, ties by
  the least global vertex index so neighbours tile) cut as a hole in a solid box over the donor sub-grid around its
  bounding box, in coordinates shifted to the sub-grid's low corner (unshifted Green sums lose 1e-12 at |x|/h = 100);
  CSR `start`, `donor`, `volume`, `weight`, `centroid`; `receptor_volume` (compensated: weights sum to 1 within
  2 eps), `covered` (1 exactly inside the closed donor box); zero-area donors listed (a side through a grid node);
  flags OVERLAP_WEIGHTS_OK / BAD_INPUT / DEGENERATE / OUTSIDE / OPEN. Tests `fossil_test_overset_patch_2d/3d`,
  `fossil_test_overlap_weights_2d/3d`. Validation `fossil_val_overset_conservative` (a passive scalar, first-order
  upwind, through a cylinder O-grid at 60^2 / 120^2 and a sphere shell at 40^3): composite mass defect A 1.4e-16 /
  2.8e-16 (3-D 1.2e-16; a square ring on grid lines, 59 degenerate pieces folded, 2.8e-16), B 2.1e-2 / 1.0e-2, I
  (bilinear interpolation) 2.2e-2 / 1.2e-2; on the exact partition (major cells by their fluid fraction outside the
  minor polygon) B 3.1e-2 / 1.0e-2, I 3.0e-2 / 1.2e-2; one transfer of a smooth field, receptor mass vs the donor-side
  mass (B 1.1e-15 / 3.7e-15, I 1.9e-4 / 5.0e-5) and vs the exact integral (B 1.2e-4 / 3.0e-5, I 7.3e-5 / 2.0e-5): B
  conserves each transfer, its composite drift is the coupling (receptors overwritten, not fluxed), as large as
  interpolation's. The step-0 probe `fossil_val_overset_conservative_proto` stays (design variants, 48^3/96^3 timings
  and untrusted counts not covered elsewhere).
- `fossil_implicit_points.f90` — implicit points and exact predicates on them (issue #65, the boolean's exact
  arrangement): `implicit_point_t` of kind `IPK_EXPLICIT` (input vertex), `IPK_LPI` (line p,q meets plane r,s,t),
  `IPK_SSI` (coplanar segments meeting, solved in the projection dropping `drop`), canonical input order (one
  construction = one identity, `ip_same`; `ip_approx` rounds by one formula: equal identities, equal bits);
  `ip_orient3d` (implicit vs explicit plane), `ip_orient2d` (three implicit points in a projection), `ip_compare`
  (coordinate order): homogeneous (X, W) polynomials, a running-error filter in doubles, then exact expansion
  arithmetic (18-bit bit-mask splits, FMA-safe); since roadmap #92 WP5 (d) a coordinate fixed by the construction is
  compared as that double (`exact_coordinate`: an explicit point; LPI/SSI/EEI whose line ends share coordinate k; an
  LPI whose three plane points share k), exact by construction, skipping the exact stage axis-aligned ties reached
  (`fossil_test_implicit_points_axis`: 750 k compares against an integer oracle, 0 mismatches). Golden signs in `src/tests/implicit_points_golden.txt` (exact
  rational reference). Also `IPK_EEI` (two exactly coplanar lines), `ip_lex_compare` (exact lexicographic order),
  `ip_identity_compare`, `ip_rank`. Also `IPK_TPI` (issue #91: three planes through three input points each, Cramer about a1,
  W degree 6, X degree 7; identity = sorted plane triple; rank 5); stage-2 intermediates are compressed (Shewchuk
  Compress: signs only, the older kinds round to the same bits). `fossil_arrangement%add_triple_points`: two cuts of
  one facet made by distinct facets that cross strictly inside it (exact `ip_orient2d`) get their crossing point in all
  three facets: the TPI for two transversal cuts; for a coplanar partner's cut against a transversal one,
  `ip_lpi(partner edge, transversal plane)` (`partner_edge` finds the partner's edge carrying the cut: both cut ends
  exactly collinear with it); for two partners' cuts `ip_eei` of their edges (roadmap #92 WP5 (a)); each is the
  identity the pair of the two facets gives, so it merges by identity. Never fires with clean operands. `fossil_cdt_implicit%cdt_facet`: constrained triangulation of a facet on implicit
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
  `fossilizer --distance other.stl`. CPU only. Performance (issue #94): `body_distance` poses each facet once per call
  (`posed_facets_t`), skips a leaf-pair row whose box is STRICTLY beyond the pruning bound from the union of the other
  leaf's boxes, and skips `pair_cut` on strictly disjoint boxes of two regular triangles (finite, exactly
  non-degenerate, safe exponents: `pair_cut` is not exact on degenerate ones); `body_distance_brute` stays the
  independent reference. Keep the strict comparisons: an equal-distance pair may still win the tie rule. Roadmap #92 WP5
  (d): `leaf_pair` evaluates the pair of least box gap first (ties: loop order), then the rest in loop order (order-free
  result: the running best never drops below the final minimum, the strict update keeps the winner), leaf id lists read
  without copying (`facet_id_count` / `facet_id_at`, `row_skipped`); `try_pair_posed` skips a pair after the best in
  (i, j) order once d2 = 0 (contact pairs 35.5 -> 3.4 per call). Optional `pairs=` (pair evaluations) on
  `body_distance`, `body_distance_brute`, `surface%distance_to`. With the `ip_compare` shortcut: two 12-facet cubes
  35.9 -> 14.8 us per call (27.3 evaluations apart, of 144; box gaps allow no fewer than 26.7), dragon-fine vs posed
  copy 3.3 -> 1.26 ms, CCD 1.13x / 1.28x (sliding), outputs bit for bit; what remains is the 15 sub-distances per
  pair. `fossil_test_body_distance_small` (multi-leaf grid part catches a broken skip; the cube parts do not). Degenerate
  facets (issue #96): `pair_cut` requires non-degenerate triangles (its header says so; a degenerate one has no plane, so
  its answer depended on the argument order); `tri_tri_distance_core` decides contact for an exactly degenerate pair
  (`facet_drop == 0`) with `degenerate_pair_meet` (public in `fossil_tritri_exact` since issue #99; segment/point vs
  triangle by orient3d sides and edges, orient2d in the plane; segment vs segment by exact coplanarity, lexicographic
  intervals when collinear), then the 15 sub-distances.
  Every `pair_cut` caller filters degenerate pairs (issue #98, plan #99): the boolean arrangement (`collect_pair`) by
  `degenerate_pair_meet` (no contact: the zero-area facet passes through uncut; contact: `BOOL_STATUS_INVALID_INPUT`,
  surface unchanged; in `resolve_self_intersections` a sliver closing a T-junction is dropped once the facet across
  it is split at its middle vertex), `find_self_intersections` by `degenerate_pair_contact` (exact SEPARATED /
  CROSSING / TOUCHING, `beyond_shared` for adjacent pairs; fans with a degenerate facet skipped).
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
- `fossil_shape_operator.f90` — the shape operator at a foot (#109 WP1 step 03; plain-real `!$acc routine seq` routines shared by host and device, no vector_R8P): `surface%build_curvature(status)` stores per vertex k1 >= k2 (`curv_k(2, nv)`), unit directions (`curv_d(3, 2, nv)`, tangent to the pseudo-normal, (d1, d2, n) right-handed) and the corner map `curv_vid(3, nf)` (public, like `thin_pay`), from `build_jet_curvature(principal=, directions=)` (the jet's H, K unchanged bit for bit; S_ij = -c_i^T II c_j / hs in the orthonormal basis u1 = r_x/|r_x|, u2 = N x u1 of the fitted graph: k1 + k2 = 2H, k1 k2 = K to round-off; eigenvectors projected on the pseudo-normal's plane; unfitted vertices H +- sqrt(max(H^2 - K, 0)) and the jet frame); `analyze` drops it, `move` rotates `curv_d`; `enter_device` maps it if built, `exit_device` releases it. `foot_principal_oac` interpolates the tensors k1 d1 d1^T + k2 d2 d2^T barycentrically (the `curvatures=` rule, near-collinear facets along the longest edge), restricts to the foot normal's tangent frame (`tangent_frame_oac`: e_x x n, or e_y x n when |n_x| >= 0.9) and `eig_sym2_oac` (closed form, sqrt only: device and host bitwise; the eigenvector from the column without cancellation; umbilic -> t1). `query_many` and `query_many_device` (second kernel on the first's feet, normals, facets) take `principal_curvatures(2, n)`, `principal_directions(3, 2, n)` (not built: STATUS_INVALID_INPUT + stderr; pose: directions rotated). `normal_curvature(k, d, t)` (Euler, re-exported). 2-D: `curve%build_curvature` (private `curv_k(nv)`, `get_vertex_curvature`): `fossil_curvature%vertex_curvature_2d` = Menger circle minus its first-order spacing term kappa_s (h2 - h1)/3 (kappa_s from a 5-point cubic jet; O(h^2) on any spacing, = the circle on uniform spacing; the plain jet was 2-3x worse on ellipses, the plain circle O(h) on alternating spacing); closed loops < 5 vertices and open-loop ends: the (one-sided) circle; `principal_curvatures(1, n)` by `curve_foot_curvature_oac`, `principal_directions(2, 1, n)` = (-n_v, n_u). Accuracy: icosphere R = 2 0.098/0.023/0.0059, cylinder 0.32/0.074/0.018, torus k O(h^2), in-plane angle O(h^2), the plane tilt O(h) (facet normal); marching-cubes vertices (position error O(h^2)) give O(1) curvature noise (spheroid 17 % vs 0.18 % on a scaled icosphere). Tests `fossil_test_principal_curvature`, `_2d`; bench lines `curvature device vs serial mismatches`. Docs: curvature.md "Principal curvatures and directions", fd-closure.md "Wall curvature" (stress-free wall d_n v_t = S v_t).
- `fossil_ghost_points.f90` — ghost points of a finite-difference wall closure (#109 WP1 step 03; free routines, re-exported by `fossil`; page `docs/guide/advanced/fd-closure.md`): `ghost_points_3d(surface, x, y, z, reach, gp, status, device)` / `ghost_points_2d(curve, x, y, reach, gp, status, device)` on the solver's own point lines (centres or nodes, >= 2 per axis); solid = the cut cells' exact perturbed node states (host `cut_cells_3d_nodes`, public since this step, = `cut%node_in` bit for bit; device `cut_cells_3d_nodes_device`; 2-D `node_flags`, `node_rows` + `node_parity`, now public); ghost = solid with a fluid point within reach(a) steps along an axis-a line inside the block; `ghost_points_t` (`n`, `ijk`, `points`, `mask` with bits GHOST_XM/XP/YM/YP/ZM/ZP = 2(a-1) for -a, 2a-1 for +a, `steps(2 dims, n)` = steps to the first fluid point, 0 if none), lexicographic order (i fastest). Device: count (`ghost_mask_oac`), `links_scan_dev`, fill (`ghost_fill_oac`), bit for bit the host; the 3-D device node states are copied back and mapped again (the closure table keeps them resident); the device path needs the SAH BVH. Tests `fossil_test_ghost_points_3d`, `_2d` (hand cases incl. the low-side rule, brute force, device twin).
- `fossil_fd_closure.f90` — the closure table of a finite-difference wall closure (#109 WP1 step 04; free routines, re-exported by `fossil`; page fd-closure.md): `fd_closure_3d(surface, x, y, z, reach, closure, status, lengths, pose, linear_velocity, angular_velocity, centre, ilw, device, resident)` / `fd_closure_2d(curve, x, y, ...)` (2-D angular velocity scalar) fill `ib_closure_t` (one type, `dims` 3 or 2; every array allocated, zero last extent when not asked): the ghosts (as `ghost_points_t`), distance, foot, normal, element, feature (3-D FEATURE_*, 2-D closest vertex), body, patch; probes at foot + lengths(k) n (lengths absolute, + outward, the same for every ghost), `visible`, `probe_cell` (low side: x(i) < p <= x(i+1) -> i, outside 0); `kappa`/`direction` if `build_curvature` ran; `velocity` = V + omega x (foot - c), c default pose%t; `ilw=.true.`: per slot with steps s > 0 the crossing nearest the fluid point between the points s-1 and s (sorted `line_crossings` in [elo - h, ehi + h], `ilw_select_oac` by exact `line_crossing_before`; 2-D row/column buckets of `node_rows`, `ilw_item_2d_oac`), its pseudo-normal and curvature by a query of the crossing (= the cut cells' edge fraction bit for bit on cut-cell nodes). Pose: ghosts and ILW on the moved geometry without moving the surface (3-D since step 06: `posed_tree_3d` places the vertices by `to_world_3d_oac` and refits a copy of the SAH BVH to them level by level (`bvh_node_refit_oac`; min/max exact, so = `move`'s vertices and boxes bit for bit), on the device on the device path; node states `cut_cells_3d_nodes_tree(_present)`, ghosts `ghost_points_nodes_3d`; an octree (host only) keeps a moved copy; 2-D vertices by `to_world`), queries in the body frame and mapped back (scalar replicas `to_body_3d_oac` etc., the same arithmetic). Host path = `ghost_points_*`, `query_many`, `probe_many`, `line_crossings`; device path = one pipeline on present arrays (`cut_cells_3d_nodes_present`, ghost kernels + `links_scan_dev`, `surface%query_points_present` (the refactored kernel of `query_many_device`), curve `query_present`/`probes_present`/`curvature_present`, `bvh_segment_blocked_oac`, `foot_principal_oac`, ILW count/scan/fill + compacted query), bit for bit the host; `resident=.true.` (implies device): built in a local table with `enter data create` per non-empty group (`if(...)`: zero-size arrays are never mapped or referenced in compute constructs), then `move_alloc`; `on_device`, `device_key` (address of `points`), `enter_device`/`update_host`/`exit_device`, `ib_closure_reset`; `ib_closure_mismatches(a, b)` (every field, reals by bits). nvfortran 26.9: no assumed-rank `select rank` (syntax error): use generic specifics. Tests `fossil_test_fd_closure_3d`, `_2d`; bench lines `closure device vs serial mismatches` (+ a closure timing table in the distance bench; `FOSSIL_CLOSURE_TIMING=N1,N2` runs only the closure timing on N^3 / N^2 lines in the distance and curve benches, step 06). Device ILW (step 06): one traversal per slot into a private buffer of `ILW_BUF` = 16 sorted crossings, overflowing slots redone by count/scan/exact fill (`ilw_3d_overflow_dev`; the slab case of the 3-D test reaches it). Curve `query_present` / `probes_present` / `curvature_present` loop in `*_core` routines on dummy arrays: pointer aliases in a `present` clause cost one descriptor upload and download each per launch (28 per curve query; resident 2-D table at 1024^2: NACA 7.6 -> 5.6 ms, disc 5.6 -> 3.9 ms). GPU (step 06, RTX 4070, nvfortran 26.9, indicative): resident 3-D table dragon-fine 128^3 (112 001 ghosts) 0.18 s vs CPU 2.6 s, 64^3 9.4x, 16^3 1.0x (too few threads); kernel-bound (probes + visibility 39 %, the two queries 28 %, ILW 23 %, node states 9 %). GPU-validated 2026-11-16 (26.9, RTX 4070): closure and every other device line 0 on bunny / cube-01 / dragon-fine 32^3, curve bench 512, `fossil_val_device_chain` 0. Accuracy with planar facets (step 05): validation `fossil_val_fd_closure` (circle, icosphere; image-line closure of a smooth field from the table's probes and cells): exact geometry gives orders m + 1 (m = 1..3 probes, Dirichlet and Neumann); facets h_s ~ h cap Dirichlet at order 2 (sag s ~ h_s^2/(8R) at the foot, weight Lambda_m 2.4/4.1/6.0); frozen facets: Dirichlet floor 0.85-0.89 of Lambda_m |d_n u| s, Neumann first order (no floor); page section "What accuracy to expect".
- `fossil_surface_loads.f90` — load samples of the wall integrals (#109 WP2 step 01; free routines, re-exported by `fossil`; page `docs/guide/advanced/surface-loads.md`, conventions L1-L3: F = int (-p n + tau) dA, n outward; facet-plane quadrature with the facet normal; the solver evaluates the flow): `load_samples_3d(surface, x, y, z, samples, status, order, lengths, pose, device, resident)` / `load_samples_2d(curve, x, y, ...)` (lines optional, all or none) fill `load_samples_t`, one row per sample, s = (e-1) q + iq: `points`, `weight` (rule weight x area/length), `normal` (from the vertices: 3-D (v2-v1)x(v3-v1)/|.|, 2-D (dv,-du)/|.|), `element`, `patch`, `body`; probes point + lengths(k) normal, `visible` (3-D `bvh_segment_blocked_oac` from the point, tol 1e-9 box diagonal, SAH BVH required (octree refused with lengths); 2-D `curve%probes_present` on the device, the new TBP `curve%path_visible` on the host, the same kernel), `probe_cell` (`fd_cell_index_oac`). Rules (`load_rule_size`): 3-D centroid / 3-point (2/3,1/6,1/6) / Dunavant 6-point degree 4 (constants solved to 40 digits from the moment equations); 2-D Gauss-Legendre 1-3 points. Degenerate element (zero area/length, missing open-loop segment): weight 0, normal 0, probes not visible. Pose: body-frame samples, probes and visibility on the given tree, mapped by `to_world_3d_oac` / `rotate_3d_oac` (made public in `fossil_fd_closure`): no refit (no decision on the moved geometry); weights and visibility = unposed bit for bit. One set of loops with `if(dev)` (`load_point_3d_oac`, `load_point_2d_oac` the shared routine-seq kernels); resident pattern as `ib_closure_t` (`load_samples_reset`, `load_samples_mismatches`). Tests `fossil_test_load_samples_3d`, `_2d` (weights, exactness to degree 1/2/4 and 1/3/5, visibility blocked by a second body, cells brute force, twin, posed, resident, refusals); bench lines `load samples device vs serial mismatches`. Loads (step 02): `integrate_loads(samples, force, moment, status, pressure, shear, about, per_patch_force, per_patch_moment, per_body_force, per_body_moment, area, device)` (generic: 3-D moment(3), 2-D scalar moment): f = w (-p n + tau), m = (x - about) x f per sample (`load_force_3d_oac` / `_2d_oac`, one `if(dev)` loop; samples, p, tau read on the device where mapped), then Neumaier sums in sample order on the host for both paths (bit for bit, no device reduction); buoyancy F = -g V exact at every order, M = V (xc - c) x (-g) only from order 2 (quadratic integrand: the centroid rule misses it, 5e-3 on a tetrahedron); tests `fossil_test_integrate_loads_3d`, `_2d` (error <= eps x sum of |terms|, measured <= 0.34); bench lines `loads device vs serial mismatches`.
- `fossil_box_predicates.f90`, `fossil_box_query.f90` — AMR tagging (issue #38): box-triangle separating-axis test (conservative: run on the box inflated by `box_tolerance`, 1e-12 of the problem size, never misses a contact), exact 2-D box-segment test; `bvh_box_contact_oac` / `curve_box_contact_oac` (shared CPU/device traversals), `cell_bands` / `cell_components` (tag-band BFS, 26-connected components). The device box path follows the #20 rules and needs no new device arrays.

## Coding conventions

- Fortran 2003+ OOP: all types use `type :: foo_object` with `contains` bound procedures
- All bound procedures use `pass(self)` explicitly
- `destroy` / `initialize` pattern: objects are reset with a fresh local instance (`self = fresh`)
- Numeric types from PENF (`I4P`, `R8P`, etc.) — never use bare `integer` or `real` for library code
- 3D geometry via VecFor `vector_R8P` — not plain arrays

## Known limitations

- **`surface%boolean` / `surface%resolve_self_intersections`**: zero-area facets that touch the other body (or
  another part of the body) are rejected with `BOOL_STATUS_INVALID_INPUT` (issue #98; `sanitize` first). Exact
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
  point (`IPK_TPI`, issue #91): the dragon (open) gives `BOOL_STATUS_NOT_CLOSED`.
  Crossing cuts inside a facet are split at their crossing point: TPI for two
  transversal cuts, LPI (partner's edge against the transversal plane) for a
  partner's edge against a transversal cut, EEI for two partners' edges (roadmap
  #92 WP5 (a), `fossil_test_boolean_coplanar_cross`: four box cases, CDT_FAILED
  before). No `CDT_FAILED` case is known. The configuration needs a
  self-intersecting input (the self-union, or an out-of-contract boolean operand,
  which now gets OK and a non-watertight result, as any self-intersecting operand).
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