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
`device vs serial max |d| err 0.000E+00`; the `rich`, `box`, `posed`, `sign bench` and `line crossings` mismatch lines must be 0 too (the 3-D bench prints the same lines). Validated on the GPU 2026-11-05 (26.9, RTX 4070, 512²).

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
                                # docs/examples/{snippets,output,images} (~1 min); never edit those by hand
cd docs && npm run docs:build   # API pages (formal) + VitePress site
```

The site follows the FLAP/VecFor/StringiFor/FoXy skeleton (issue #43): Start here (`docs/guide/index.md`,
`installation.md`) → Tutorial (`docs/manual/tutorial/NN-*.md`, 12 chapters) → Cookbook (`docs/manual/cookbook.md`) →
Reference (`docs/guide/features.md` is the feature map, then `usage.md`, `advanced/*.md`, `api/*.md`) → Project. One
sidebar in reading order in `docs/.vitepress/config.mts`. Every tutorial/cookbook snippet and output comes from a
program in `docs/examples/src/` (markers `!run`, `!region`, `!cast`, `!anim`; see the script header); explicit output
formats only. `scripts/anim2svg.py` draws the SMIL geometry animations (`!anim`), `scripts/ansi2svg.py` the terminal
casts, `scripts/mesh2svg.py` the figures (`!fig NAME`: the program writes NAME.fig through
`docs/examples/lib/fossil_fig.f90`; issue #47). The landing gallery and the README's are generated by
`scripts/gallery.py` (regions `<!-- gallery:start/end -->`; never edit them by hand; `--check` runs in
`check_doc_snippets.sh`). Brand: `docs/brand/marks.py` generates `docs/public/{logo,logo-dark,favicon}.svg` (the carved ammonite);
palette in `docs/.vitepress/theme/custom.css`. `scripts/check_doc_snippets.sh` still compiles the fortran blocks of
the hand-written reference pages.

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
- `fossil_aabb_object.f90` / `fossil_aabb_node_object.f90` — AABB primitives
- `fossil_block_object.f90` / `fossil_block_aabb_object.f90` — structured mesh block with AABB (`extrapolate_distances` deprecated, issue #32)
- `fossil_distance_field.f90` — `distance_field_object` (re-exported by `fossil`, issue #32): signed-distance field on a rectilinear grid; narrow band via `distance_many(_device)`, flood-fill signs, fast sweeping in hyperplane order (OpenMP bit-exact). It imports the surface/facet modules, never `fossil` (circular). `fossil_ib_generator` is built on it. `initialize` without `z` makes a 2-D field for a curve source; `compute_band` / `compute_signs` / `compute` / `default_band` are generics over `surface=` / `curve=`.
- `fossil_curve_object.f90`, `fossil_segment_bvh.f90`, `fossil_curve_slice.f90` — 2-D geometry (issue #36): `curve_object` stores (u, v) as SoA arrays in the plane normal to `normal_axis` (u, v = the other axes in cyclic order); the solid is on the left of every segment (outer loops CCW, holes CW), signed distance negative inside; open loops allowed (`# open` in the ASCII format). The segment BVH's traversals are `_oac` kernels shared by the CPU and the device. `slice_surface` cuts an STL by an axis-aligned plane (canonical crossings, no tolerances). `fossil_predicates%orient2d` and `orient3d` (`orient3d_xyz`) are exact (filter + exact expansion with a bit-mask split, safe under nvfortran FMA contraction; orient3d: 24 triple products of input coordinates, 17-bit pieces; no `predicates_initialize` needed).
- Issue #42 (WP6): `distance_field_object%gradient_diagnostic` (eikonal q per node), `wall_distance_many(_device)` (BVH started at cutoff^2; curve twin). `query_many(curvatures=)` (foot curvature: (H, K) from `jet_curvature` — cubic fit over the 2-ring, O(h^2) pointwise; the #18 cotangent / angle-defect operators stall near 15 % at valence-5 vertices — by barycentric weights; 2-D turning angle; CPU only), `surface_markers` (remesh to sqrt(3) dx + one calibration pass) / `curve_markers` (exact arc-length resampling). `isotropic_remesh(preserve_features)` re-detects features after each split (else crease midpoints were relaxed off and edges chamfered).
- `fossil_rigid_pose.f90` — moving bodies (issue #41): `rigid_pose_t` / `rigid_pose_2d_t`; optional `pose=` on the point queries (body-frame, CPU and device), `move(pose)` with `aabb_tree%refit` / `segment_bvh%refit_vertices` (O(N); `aabb_node%set_box` writes the box only — `node%initialize` wipes children), `state_changes`; `fossil_swept_cells.f90`: `swept_cells_2d(curve, ...)` / `swept_cells_3d(surface, ...)` (free routines, not TBPs: the module uses curve/surface/cut_cells_3d) — swept volume = cut-cell volume change (exact), face apertures time-averaged (translation exact via breakpoints: 2-D per face, 3-D per grid plane — mesh vertex on the moving plane, mesh edge on a moving grid line, grid node on a triangle; area quadratic per piece, 2-point Gauss; `exact=.false.` forces Gauss. Rotation: Gauss-Legendre `nq`, nq-vs-2nq estimate), `wall` = -sum A_f n_f, `residual` = rigid-field volume balance; `pose_interpolate` (screw motion), `gauss_legendre_01`.
- `fossil_mass_properties.f90` — mass properties (issue #59, roadmap #58 WP1): pure free routines (compensated moment accumulator over signed tetrahedra, cyclic Jacobi `symmetric3_eigen`); `surface%mass_properties(mass, center, inertia, principal_moments, principal_axes, density, about, volume, status)` sums about the bounding-box centre (a body at 1e6 keeps its tensor), vertices only, rejects non-watertight or V <= 0 (`STATUS_INVALID_INPUT`). `compute_centroid` uses the same vertex-only tetrahedra (it scaled with the facet size before #59). `inertia_to_world` (R I R^T), `principal_frame` (rigid_pose_t); 2-D twin `curve%mass_properties` (Green sums per unit depth: `second_moments` = [int u^2, int v^2, int uv] about the centroid, `polar_moment`, `principal_angle` of the long axis).
- `fossil_line_query.f90`, `fossil_cut_cells_2d.f90`, `fossil_cut_cells_3d.f90` — cut geometry (issue #40): exact axis-aligned line crossings (orient2d on the projected triangle, symbolic tie-break; exact facet vertices in `surface%tri_vertices`, device-mapped), 2-D cut cells (exact grid walk, Green, union-find pieces), 3-D faces from `slice_surface(on_plane_low=.true.)` sections and cell volumes by the divergence theorem. Rule everywhere: a point on a grid line/plane belongs to the low side (cells (x_i, x_i+1]). 3-D pieces (`cut%pieces`, -1 = untrusted cell): per-cell union-find over flagged clipped wall polygons + cell-edge fluid intervals (piercings = polygon vertices on two planes), joined by shared vertices, piercings, fluid corners, face islands (ray); touching sets count, closed sets only if a void. Per-piece measures (`np`, `first_piece`, `piece_cell`, `piece_volume`, `piece_centroid`, `piece_wet_area/_centroid/_normal`, `piece_face(6, np)`): wall terms per set, face areas by Green over edge intervals + traces oriented by the polygon winding (high face as wound, low reversed), inclusions attached by a +x ray; single-piece cells take the totals. Cross-axis on-surface agreement (one global symbolic perturbation) still not done.
- `fossil_box_predicates.f90`, `fossil_box_query.f90` — AMR tagging (issue #38): box-triangle separating-axis test (conservative: run on the box inflated by `box_tolerance`, 1e-12 of the problem size, never misses a contact), exact 2-D box-segment test; `bvh_box_contact_oac` / `curve_box_contact_oac` (shared CPU/device traversals), `cell_bands` / `cell_components` (tag-band BFS, 26-connected components). The device box path follows the #20 rules and needs no new device arrays.

## Coding conventions

- Fortran 2003+ OOP: all types use `type :: foo_object` with `contains` bound procedures
- All bound procedures use `pass(self)` explicitly
- `destroy` / `initialize` pattern: objects are reset with a fresh local instance (`self = fresh`)
- Numeric types from PENF (`I4P`, `R8P`, etc.) — never use bare `integer` or `real` for library code
- 3D geometry via VecFor `vector_R8P` — not plain arrays

## Known limitations

- **`surface%boolean` / `surface%resolve_self_intersections`**: coplanar
  overlaps, shared edges and contacts are handled since issue #33 (exact
  orient3d pair classification in `fossil_tritri_exact`, coplanar cuts,
  T-junction closing). Residual: a 1e-10 relative on-plane snap (features
  below it merge by design); validated on box pairs only
  (`fossil_test_boolean_coplanar`, `fossil_test_boolean_stress`); contact-only
  and touching-symdiff results are non-manifold by definition. Intersection
  points must stay canonical constructions (same bits from every facet pair),
  because the vertex weld compares coordinates exactly (`EPS` is 0).
- **`find_self_intersections`**: pairs sharing one vertex are classified exactly (issue #46); pairs sharing an edge
  are skipped, so two sheets crossing along a shared non-manifold edge are not detected, nor coplanar overlaps around
  a shared vertex (unit normals parallel to 1e-12 rad count as coplanar).

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