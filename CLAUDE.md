# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is FOSSIL?

FOSSIL (FOrtran Stereo Litography parser) is a pure Fortran 2003+ library for reading, writing, and manipulating STL (stereolithography) mesh files. It is OOP-designed and TDD-designed. The companion CLI app `fossilizer` (`src/app/fossilizer.f90`) wraps the library for command-line STL processing.

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

### Documentation

```bash
fobis rule --ex makedoc   # Requires 'formal' and npm/VitePress
fobis rule --ex deldoc    # Delete doc build artifacts
```

### Alternative builds (via install.sh)

```bash
./scripts/install.sh --build make     # GNU Make
./scripts/install.sh --build cmake    # CMake
```

## Source Architecture

```
src/
├── lib/           # Core library modules
├── tests/         # Test programs (one per feature)
├── app/           # fossilizer CLI application
└── third_party/   # Git submodules
```

### Core object hierarchy

The public API is accessed via `use fossil`, which re-exports two types plus the named constants:

| Type | Module file | Purpose |
|------|-------------|---------|
| `surface_stl_object` | `fossil_surface_stl.f90` | Surface geometry: I/O, analysis, sanitization, manipulation, distance queries |
| `facet_object` | `fossil_facet_object.f90` | Individual triangular facet — vertices, normal, connectivity, pseudo-normals |

`fossil.f90` also re-exports:

- `SIGN_PSEUDO_NORMAL` (default), `SIGN_RAY_INTERSECTIONS`, `SIGN_SOLID_ANGLE` — sign algorithms for signed-distance queries.
- `AABB_TREE_SAH_BVH` (default), `AABB_TREE_OCTREE` — AABB tree kinds.
- `AABB_AUTO_REFINEMENT` — sentinel for octree depth auto-tuning.
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
2. `surface%sanitize` — orchestrator: `remove_degenerate_facets` → (`connect_nearby_vertices` if needed) → `analyze` → `remove_duplicate_facets` → `analyze` → `sanitize_normals`. Emits one stderr warning per defect class.
3. `surface%analyze(aabb_refinement_levels=, aabb_tree_kind=)` — invoked automatically by load/sanitize/clip. Rebuilds connectivity and the AABB tree. The optional args are passed through to `aabb%initialize`.
4. Manipulation: `clip`, `rotate`, `translate`, `mirror`, `resize`, `merge_solids`.
5. Queries: `distance`, `is_point_inside`, plus `compute_distance` (richer form with closest facet/edge/vertex outputs).
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
- `SIGN_RAY_INTERSECTIONS` — axis-aligned ray casts; odd intersections = inside.
- `SIGN_SOLID_ANGLE` — sum of projected solid angles ≈ ±4π.

### Third-party submodules (`src/third_party/`)

| Submodule | Key exports used in FOSSIL |
|-----------|---------------------------|
| **PENF** | Portable numeric types: `I4P` (int32), `R8P` (real64), `str()` |
| **VecFor** | 3D vectors: `vector_R8P`, unit vectors `ex_R8P/ey_R8P/ez_R8P`, `rotation_matrix_R8P`, `mirror_matrix_R8P` |
| **FLAP** | CLI argument parsing (`command_line_interface`) used by `fossilizer` |
| **FACE** | ANSI color/style output |
| **BeFoR64** | Base64 encoding |
| **FoXy** | XML/XDMF support |
| **StringiFor** | String utilities |
| **VTKFortran** | VTK file I/O |

After cloning, initialize submodules:
```bash
git submodule update --init
```

### Key internal modules

- `fossil_utils.f90` — constants (`EPS`, `PI`, `FRLEN`) and helpers (`is_inside_bb`)
- `fossil_list_id_object.f90` — dynamic integer ID list used for connectivity tracking
- `fossil_aabb_object.f90` / `fossil_aabb_node_object.f90` — AABB primitives
- `fossil_block_object.f90` / `fossil_block_aabb_object.f90` — structured mesh block with AABB

## Coding conventions

- Fortran 2003+ OOP: all types use `type :: foo_object` with `contains` bound procedures
- All bound procedures use `pass(self)` explicitly
- `destroy` / `initialize` pattern: objects are reset with a fresh local instance (`self = fresh`)
- Numeric types from PENF (`I4P`, `R8P`, etc.) — never use bare `integer` or `real` for library code
- 3D geometry via VecFor `vector_R8P` — not plain arrays

## Known limitations

- **`surface%boolean` / `surface%resolve_self_intersections`** (issue #18 §1.1):
  bit-exact correct on inputs where no facets of A and B are coplanar (the
  generic case). On inputs with axis-aligned coplanar face overlaps (e.g.
  cubes offset along a single axis only), the call returns `BOOL_STATUS_OK`
  but the result volume is wrong. Workaround: perturb one input by
  `eps = 1e-6 * bbox_diagonal` along each axis before the boolean. Full
  details and the upgrade path in the `surface%boolean` docstring's
  "Limitations" section in `src/lib/fossil_surface_stl.f90`.

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