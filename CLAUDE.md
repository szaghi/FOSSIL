# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is FOSSIL?

FOSSIL (FOrtran Stereo Litography parser) is a pure Fortran 2003+ library for reading, writing, and manipulating STL (stereolithography) mesh files. It is OOP-designed and TDD-designed. The companion CLI app `fossilizer` (`src/app/fossilizer.f90`) wraps the library for command-line STL processing.

## Build System

The project uses **FoBiS.py** (Fortran Build System) with a `fobos` configuration file. FoBiS.py must be installed (`pip install FoBiS.py`).

### Common build commands

```bash
# Build static library (GNU compiler)
FoBiS.py build -mode static-gnu

# Build shared library (GNU compiler)
FoBiS.py build -mode shared-gnu

# Build all tests (GNU compiler)
FoBiS.py build -mode tests-gnu

# Build with debug flags (GNU compiler)
FoBiS.py build -mode tests-gnu-debug

# Clean build artifacts
FoBiS.py rule -ex clean
```

Intel compiler variants are available by replacing `-gnu` with `-intel`.

### Running tests

```bash
# Build tests first, then run all
FoBiS.py build -mode tests-gnu
./scripts/run_tests.sh

# Run a single test executable
./exe/fossil_test_clip

# Full coverage analysis (build debug + run tests + gcov)
FoBiS.py rule -ex makecoverage
```

Test programs output `Are all tests passed? T` on success or `Are all tests passed? F` on failure.

### Documentation

```bash
FoBiS.py rule -ex makedoc   # Requires 'formal' and npm/VitePress
FoBiS.py rule -ex deldoc    # Delete doc build artifacts
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

The public API is accessed via `use fossil`, which re-exports three types:

| Type | Module file | Purpose |
|------|-------------|---------|
| `file_stl_object` | `fossil_file_stl_object.f90` | File I/O (open, load, save, close) |
| `surface_stl_object` | `fossil_surface_stl.f90` | Surface geometry: analysis, manipulation, distance queries |
| `facet_object` | `fossil_facet_object.f90` | Individual triangular facet with vertices, normal, connectivity |

### Typical usage flow

1. `file_stl%load_from_file(facet=surface%facet, ...)` — reads STL (ASCII or binary, auto-detected)
2. `surface%analize` — computes bounding box, volume, centroid, connectivity, AABB octree
3. Manipulate: `surface%clip`, `surface%rotate`, `surface%translate`, `surface%mirror`, `surface%resize`, `surface%merge_solids`, `surface%sanitize_normals`, `surface%connect_nearby_vertices`
4. Query: `surface%distance`, `surface%is_point_inside`, `surface%compute_mesh_distance`
5. `file_stl%save_into_file(facet=surface%facet, ...)` — writes result

### AABB acceleration

`fossil_aabb_tree_object.f90` implements an **octree** (8-child refinement). It is embedded in `surface_stl_object%aabb` and built during `analize`. Distance and point-in-polyhedron queries use AABB traversal before brute-force per-facet computation.

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
