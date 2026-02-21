---
title: Installation
---

# Installation

## Prerequisites

A Fortran 2003+ compliant compiler is required. The following compilers are known to work:

| Compiler | Minimum version |
|----------|----------------|
| GNU gfortran | ≥ 5.3.0 |
| Intel Fortran (ifort / ifx) | ≥ 16.x |

FOSSIL is developed on GNU/Linux. Windows should work but is not officially tested.

## Download

FOSSIL uses **git submodules** for its third-party dependencies. Clone recursively:

```bash
git clone https://github.com/szaghi/FOSSIL --recursive
cd FOSSIL
```

If you already have a non-recursive clone:

```bash
git submodule update --init
```

### Third-Party Dependencies

The submodules live under `src/third_party/`:

| Library | Purpose |
|---------|---------|
| [PENF](https://github.com/szaghi/PENF) | Portable numeric kind parameters (`I4P`, `R8P`, etc.) |
| [VecFor](https://github.com/szaghi/VecFor) | 3D vector type and operations (`vector_R8P`, unit vectors, rotation/mirror matrices) |
| [FLAP](https://github.com/szaghi/FLAP) | Fortran command-line argument parser (used by `fossilizer`) |
| [FACE](https://github.com/szaghi/FACE) | ANSI terminal color/style support |
| [BeFoR64](https://github.com/szaghi/BeFoR64) | Base64 encoding |
| [FoXy](https://github.com/szaghi/FoXy) | XML/XDMF support |
| [StringiFor](https://github.com/szaghi/StringiFor) | String utilities |
| [VTKFortran](https://github.com/szaghi/VTKFortran) | VTK file I/O |

## Build with FoBiS.py

[FoBiS.py](https://github.com/szaghi/FoBiS) is the primary build system. Install it once with pip:

```bash
pip install FoBiS.py
```

### List all build modes

```bash
FoBiS.py build -lmodes
```

Available modes:

| Mode | Description |
|------|-------------|
| `static-gnu` | Static library with gfortran (release) |
| `shared-gnu` | Shared library with gfortran (release) |
| `static-gnu-debug` | Static library with gfortran (debug + checks) |
| `shared-gnu-debug` | Shared library with gfortran (debug) |
| `tests-gnu` | Build all tests with gfortran (release) |
| `tests-gnu-debug` | Build all tests with gfortran (debug) |
| `static-intel` | Static library with ifort (release) |
| `shared-intel` | Shared library with ifort (release) |
| `tests-intel` | Build all tests with ifort (release) |
| `tests-intel-debug` | Build all tests with ifort (debug) |

### Build the library

```bash
# Static library (GNU gfortran)
FoBiS.py build -mode static-gnu

# Shared library (GNU gfortran)
FoBiS.py build -mode shared-gnu

# Static library (Intel Fortran)
FoBiS.py build -mode static-intel
```

The library is placed in `./static/` or `./shared/` respectively. Module files go to `./mod/`.

### Build and run tests

```bash
FoBiS.py build -mode tests-gnu
./scripts/run_tests.sh
```

Compiled test executables are placed in `./exe/`. Each test prints `"Are all tests passed? T"` on success.

To run a single test:

```bash
./exe/fossil_test_clip
```

### Coverage and documentation

```bash
FoBiS.py rule -ex makecoverage   # build + run tests + gcov report
FoBiS.py rule -ex makedoc        # build ford API docs + VitePress site
```

## Alternative builds

The `scripts/install.sh` helper also supports `make` and `cmake` backends:

```bash
./scripts/install.sh --build make
./scripts/install.sh --build cmake
```
