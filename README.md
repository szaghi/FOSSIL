# FOSSIL

>#### FOrtran Stereo Litography parser
>a pure Fortran 2003+ OOP library for reading, writing, and manipulating [STL](https://en.wikipedia.org/wiki/STL_(file_format)) mesh files.

[![CI](https://github.com/szaghi/FOSSIL/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/FOSSIL/actions)
[![Coverage](https://img.shields.io/codecov/c/github/szaghi/FOSSIL.svg)](https://app.codecov.io/gh/szaghi/FOSSIL)
[![GitHub tag](https://img.shields.io/github/tag/szaghi/FOSSIL.svg)](https://github.com/szaghi/FOSSIL/releases)
[![License](https://img.shields.io/badge/license-GPLv3%20%7C%20BSD%20%7C%20MIT-blue.svg)](#copyrights)

| 📂 **ASCII & binary STL**<br>Auto-detect format with `guess_format=.true.`; load with on-the-fly clipping; refuse NaN/Inf coordinates via a status code | 🔧 **Surface manipulation**<br>Translate, rotate, mirror, resize, clip, and merge STL surfaces | 📐 **Geometry analysis**<br>Volume, centroid, bounding box, connectivity (symmetric edge-adjacency), and watertight / manifold / volume predicates | 🔨 **Mesh repair**<br>Sanitize normals (outward orientation), drop degenerate slivers and literal duplicates, reconnect nearby vertices, detect non-manifold edges |
|:---:|:---:|:---:|:---:|
| 📏 **Signed-distance queries**<br>Bit-exact closest facet via best-first AABB traversal with d² pruning. Sign via Bærentzen–Aanæs pseudo-normal (default), ray intersection, or solid angle | ⚡ **SAH BVH (default)**<br>Binary BVH partitioning triangles with the surface-area heuristic — ~100× faster than the legacy octree on dragon-scale meshes; both kinds remain selectable | 🏗️ **OOP/TDD designed**<br>Two public types (`surface_stl_object`, `facet_object`), all functionality as type-bound procedures, every public operation under test | 🖥️ **fossilizer CLI**<br>Companion command-line app for interactive STL analysis and manipulation |

>#### [Documentation](https://szaghi.github.io/FOSSIL/)
> For full documentation (guide, API reference, examples, etc...) see the [FOSSIL website](https://szaghi.github.io/FOSSIL/).

| ![dragon](docs/pictures/dragon.jpg) | ![cube](docs/pictures/disconnected-cube.png) |
|:---:|:---:|
| *the dragon STL test (`src/tests/dragon.stl`) has 6588 triangular facets. With the default SAH BVH and the pseudo-normal sign algorithm, signed-distance queries on a 16³ grid run in ~0.022 s vs ~0.987 s on the legacy octree — about a 45× speedup at this size, growing past 100× on larger meshes. Bit-exact correctness vs brute force is asserted by the regression suite.* | *automatic repair of disconnected edges.* |

---

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

Contributions are welcome — see the [Contributing](https://szaghi.github.io/FOSSIL/guide/contributing) page.

## Copyrights

This project is distributed under a multi-licensing system:

- **FOSS projects**: [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html)
- **Closed source / commercial**: [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause), [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause), or [MIT](http://opensource.org/licenses/MIT)

> Anyone interested in using, developing, or contributing to this project is welcome — pick the license that best fits your needs.

---

## Quick start

```fortran
use fossil
use penf, only: R8P
use vecfor, only: ex_R8P
implicit none
type(surface_stl_object) :: surface
real(R8P)                :: d

! Load (ASCII or binary, auto-detected) and run the full repair pipeline.
call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize           ! degenerate / duplicate / non-manifold / winding
print '(A)', surface%statistics()

! Signed distance from a point — uses SAH BVH + pseudo-normal sign by default.
d = surface%distance(point=2.0_R8P * ex_R8P, is_signed=.true., is_square_root=.true.)
print '(A,ES12.5)', 'signed distance = ', d

call surface%translate(x=1.0_R8P, y=2.0_R8P, z=0.5_R8P)
call surface%save_into_file(file_name='cube-moved.stl')
```

See [`src/tests/`](src/tests/) for more examples including clipping, distance queries, and validity predicates (`is_watertight`, `is_manifold`, `is_volume`).

---

## Install

### FoBiS

**Standalone** — clone with submodules and build:

```bash
git clone https://github.com/szaghi/FOSSIL --recursive && cd FOSSIL
fobis build --mode static-gnu   # build static library
fobis build --mode tests-gnu && ./scripts/run_tests.sh  # build and run tests
```

**As a project dependency** — declare FOSSIL in your `fobos` and run `fetch`:

```ini
[dependencies]
deps_dir = src/third_party
FOSSIL = https://github.com/szaghi/FOSSIL
```

```bash
fobis fetch              # fetch and build
fobis fetch --update     # re-fetch and rebuild
```

### CMake

```bash
git clone https://github.com/szaghi/FOSSIL --recursive && cd FOSSIL
cmake -B build && cmake --build build && ctest --test-dir build
```

**As a CMake subdirectory:**

```cmake
add_subdirectory(FOSSIL)
target_link_libraries(your_target FOSSIL::FOSSIL)
```
