# FOSSIL

**FOrtran Stereo Litography parser** — a pure Fortran 2003+ OOP library for reading, writing, and manipulating [STL](https://en.wikipedia.org/wiki/STL_(file_format)) mesh files.

[![CI](https://github.com/szaghi/FOSSIL/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/FOSSIL/actions)
[![Coverage](https://img.shields.io/codecov/c/github/szaghi/FOSSIL.svg)](https://app.codecov.io/gh/szaghi/FOSSIL)
[![GitHub tag](https://img.shields.io/github/tag/szaghi/FOSSIL.svg)](https://github.com/szaghi/FOSSIL/releases)
[![License](https://img.shields.io/badge/license-GPLv3%20%7C%20BSD%20%7C%20MIT-blue.svg)](#copyrights)

---

## Features

- Load ASCII or binary STL files — format auto-detected with `guess_format=.true.`
- Translate, rotate, mirror, resize, clip, and merge STL surfaces
- Compute volume, centroid, bounding box, connectivity, and disconnected edges
- Sanitize and reverse facet normals; repair disconnected edges automatically
- Signed distance and point-in-polyhedron queries (solid angle or ray intersection)
- AABB octree acceleration — up to 7× faster distance queries over brute force
- OOP/TDD designed — three types (`file_stl_object`, `surface_stl_object`, `facet_object`), all functionality as type-bound procedures
- `fossilizer` CLI for interactive STL analysis and manipulation

**[Documentation](https://szaghi.github.io/FOSSIL/)** | **[API Reference](https://szaghi.github.io/FOSSIL/api/)**

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

Load an STL file, print its statistics, and translate it:

```fortran
use fossil
use penf, only: R8P
use vecfor, only: vector_R8P
implicit none
type(file_stl_object)    :: file_stl
type(surface_stl_object) :: surface

call file_stl%load_from_file(facet=surface%facet, file_name='cube.stl', guess_format=.true.)
call surface%analize
print '(A)', surface%statistics()

call surface%translate(x=1.0_R8P, y=2.0_R8P, z=0.5_R8P)
call file_stl%save_into_file(facet=surface%facet, file_name='cube-moved.stl')
```

Sanitize normals and compute volume:

```fortran
call surface%sanitize_normals
call surface%compute_volume
print *, 'volume =', surface%volume
```

---

## Install

```sh
git clone https://github.com/szaghi/FOSSIL --recursive
cd FOSSIL
```

| Tool | Command |
|------|---------|
| FoBiS.py | `FoBiS.py build -mode static-gnu` |
| FoBiS.py (tests) | `FoBiS.py build -mode tests-gnu && ./scripts/run_tests.sh` |
| make | `./scripts/install.sh --build make` |
| cmake | `./scripts/install.sh --build cmake` |
