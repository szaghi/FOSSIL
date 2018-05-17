<a name="top"></a>

# FOSSIL [![GitHub tag](https://img.shields.io/github/tag/szaghi/FOSSIL.svg)]()

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3%20,%20GPLv3-blue.svg)]()
[![License](https://img.shields.io/badge/license-BSD2-red.svg)]()
[![License](https://img.shields.io/badge/license-BSD3-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-red.svg)]()

[![Status](https://img.shields.io/badge/status-stable-brightgreen.svg)]()
[![Build Status](https://travis-ci.org/szaghi/FOSSIL.svg?branch=master)](https://travis-ci.org/szaghi/FOSSIL)
[![Coverage Status](https://img.shields.io/codecov/c/github/szaghi/FOSSIL.svg)](http://codecov.io/github/szaghi/FOSSIL?branch=master)

### FOSSIL, FOrtran Stereo (si) Litography parser

+ FOSSIL is a pure Fortran (KISS) library for IO and manipulation of STL (Stereo Litography) files for modern (2003+) Fortran projects;
+ FOSSIL is Fortran 2003+ standard compliant;
- FOSSIL is OOP designed;
- FOSSIL is TDD designed;
+ FOSSIL is a Free, Open Source Project.

#### Issues
[![GitHub issues](https://img.shields.io/github/issues/szaghi/FOSSIL.svg)]()
[![Ready in backlog](https://badge.waffle.io/szaghi/FOSSIL.png?label=ready&title=Ready)](https://waffle.io/szaghi/FOSSIL)
[![In Progress](https://badge.waffle.io/szaghi/FOSSIL.png?label=in%20progress&title=In%20Progress)](https://waffle.io/szaghi/FOSSIL)
[![Open bugs](https://badge.waffle.io/szaghi/FOSSIL.png?label=bug&title=Open%20Bugs)](https://waffle.io/szaghi/FOSSIL)

#### Compiler Support

[![Compiler](https://img.shields.io/badge/GNU-v5.3.0+-orange.svg)]()
[![Compiler](https://img.shields.io/badge/Intel-v16.x+-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/IBM%20XL-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/g95-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/NAG-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/PGI-not%20tested-yellow.svg)]()

---

[What is FOSSIL?](#what-is-fossil?) | [Main features](#main-features) | [Copyrights](#copyrights) | [Documentation](#documentation) | [A Taste of FOSSIL](#a-taste-of-fossil)

---

## What is FOSSIL?

FOSSIL is a pure Fortran (KISS) library for IO and manipulation of STL (Stereo Litography) files for modern (2003+) Fortran projects.

FOSSIL provides a simple API to IO STL files and also to manipulate the triangulated surface contained into the STL file.

![dragon](pre_docs/dragon.jpg)

> the dragon STL test (src/tests/dragon.stl) is composed by 6588 triangular facets. The signed distance computation on a uniform
> grid of `64^3` is accelerated by a factor of 7x using AABB algorithm with respect the simple brute force.

![disconnected-cube](pre_docs/disconnected-cube.png)

> automatic repair of disconnected edges.

Go to [Top](#top)

## Main features

* [X] User-friendly methods for IO STL files:
    * [x] input:
        * [x] automatic guessing of file format (ASCII or BINARY);
        * [x] load STL file effortless;
    * [x] output:
        * [x] save STL file effortless;
* [x] powerful surface analysis and manipulation:
    * [x] build facets connectivity;
    * [x] sanitize normals:
        * [x] reverse normals:
        * [x] make normals consistent:
    * [x] compute volume;
    * [x] clip surface outside a bounding box;
    * [x] merge STL files;
    * [x] rotate facets;
    * [x] translate facets;
    * [x] mirror facets;
    * [x] resize (scale) facets;
    * [x] compute minimal distance:
        * [x] square distance;
        * [x] square root distance;
        * [x] signed distance:
            * [x] by means of solid angle computation;
            * [x] by means of rays intersection count;
        * [x] AABB (Axis-Aligned Bounding Box) tree acceleration with user defined refinement levels;
    * [x] point-in-polyhedra test:
        * [x] by means of solid angle computation;
        * [x] by means of rays intersection count;
    * [ ] fill holes;
    * [x] check surface watertight:
        * [x] identify disconected edges;
    * [x] connect nearby facets;
* [ ] errors trapping mechanism.

Any feature request is welcome.

Go to [Top](#top)

## Copyrights

FOSSIL is an open source project, it is distributed under a multi-licensing system:

+ for FOSS projects:
  - [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html);
+ for closed source/commercial projects:
  - [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause);
  - [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause);
  - [MIT](http://opensource.org/licenses/MIT).

Anyone is interest to use, to develop or to contribute to FOSSIL is welcome, feel free to select the license that best matches your soul!

More details can be found on [wiki](https://github.com/szaghi/FOSSIL/wiki/Copyrights).

Go to [Top](#top)

## Documentation

Besides this README file the FOSSIL documentation is contained into its own [wiki](https://github.com/szaghi/FOSSIL/wiki). Detailed documentation of the API is contained into the [GitHub Pages](http://szaghi.github.io/FOSSIL/) that can also be created locally by means of [ford tool](https://github.com/cmacmackin/ford).

### A Taste of FOSSIL

FOSSIL is an KISS library:

#### simple load
> effortless load of file (with STL surface analysis)
```fortran
use fossil
type(file_stl_object) :: file_stl ! STL file handler.
call file_stl%initialize(file_name='cube.stl')
call file_stl%load_from_file(guess_format=.true.)
```

#### print STL statistics
> print main informations of STL
```fortran
print '(A)', file_stl%statistics()
```
> upon exection will print something like
```bash
Mesh_1
file name:   src/tests/cube.stl
file format: ascii
X extents: [ 0.000000000000000E+000, +0.100000000000000E+001]
Y extents: [ 0.000000000000000E+000, +0.100000000000000E+001]
Z extents: [ 0.000000000000000E+000, +0.100000000000000E+001]
volume: -0.100000000000000E+001
number of facets: +12
```

#### sanitiza normals
> make normals consistent
```fortran
call file_stl%sanitize_normals
```

#### simple manipulations
> simply manipulate geometry
```fortran
call file_stl%resize(factor=3.4*ex + 2*ey + 0.5*ez) ! ex, ey, ez being axis versors
call file_stl%resize(x=0.5, z=1.2)                  ! scale only x and z axis
call file_stl%mirror(normal=ex)                     ! mirror respect yz-plane
call file_stl%mirror(normal=ex+ey)                  ! mirror respect plane with normal ex+ey
call file_stl%mirror(matrix=matrix)                 ! mirror by a given mirroring matrix
call file_stl%rotate(axis=ex, angle=1.57)           ! rotate around x axis by pi/2
call file_stl%rotate(matrix=matrix)                 ! rotati by a given rotating matrix
call file_stl%translate(delta=3*ex + 2*ey + 0.5*ez) ! translate by a vectorial delta
call file_stl%translate(x=0.5, z=1.2)               ! translate by only x and z delta
```

Go to [Top](#top)
