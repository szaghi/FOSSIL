---
layout: home

hero:
  name: FOSSIL
  text: FOrtran Stereo Litography parser
  tagline: A pure Fortran 2003+ OOP library for reading, writing, and manipulating STL mesh files.
  actions:
    - theme: brand
      text: Guide
      link: /guide/
    - theme: alt
      text: API Reference
      link: /guide/api-reference
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/FOSSIL

features:
  - icon: 📂
    title: STL I/O
    details: Load and save ASCII or binary STL files. Format is auto-detected — no flags required.
  - icon: 🔧
    title: Surface Manipulation
    details: Translate, rotate, mirror, resize, clip, and merge STL surfaces with a clean type-bound API.
  - icon: 🩺
    title: Surface Analysis
    details: Compute volume, centroid, connectivity, bounding box, watertightness, and disconnected edges.
  - icon: ⚡
    title: AABB Acceleration
    details: Octree Axis-Aligned Bounding Box tree speeds up distance and point-in-polyhedron queries by up to 7×.
  - icon: 🧪
    title: OOP / TDD Designed
    details: Three clean types — file_stl_object, surface_stl_object, facet_object — each tested individually.
  - icon: 🆓
    title: Free & Open Source
    details: Multi-licensed — GPLv3 for FOSS projects, BSD 2/3-Clause or MIT for commercial use. Fortran 2003+ standard compliant.
---

## Quick start

Load an STL file, print its statistics, and translate it:

```fortran
use fossil
use penf, only: R8P
use vecfor, only: vector_R8P

type(file_stl_object)    :: file_stl
type(surface_stl_object) :: surface

! Load (ASCII or binary, auto-detected)
call file_stl%load_from_file(facet=surface%facet, file_name='cube.stl', guess_format=.true.)
call surface%analize

! Print statistics
print '(A)', file_stl%statistics()
print '(A)', surface%statistics()

! Translate and save
call surface%translate(x=1.0_R8P, y=2.0_R8P, z=0.5_R8P)
call file_stl%save_into_file(facet=surface%facet, file_name='cube-moved.stl')
```

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

Contributions are welcome — see the [Contributing](/guide/contributing) page.

## Copyrights

FOSSIL is distributed under a multi-licensing system:

| Use case | License |
|----------|---------|
| FOSS projects | [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html) |
| Closed source / commercial | [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause) |
| Closed source / commercial | [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause) |
| Closed source / commercial | [MIT](http://opensource.org/licenses/MIT) |

> Anyone interested in using, developing, or contributing to FOSSIL is welcome — pick the license that best fits your needs.
