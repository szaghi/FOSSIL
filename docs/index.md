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
    details: Load and save ASCII or binary STL files. Format is auto-detected; bad input (NaN/Inf coordinates) is rejected via a status code.
  - icon: 🔧
    title: Surface Manipulation
    details: Translate, rotate, mirror, resize, clip, and merge STL surfaces with a clean type-bound API.
  - icon: 🩺
    title: Surface Analysis & Repair
    details: Sanitize pipeline detects and fixes degenerate slivers, literal duplicates, disconnected edges, non-manifold edges, and inward-pointing normals. Predicates is_watertight / is_manifold / is_volume summarise the result.
  - icon: 📏
    title: Signed Distance
    details: Bit-exact closest-facet lookup via best-first traversal with d² pruning. Sign via Bærentzen–Aanæs pseudo-normal (default), ray intersection, or solid angle.
  - icon: ⚡
    title: SAH BVH Acceleration
    details: Binary BVH built with the surface-area heuristic — adapts to triangle density. ~100× faster distance queries than the legacy octree on dragon-scale meshes. Octree remains selectable.
  - icon: 🧪
    title: OOP / TDD Designed
    details: Two public types — surface_stl_object and facet_object — every public operation exercised by the regression suite.
  - icon: 🆓
    title: Free & Open Source
    details: Multi-licensed — GPLv3 for FOSS projects, BSD 2/3-Clause or MIT for commercial use. Fortran 2003+ standard compliant.
---

## Quick start

Load an STL file, print its statistics, and translate it:

```fortran
use fossil
use penf, only: R8P
use vecfor, only: ex_R8P

type(surface_stl_object) :: surface
real(R8P)                :: d

! Load (ASCII or binary, auto-detected) and run the full repair pipeline.
call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%sanitize
print '(A)', surface%statistics()

! Signed distance — SAH BVH + pseudo-normal sign by default.
d = surface%distance(point=2.0_R8P * ex_R8P, is_signed=.true., is_square_root=.true.)
print '(A,ES12.5)', 'signed distance = ', d

! Translate and save.
call surface%translate(x=1.0_R8P, y=2.0_R8P, z=0.5_R8P)
call surface%save_into_file(file_name='cube-moved.stl')
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
