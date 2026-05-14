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
  - icon: 🚀
    title: Advanced Geometry Processing
    details: Booleans, winding number, alpha wrap, isotropic remesh, decimation, marching cubes, SDF segmentation, ray queries, cotangent Laplacian, per-vertex curvature, Taubin smoothing — twelve libigl/CGAL-class primitives, every one a single TBP. See the Advanced Features guide.
  - icon: 🆓
    title: Free & Open Source
    details: Multi-licensed — GPLv3 for FOSS projects, BSD 2/3-Clause or MIT for commercial use. Fortran 2003+ standard compliant.
---

## Advanced features

FOSSIL ships a focused subset of the libigl / CGAL geometry-processing
toolkit, drawn directly from the standard references in the field:

- **[Boolean operations](/guide/advanced/booleans)** — union, intersection,
  difference, symmetric difference. Zhou et al. 2016.
- **[Self-intersection detection / resolution](/guide/advanced/self-intersection)** —
  Möller tri-tri broad-and-narrow phase, with one-line resolution via the
  boolean engine.
- **[Mesh decimation](/guide/advanced/decimate)** — quadric edge collapse
  with normal-flip / non-manifold safety. Garland & Heckbert 1997.
- **[Generalized winding number](/guide/advanced/winding-number)** — robust
  inside/outside on dirty STL, hierarchical Barnes-Hut traversal.
  Jacobson 2013 + Barill 2018.
- **[Marching cubes](/guide/advanced/marching-cubes)** — isosurface
  extraction from a regular scalar field, with the SDF→STL roundtrip
  TBP. Lorensen-Cline 1987 + Bourke tables.
- **[Alpha wrapping](/guide/advanced/alpha-wrap)** — guaranteed watertight
  2-manifold surrogate from any triangle soup. Portaneri et al. 2022.
- **[Isotropic remeshing](/guide/advanced/isotropic-remesh)** — uniform-edge
  re-tessellation with optional sharp-feature preservation.
  Botsch & Kobbelt 2004.
- **[SDF segmentation](/guide/advanced/sdf-segmentation)** — per-facet
  semantic labels via Shape Diameter Function + Gaussian-mixture
  clustering. Shapira, Shamir & Cohen-Or 2008.
- **[Ray-mesh intersection queries](/guide/advanced/ray-queries)** — closest
  hit, all hits, any-hit early-exit. Möller-Trumbore + AABB tree.
- **[Cotangent Laplacian + barycentric mass](/guide/advanced/cotangent-laplacian)** —
  sparse SPD operator over the surface; foundation for curvature,
  smoothing, and (when a sparse solver lands) heat-method geodesics.
  Pinkall & Polthier 1993; Meyer et al. 2003.
- **[Per-vertex Gaussian + mean curvature](/guide/advanced/curvature)** —
  angle-defect Gaussian and signed mean curvature from `H n = (1/2) M⁻¹ L V`.
  Meyer et al. 2003.
- **[Mesh smoothing (explicit + Taubin λ\|μ)](/guide/advanced/smoothing)** —
  the right production-grade denoiser for CFD-grade STL. Taubin 1995.

Each algorithm has a dedicated page with a CFD-relevant motivation, a
worked Fortran example, and an honest list of known limitations.

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

**[Stefano Zaghi](https://github.com/szaghi)** · stefano.zaghi@cnr.it
> *Chief Yak Shaver, Accidental Research Scientist, and HPC Farmer* — CFD researcher who decided that one more day debugging Fortran build systems was one day too many, opened a Python REPL "just to prototype," and now finds himself maintaining a meshing library, a chimera assembler, an MPI load balancer, and the seven blog tabs he keeps meaning to read.

**[Claude](https://claude.ai)** (Anthropic)
> *Omniscient Code Oracle & Tireless Rubber Duck* — AI pair programmer, responsible for writing the boring parts so humans don't have to.

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
