<div align="center">

<img src="docs/public/logo.svg" alt="FOSSIL logo: an ammonite meshed in triangles, carved by light" width="128">

# FOSSIL
#### FOrtran Surfaces, Solids, Immersed boundaries and Level sets

[![GitHub tag](https://img.shields.io/github/v/tag/szaghi/FOSSIL)](https://github.com/szaghi/FOSSIL/tags)
[![GitHub issues](https://img.shields.io/github/issues/szaghi/FOSSIL)](https://github.com/szaghi/FOSSIL/issues)
[![CI](https://github.com/szaghi/FOSSIL/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/FOSSIL/actions/workflows/ci.yml)
[![Coverage](https://img.shields.io/codecov/c/github/szaghi/FOSSIL.svg)](https://app.codecov.io/gh/szaghi/FOSSIL)
[![License](https://img.shields.io/badge/license-GPLv3%20%7C%20BSD%20%7C%20MIT-blue.svg)](#copyrights)

> The geometry of complex bodies for simulation codes, in pure Fortran: triangulated surfaces and 2-D curves read and
> repaired, signed distances and level sets, immersed-boundary queries, AMR tagging, exact cut cells and moving bodies,
> on the CPU or bit-exact on the GPU.

</div>

| 📂 **STL in, STL out**<br>ASCII or binary, format detected, bad input rejected; box, volume, connectivity, shells and bodies. [I/O](https://szaghi.github.io/FOSSIL/guide/usage) | 🩹 **Repair dirty meshes**<br>Slivers, duplicates, gaps, flipped normals; holes patched and faired. [Repair](https://szaghi.github.io/FOSSIL/manual/tutorial/02-repair) | 📏 **Signed distance**<br>SAH BVH, bit-exact against brute force; four sign algorithms, the fast winding number for holed meshes. [Signs](https://szaghi.github.io/FOSSIL/guide/advanced/sign-benchmark) | 🌊 **Level sets**<br>Signed-distance fields on stretched grids: exact band, flood-fill signs, fast sweeping. [Fields](https://szaghi.github.io/FOSSIL/guide/advanced/distance-field) |
|:---:|:---:|:---:|:---:|
| 〰️ **2-D curves**<br>Polylines for quadtree solvers with the same queries; sections of any STL. [2-D](https://szaghi.github.io/FOSSIL/guide/advanced/2d-geometry) | 🎯 **Immersed-boundary queries**<br>Closest point, feature, normal and body per point; image points, probes, wall distance. [IB queries](https://szaghi.github.io/FOSSIL/guide/advanced/ib-queries) | 🧭 **AMR tagging**<br>Blocks and cells cut, near, inside or outside from the tree; signed tag bands. [Tagging](https://szaghi.github.io/FOSSIL/guide/advanced/amr-tagging) | 🧱 **Exact cut cells**<br>Exact crossings; 2-D and 3-D fractions, apertures and fluid pieces with their own volumes. [Cut cells](https://szaghi.github.io/FOSSIL/guide/advanced/cut-cells) |
| 🔄 **Moving bodies**<br>Body-frame queries, O(N) tree refit, exact swept volumes of a time step. [Motion](https://szaghi.github.io/FOSSIL/guide/advanced/moving-bodies) | 🛠️ **Geometry processing**<br>Booleans, self-intersections, decimation, remeshing, smoothing, alpha wrap, marching cubes, curvature. [Surfaces](https://szaghi.github.io/FOSSIL/guide/features#surfaces) | 🚀 **On the GPU, bit for bit**<br>OpenACC kernels shared with the CPU path, equal to it bit for bit. [GPU](https://szaghi.github.io/FOSSIL/guide/advanced/device-distance) | 🆓 **Free and portable**<br>Pure Fortran 2003+, no C; gfortran, ifx, nvfortran; FoBiS or fpm; GPLv3, BSD or MIT. [Install](https://szaghi.github.io/FOSSIL/guide/installation) |

<!-- gallery:start -->

<p align="center">
<img src="docs/examples/images/airfoil-cut-cells.svg" alt="animation: the exact cut cells of a NACA 0012 section pitching and translating through a 48 by 24 grid, each cell shaded by its solid fraction" width="420">
<img src="docs/examples/images/airfoil-level-sets.svg" alt="animation: iso-distance contours from 0.05 to 0.6 around the same pitching section" width="420">
<br><sub>A NACA 0012 section pitching and translating through a Cartesian grid, every frame computed by FOSSIL: its exact cut cells and the level sets of its signed distance, with the pose passed to the query instead of rebuilding the tree.</sub></p>

<table>
<tr>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/01-first-surface"><img src="docs/examples/images/tut01-dragon.svg" alt="STL in, STL out: the dragon loaded, 6588 facets, and the box analyze computes" width="260"></a><br><b>STL in, STL out</b><br><sub>the dragon loaded, 6588 facets, and the box analyze computes</sub></td>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/02-repair"><img src="docs/examples/images/anim-holes.svg" alt="Repair dirty meshes: holes closed one by one by fill_holes, refined and faired" width="260"></a><br><b>Repair dirty meshes</b><br><sub>holes closed one by one by fill_holes, refined and faired</sub></td>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/04-distance"><img src="docs/examples/images/tut04-probes.svg" alt="Signed distance: probes and their closest points on the dragon" width="260"></a><br><b>Signed distance</b><br><sub>probes and their closest points on the dragon</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/06-distance-field"><img src="docs/examples/images/tut06-field.svg" alt="Level sets: the signed-distance field on a stretched grid" width="260"></a><br><b>Level sets</b><br><sub>the signed-distance field on a stretched grid</sub></td>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/07-curves"><img src="docs/examples/images/tut07-field.svg" alt="2-D curves: a section of the dragon and its 2-D signed distance" width="260"></a><br><b>2-D curves</b><br><sub>a section of the dragon and its 2-D signed distance</sub></td>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/08-ib-queries"><img src="docs/examples/images/tut08-features.svg" alt="Immersed-boundary queries: closest points on a face, an edge or a vertex" width="260"></a><br><b>Immersed-boundary queries</b><br><sub>closest points on a face, an edge or a vertex</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/09-amr-tagging"><img src="docs/examples/images/tut09-tags.svg" alt="AMR tagging: cut, inside and outside cells of a block, in signed bands" width="260"></a><br><b>AMR tagging</b><br><sub>cut, inside and outside cells of a block, in signed bands</sub></td>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/10-cut-cells"><img src="docs/examples/images/tut10-cells2d.svg" alt="Exact cut cells: the fluid fraction of the cut cells of a section" width="260"></a><br><b>Exact cut cells</b><br><sub>the fluid fraction of the cut cells of a section</sub></td>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/11-moving-bodies"><img src="docs/examples/images/tut11-swept.svg" alt="Moving bodies: the area a moving section sweeps in a time step" width="260"></a><br><b>Moving bodies</b><br><sub>the area a moving section sweeps in a time step</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/cookbook#unite-two-bodies"><img src="docs/examples/images/anim-boolean.svg" alt="Geometry processing: two cubes, then their union, intersection and differences" width="260"></a><br><b>Geometry processing</b><br><sub>two cubes, then their union, intersection and differences</sub></td>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/12-gpu"><img src="docs/examples/images/tut12-slice.svg" alt="On the GPU, bit for bit: the device code path on a slice: 0 differences from the host" width="260"></a><br><b>On the GPU, bit for bit</b><br><sub>the device code path on a slice: 0 differences from the host</sub></td>
<td align="center" width="33%"><a href="https://szaghi.github.io/FOSSIL/manual/tutorial/01-first-surface"><img src="docs/examples/images/tour.svg" alt="Free and portable: the tutorial programs, run one after the other" width="260"></a><br><b>Free and portable</b><br><sub>the tutorial programs, run one after the other</sub></td>
</tr>
</table>

<!-- gallery:end -->

>#### [Documentation](https://szaghi.github.io/FOSSIL/)
> The [tutorial](https://szaghi.github.io/FOSSIL/manual/), the [cookbook](https://szaghi.github.io/FOSSIL/manual/cookbook) and the [reference](https://szaghi.github.io/FOSSIL/guide/features) are on the [FOSSIL website](https://szaghi.github.io/FOSSIL/). Every code sample there is a program compiled and run to produce the output shown.

---

## Authors

**[Stefano Zaghi](https://github.com/szaghi)** · stefano.zaghi@cnr.it
> *Chief Yak Shaver, Accidental Research Scientist, and HPC Farmer* — CFD researcher who decided that one more day debugging Fortran build systems was one day too many, opened a Python REPL "just to prototype," and now finds himself maintaining a meshing library, a chimera assembler, an MPI load balancer, and the seven blog tabs he keeps meaning to read.

**[Claude](https://claude.ai)** (Anthropic)
> *Omniscient Code Oracle & Tireless Rubber Duck* — AI pair programmer, responsible for writing the boring parts so humans don't have to.

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

For interactive use without writing Fortran, the companion **`fossilizer`** CLI (`src/app/fossilizer.f90`) wraps the library for command-line STL analysis and manipulation.

---

## Install

### FoBiS

**Standalone** — clone, fetch the dependencies into `src/third_party/`, and build:

```bash
git clone https://github.com/szaghi/FOSSIL && cd FOSSIL
fobis fetch                     # fetch the fobos [dependencies]
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

### fpm

```bash
git clone https://github.com/szaghi/FOSSIL && cd FOSSIL
fpm build --profile release     # fpm resolves the dependencies from fpm.toml
```

**As an fpm dependency:**

```toml
[dependencies]
FOSSIL = { git = "https://github.com/szaghi/FOSSIL.git" }
```

FOSSIL ships no CMake or Make build files.
