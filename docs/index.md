---
layout: home

hero:
  name: FOSSIL
  text: FOrtran Surfaces, Solids, Immersed boundaries and Level sets
  tagline: "The geometry of complex bodies for simulation codes, in pure Fortran: triangulated surfaces and 2-D curves read and repaired, signed distances and level sets, immersed-boundary queries, AMR tagging, exact cut cells and moving bodies, on the CPU or bit-exact on the GPU."
  image:
    light: /logo.svg
    dark: /logo-dark.svg
    alt: the FOSSIL logo, an ammonite meshed in triangles and carved by light
  actions:
    - theme: brand
      text: Tutorial
      link: /manual/tutorial/01-first-surface
    - theme: alt
      text: Cookbook
      link: /manual/cookbook
    - theme: alt
      text: Reference
      link: /guide/features
    - theme: alt
      text: API
      link: /api/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/FOSSIL

features:
  - icon: 📂
    title: STL in, STL out
    details: "ASCII or binary, format detected, NaN and Inf rejected with a status; analyze gives the box, volume, centroid, connectivity, shells and bodies."
    link: /guide/usage
    linkText: I/O and analysis
  - icon: 🩹
    title: Repair dirty meshes
    details: "sanitize removes slivers and duplicates, closes gaps and orients the normals outward; fill_holes patches and fairs the holes left."
    link: /manual/tutorial/02-repair
    linkText: Repair
  - icon: 📏
    title: Signed distance
    details: "Closest facet on a SAH BVH, bit-exact against brute force; four sign algorithms, the fast winding number for surfaces with holes."
    link: /guide/advanced/sign-benchmark
    linkText: Distance and signs
  - icon: 🌊
    title: Level sets
    details: "A signed-distance field on a stretched grid: exact narrow band, flood-fill signs, fast sweeping, an eikonal diagnostic."
    link: /guide/advanced/distance-field
    linkText: Distance fields
  - icon: 〰️
    title: 2-D curves
    details: "Closed and open polylines for quadtree solvers, with the same queries as the surfaces; sections of any STL by a grid plane."
    link: /guide/advanced/2d-geometry
    linkText: 2-D geometry
  - icon: 🎯
    title: Immersed-boundary queries
    details: "Per point: closest point, feature, pseudo-normal and body; image points and probes; wall distance with cutoff; wall curvature."
    link: /guide/advanced/ib-queries
    linkText: IB queries
  - icon: 🧭
    title: AMR tagging
    details: "Blocks and cells classified as cut, near, inside or outside from the tree alone; signed tag bands on a block."
    link: /guide/advanced/amr-tagging
    linkText: AMR tagging
  - icon: 🧱
    title: Exact cut cells
    details: "Grid-line crossings decided by exact predicates; 2-D and 3-D fluid fractions, apertures and fluid pieces with their own volumes."
    link: /guide/advanced/cut-cells
    linkText: Cut cells
  - icon: 🔄
    title: Moving bodies
    details: "Body-frame queries without a rebuild, an O(N) tree refit, and the swept geometry of a time step: exact swept volumes."
    link: /guide/advanced/moving-bodies
    linkText: Moving bodies
  - icon: 🛠️
    title: Geometry processing
    details: "Booleans, self-intersections, decimation, remeshing, smoothing, alpha wrapping, marching cubes, curvature, segmentation."
    link: /guide/features#surfaces
    linkText: Surfaces
  - icon: 🚀
    title: On the GPU, bit for bit
    details: "OpenACC kernels shared with the CPU path: distances, rich queries, tagging and posed queries equal the CPU results bit for bit."
    link: /guide/advanced/device-distance
    linkText: GPU queries
  - icon: 🆓
    title: Free and portable
    details: "Pure Fortran 2003+, no C; gfortran, ifx, nvfortran; FoBiS or fpm; GPLv3, BSD or MIT."
    link: /guide/installation
    linkText: Installation
---

## In a nutshell

Load a body, read what FOSSIL computed about it, save it (the first chapter of the
[tutorial](/manual/tutorial/01-first-surface)):

<<< @/examples/snippets/tut01_first_surface.f90

<<< @/examples/output/tut01.ansi{ansi}

The dragon is neither watertight nor manifold as loaded: [chapter 2](/manual/tutorial/02-repair) repairs it. Start
with the [tutorial](/manual/tutorial/01-first-surface), or look up a task in the [cookbook](/manual/cookbook).

## See it run

A NACA 0012 section pitching and translating through a Cartesian grid, every frame computed by FOSSIL
([`anim_airfoil.f90`](https://github.com/szaghi/FOSSIL/blob/master/docs/examples/src/anim_airfoil.f90)): on the
left the exact cut cells, shaded by their solid fraction ([chapter 10](/manual/tutorial/10-cut-cells)); on the right
the level sets of the signed distance ([chapter 6](/manual/tutorial/06-distance-field)), with the body's pose passed to
the query instead of rebuilding its tree ([chapter 11](/manual/tutorial/11-moving-bodies)).

<div class="fossil-casts">
<img src="./examples/images/airfoil-cut-cells.svg" alt="animation: the cut cells of a NACA 0012 section pitching and translating through a 48 by 24 grid, each cell shaded by its solid fraction">
<img src="./examples/images/airfoil-level-sets.svg" alt="animation: iso-distance contours from 0.05 to 0.6 around the same pitching section">
</div>

The tutorial programs, run one after the other on the dragon test body:

<div class="fossil-casts">
<img src="./examples/images/tour.svg" alt="a terminal session running the tutorial programs of FOSSIL on the dragon test body">
</div>

Every program of the tutorial and of the cookbook is compiled and run by
[`scripts/docs_examples.sh`](https://github.com/szaghi/FOSSIL/blob/master/scripts/docs_examples.sh): the outputs above
are captured, not typed.

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

## Copyrights

FOSSIL is distributed under a multi-licensing system: [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html) for FOSS
projects, [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause),
[BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause) or [MIT](http://opensource.org/licenses/MIT) for
closed-source or commercial ones.
