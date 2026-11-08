---
layout: home

hero:
  name: FOSSIL
  text: FOrtran Surfaces, Solids, Immersed boundaries and Level sets
  tagline: "Geometry for simulation codes, in pure Fortran: surfaces and 2-D curves read, repaired, processed and measured; signed distances and level sets; immersed-boundary queries and AMR tagging; exact cut cells and their solver support; curvilinear blocks and overset assembly; body-fitted grids; moving bodies and collision; on the CPU or bit-exact on the GPU."
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
    link: /manual/tutorial/01-first-surface
    linkText: First surface
  - icon: 🩹
    title: Repair dirty meshes
    details: "sanitize removes slivers and duplicates, closes gaps and orients the normals outward; fill_holes patches and fairs the holes left."
    link: /manual/tutorial/02-repair
    linkText: Repair
  - icon: 🛠️
    title: Geometry processing
    details: "Exact booleans and self-intersections, decimation, isotropic remeshing, smoothing, alpha wrapping, marching cubes, segmentation."
    link: /manual/tutorial/12-geometry
    linkText: Processing and analysis
  - icon: 🔬
    title: Geometry analysis
    details: "Mass properties, certified deviation and Hausdorff distance, feature edges, local thickness and gaps, curvature, named patches."
    link: /manual/tutorial/12-geometry
    linkText: Processing and analysis
  - icon: 📏
    title: Signed distance
    details: "Closest facet on a SAH BVH, bit-exact against brute force; four sign algorithms, the fast winding number for surfaces with holes."
    link: /manual/tutorial/04-distance
    linkText: Distance
  - icon: 🌊
    title: Level sets
    details: "A signed-distance field on a stretched grid: exact narrow band, flood-fill signs, fast sweeping, an eikonal diagnostic."
    link: /manual/tutorial/06-distance-field
    linkText: Distance fields
  - icon: 〰️
    title: 2-D curves
    details: "Closed and open polylines for quadtree solvers, with the same queries as the surfaces; sections of any STL by a grid plane."
    link: /manual/tutorial/07-curves
    linkText: Curves
  - icon: 🎯
    title: Immersed-boundary queries
    details: "Per point: closest point, feature, pseudo-normal and body; image points and probes; wall distance with cutoff; wall curvature."
    link: /manual/tutorial/08-ib-queries
    linkText: IB queries
  - icon: 🧭
    title: AMR tagging
    details: "Blocks and cells classified as cut, near, inside or outside from the tree alone; signed tag bands on a block."
    link: /manual/tutorial/09-amr-tagging
    linkText: AMR tagging
  - icon: 🧱
    title: Exact cut cells
    details: "Grid-line crossings decided by exact predicates; 2-D and 3-D fluid fractions, apertures and fluid pieces with their own volumes."
    link: /manual/tutorial/10-cut-cells
    linkText: Cut cells
  - icon: 🧩
    title: Cut-cell solver support
    details: "Links between pieces, state- and flux-redistribution neighbourhoods, merge targets, moments, quadrature rules and coarse levels."
    link: /guide/advanced/cut-solver
    linkText: Solver support
  - icon: 🔄
    title: Moving bodies
    details: "Body-frame queries without a rebuild, an O(N) tree refit, and the swept geometry of a time step: exact swept volumes."
    link: /manual/tutorial/11-moving-bodies
    linkText: Moving bodies
  - icon: 💥
    title: Collision detection
    details: "Body-body distance with the closest pair; the time of first contact by conservative advancement, with no tunnelling; the first impact among many bodies."
    link: /manual/tutorial/14-contact
    linkText: Contact
  - icon: 🕸️
    title: Structured blocks
    details: "Curvilinear blocks: topology detection, certified cell validity, point location by tree or stencil walk, interpolation stencils, quality."
    link: /manual/tutorial/13-blocks-overset
    linkText: Blocks and overset
  - icon: 🔗
    title: Overset assembly
    details: "Holes by winding number, hole map or X-rays; donor candidates with their quality; explicit or implicit classification, orphans by cause, connectivity."
    link: /manual/tutorial/13-blocks-overset
    linkText: Blocks and overset
  - icon: 📐
    title: Body-fitted grids
    details: "Projection constrained to patches, feature curves and corners; 2-D hyperbolic extrusion; a wall-spacing and quality report with y1+."
    link: /manual/tutorial/15-body-fitted
    linkText: Body-fitted grids
  - icon: 🚀
    title: On the GPU, bit for bit
    details: "OpenACC kernels shared with the CPU path: distances, rich queries, tagging and posed queries equal the CPU results bit for bit."
    link: /manual/tutorial/16-gpu
    linkText: GPU
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

<!-- gallery:start -->

<div class="fossil-casts">
<img src="./examples/images/airfoil-cut-cells.svg" alt="animation: the exact cut cells of a NACA 0012 section pitching and translating through a 48 by 24 grid, each cell shaded by its solid fraction">
<img src="./examples/images/airfoil-level-sets.svg" alt="animation: iso-distance contours from 0.05 to 0.6 around the same pitching section">
</div>

<p class="fossil-gallery-caption">A NACA 0012 section pitching and translating through a Cartesian grid, every frame computed by FOSSIL: its exact cut cells and the level sets of its signed distance, with the pose passed to the query instead of rebuilding the tree.</p>

<div class="fossil-casts">
<img src="./examples/images/anim-overset.svg" alt="animation: an O-grid around a cylinder crossing a Cartesian background, its holes, receptors and orphans recomputed in every frame">
<img src="./examples/images/anim-hyperbolic.svg" alt="animation: an O-grid marched out of a NACA 0012 section, three layers per frame, up to 24 layers">
</div>

<p class="fossil-gallery-caption">Grids around bodies: an overset assembly recomputed as the cylinder's O-grid moves through the background, with the orphans that appear where it leaves it, and a body-fitted O-grid marched out of the section by hyperbolic extrusion.</p>

<div class="fossil-gallery">
<a class="fossil-card" href="/FOSSIL/manual/tutorial/01-first-surface"><img src="./examples/images/tut01-dragon.svg" alt="STL in, STL out: the dragon loaded, 6588 facets, and the box analyze computes" loading="lazy"><strong>STL in, STL out</strong><span>the dragon loaded, 6588 facets, and the box analyze computes</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/02-repair"><img src="./examples/images/anim-holes.svg" alt="Repair dirty meshes: holes closed one by one by fill_holes, refined and faired" loading="lazy"><strong>Repair dirty meshes</strong><span>holes closed one by one by fill_holes, refined and faired</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/12-geometry"><img src="./examples/images/anim-boolean-csg.svg" alt="Geometry processing: a block drilled nine times, each step checked by an independent boolean" loading="lazy"><strong>Geometry processing</strong><span>a block drilled nine times, each step checked by an independent boolean</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/12-geometry"><img src="./examples/images/anim-deviation.svg" alt="Geometry analysis: a decimation sweep coloured by its certified distance from the original" loading="lazy"><strong>Geometry analysis</strong><span>a decimation sweep coloured by its certified distance from the original</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/04-distance"><img src="./examples/images/tut04-probes.svg" alt="Signed distance: probes and their closest points on the dragon" loading="lazy"><strong>Signed distance</strong><span>probes and their closest points on the dragon</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/06-distance-field"><img src="./examples/images/tut06-field.svg" alt="Level sets: the signed-distance field on a stretched grid" loading="lazy"><strong>Level sets</strong><span>the signed-distance field on a stretched grid</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/07-curves"><img src="./examples/images/tut07-field.svg" alt="2-D curves: a section of the dragon and its 2-D signed distance" loading="lazy"><strong>2-D curves</strong><span>a section of the dragon and its 2-D signed distance</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/08-ib-queries"><img src="./examples/images/tut08-features.svg" alt="Immersed-boundary queries: closest points on a face, an edge or a vertex" loading="lazy"><strong>Immersed-boundary queries</strong><span>closest points on a face, an edge or a vertex</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/09-amr-tagging"><img src="./examples/images/tut09-tags.svg" alt="AMR tagging: cut, inside and outside cells of a block, in signed bands" loading="lazy"><strong>AMR tagging</strong><span>cut, inside and outside cells of a block, in signed bands</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/10-cut-cells"><img src="./examples/images/tut10-cells2d.svg" alt="Exact cut cells: the fluid fraction of the cut cells of a section" loading="lazy"><strong>Exact cut cells</strong><span>the fluid fraction of the cut cells of a section</span></a>
<a class="fossil-card" href="/FOSSIL/guide/advanced/cut-solver"><img src="./examples/images/anim-cut-solver.svg" alt="Cut-cell solver support: the small pieces of a moving ellipse and their merge targets, every frame" loading="lazy"><strong>Cut-cell solver support</strong><span>the small pieces of a moving ellipse and their merge targets, every frame</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/11-moving-bodies"><img src="./examples/images/tut11-swept.svg" alt="Moving bodies: the area a moving section sweeps in a time step" loading="lazy"><strong>Moving bodies</strong><span>the area a moving section sweeps in a time step</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/14-contact"><img src="./examples/images/anim-ccd.svg" alt="Collision detection: conservative advancement to the time of impact: plain bound, 25 steps; directional, 5" loading="lazy"><strong>Collision detection</strong><span>conservative advancement to the time of impact: plain bound, 25 steps; directional, 5</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/13-blocks-overset"><img src="./examples/images/anim-block-walk.svg" alt="Structured blocks: points tracked through an O-grid by stencil walks, across its seam" loading="lazy"><strong>Structured blocks</strong><span>points tracked through an O-grid by stencil walks, across its seam</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/13-blocks-overset"><img src="./examples/images/ck-overset.svg" alt="Overset assembly: an O-grid over a background: the holes, the receptors around them, the O-grid's receptors" loading="lazy"><strong>Overset assembly</strong><span>an O-grid over a background: the holes, the receptors around them, the O-grid's receptors</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/15-body-fitted"><img src="./examples/images/ck-wall-spacing.svg" alt="Body-fitted grids: an O-grid marched out of a NACA 0012 section by hyperbolic extrusion" loading="lazy"><strong>Body-fitted grids</strong><span>an O-grid marched out of a NACA 0012 section by hyperbolic extrusion</span></a>
<a class="fossil-card" href="/FOSSIL/manual/tutorial/16-gpu"><img src="./examples/images/tut16-slice.svg" alt="On the GPU, bit for bit: the device code path on a slice: 0 differences from the host" loading="lazy"><strong>On the GPU, bit for bit</strong><span>the device code path on a slice: 0 differences from the host</span></a>
<a class="fossil-card" href="/FOSSIL/guide/installation"><img src="./examples/images/tour.svg" alt="Free and portable: the tutorial programs, run one after the other" loading="lazy"><strong>Free and portable</strong><span>the tutorial programs, run one after the other</span></a>
</div>

<!-- gallery:end -->

Every figure and every program of the tutorial and of the cookbook is computed and run by
[`scripts/docs_examples.sh`](https://github.com/szaghi/FOSSIL/blob/master/scripts/docs_examples.sh): the outputs above
are captured, not typed.

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

## Copyrights

FOSSIL is distributed under a multi-licensing system: [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html) for FOSS
projects, [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause),
[BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause) or [MIT](http://opensource.org/licenses/MIT) for
closed-source or commercial ones.
