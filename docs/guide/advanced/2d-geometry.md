---
title: 2-D geometry (curves)
---

# 2-D geometry (curves)

## What it does

A 2-D solver on a quadtree needs the same geometry queries as a 3-D one on an
octree: distance, inside/outside, a signed-distance field. FOSSIL answers them
natively in 2-D with `curve_object`, a set of closed polylines (loops) in an
axis-aligned plane.

It does not extrude the curve into an STL. That would make every 2-D query
pay for 3-D triangles. It would also put the query plane along the
extrusion's facet edges, where sign and closest-feature decisions are most
fragile.

What it covers:

- **curves**: loops, metrics, validity, an ASCII format;
- **queries**: signed distance with closest point, segment and vertex;
  generalized winding number; inside test; exact self-intersections. They run
  on a segment BVH;
- **slicing** an STL by an axis-aligned plane, so a 2-D and a 3-D run share
  one body;
- **repair**: joining and closing open chains, merging vertices, orienting
  loops;
- **batches** of distances on the CPU (OpenMP) or the device (OpenACC);
- **fields**: 2-D signed-distance fields with
  [`distance_field_object`](/guide/advanced/distance-field), one quadtree block
  at a time.

![an ellipse with a square hole and a triangle beside it, the curve in gold over a grid of 5000 cells coloured by signed distance, teal inside the solid and ember outside; six query points, white, each joined to its closest point on the curve, one in the hole, two in the solid, three outside](../../examples/images/ref-2d-geometry-queries.svg){.fossil-figure}

Two bodies, three loops, area 4.27099: `distance_many` gives the 5000 cell-centre distances, from −0.575 to 1.330, 1712
of them inside; `compute_distance` the six feet, e.g. the point (0.10, 0.05) in the hole at d = 0.25 from the hole's
top side (`fig_ref_2d_geometry_queries.f90`).

## Conventions

- **Plane.** A curve lives in the plane normal to `normal_axis` (1, 2 or 3;
  default 3). Its in-plane coordinates (u, v) are the other two axes in
  cyclic order: (y, z), (z, x) or (x, y), so u × v points along
  +`normal_axis`. Points are passed as `p(2) = [u, v]`.
- **Orientation.** The solid is on the left of every segment. Outer loops run
  counter-clockwise, holes clockwise, and islands inside holes
  counter-clockwise again. The signed distance is negative inside, as in 3-D.
- **Open loops.** A loop may be open: a polyline, such as a slice of a mesh
  with holes. An open loop has no closing segment and no area, and
  `check_validity` reports it. The distance and the generalized winding
  number still handle it, and the fast winding number sign is meant for it.
- **Validity.** `check_validity` returns the first defect loop by loop (open
  loop, too few vertices, zero-length segment, zero area, wrong orientation for
  the nesting depth). By default it does not look for crossings;
  `check_validity(check_self_intersections=.true.)` also returns
  `CURVE_SELF_INTERSECTING` when two segments cross, touch or overlap
  (`self_intersections` lists them, decided exactly). A curve can pass the
  default check and still cross itself, a section of a non-manifold mesh for
  example; the cut cells and the swept geometry assume it does not.
- **Exactness.** Decisions that must agree from one query to the next use the
  exact `orient2d`: nesting depth, self-intersections, slicing, repair. The
  point queries are plain floating point, as in 3-D, so they run unchanged
  on the device.

![four loops in the plane of normal_axis 3, with u = x and v = y drawn as axes: an outer square counter-clockwise, a square hole clockwise, a triangle island in the hole counter-clockwise, and an open C-shaped loop beside them; chevrons show the direction of every segment and the solid, to their left, is filled teal; the open loop is filled up to the straight chord between its two ends](../../examples/images/ref-2d-geometry-conventions.svg){.fossil-figure}

The nesting depths are 0, 1, 2 and the signed areas 12, −6 and 0.5: the hole's is negative because it runs
clockwise. `check_validity` returns `CURVE_OPEN_LOOP` at loop 4, whose area is 0. The cells are filled by
`is_point_inside(sign_algorithm=SIGN_FAST_WINDING_NUMBER)`, which closes the open loop by the chord between its ends
(`fig_ref_2d_geometry_conventions.f90`).

## Example

[Tutorial 7](/manual/tutorial/07-curves) slices the dragon and computes a field from the section; the cookbook
[cuts a section of a wing](/manual/cookbook#cut-a-2-d-section-of-an-stl) and
[measures the distance to a polygon](/manual/cookbook#compute-the-signed-distance-to-a-polygon). This program builds
a curve from vertex arrays with two defects, checks it, repairs it and queries it:

<<< @/examples/snippets/ref_2d_geometry-build.f90

<<< @/examples/snippets/ref_2d_geometry-check.f90

<<< @/examples/snippets/ref_2d_geometry-query.f90

<<< @/examples/output/ref_2d_geometry.ansi{ansi}

- **Building.** `adopt_loops` takes the vertices of every loop one after the other, `loop_sizes` and, per loop,
  `closed`. The third loop is passed open: its last vertex lies 1e-7 from its first, as a chain that missed closing
  by round-off would.
- **Checking.** `check_validity` returns the first defect, here `CURVE_WRONG_ORIENTATION` (4) at loop 2: the octagon
  lies inside the plate (nesting depth 1), so it is a hole and must run clockwise. The open loop 3 is a second defect,
  not reported: only the first is.
- **Repairing.** `repair(tolerance=1.e-6)` closes the open chain onto its own start (one gap) and reverses the hole.
  The curve is then valid, self-intersections included, with two bodies: the plate with its hole and the triangle.
  The hole's signed area is negative, −0.17678 = −2√2 r² for the octagon of radius 0.25, and the total 2.32322 is
  2 − 0.17678 + 0.5.
- **Querying.** `compute_distance` returns the foot point, the segment and, if the foot is a vertex, the vertex (0 in
  a segment interior). The point in the hole is outside (d > 0). The fast winding number is 1 or 0 to round-off
  here, except at the last point, 0.3 from the plate's edge, where the far-field expansion leaves 0.00096; the
  brute-force sum (`use_brute_force=.true.`) is exact to round-off, and `winding_number_exact` is the integer from
  exact `orient2d` decisions. Both inside tests agree with it.

## How it works

- **Segment BVH** (`fossil_segment_bvh`). The 2-D twin of the SAH BVH: binary,
  16 buckets per axis on box half-perimeters, with median splits below depth
  48 so the fixed traversal stack cannot overflow. Everything is stored in
  flat arrays. Queries are iterative and deterministic: the nearer child
  first, pruning on squared distance, and a candidate wins only if strictly
  closer.
- **Sign.** The sign comes from the closest feature: the side of the segment,
  or the vertex pseudo-normal. In 2-D the angle-weighted pseudo-normal of
  Bærentzen and Aanæs reduces to the bisector of the two segment normals.
  At a vertex the distance is computed from the vertex itself, so the two
  segments sharing it give the same bits.

  Near an open end the two signs differ. The fast winding number closes the gap by the straight chord between the
  ends, where w crosses 0.5; the pseudo-normal takes the side of the nearest feature, an end vertex for the points
  beyond the gap. On the open dragon section they disagree on 39 cells of 0.05, all in a wedge outside the chord:

  ![a close-up of the gap of the open dragon section: the cells inside by the fast winding number in teal, outside in ember, the boundary between them following the curve and the straight chord between its two open ends; a small wedge of gold cells just outside the chord marks where the pseudo-normal sign says inside instead](../../examples/images/ref-2d-geometry-signs.svg){.fossil-figure}
- **Winding number.** w = (1/2π) Σ (angle subtended by each segment): 1
  inside, 0 outside or in a hole, fractional near an open end.
  `winding_number` evaluates it hierarchically, the 2-D analogue of the fast
  winding number (Barill et al. 2018). A node far enough away contributes
  the second-order expansion of (1/2π) ∫ n·(x − p)/|x − p|² ds from its
  dipole and second moments, which are exact for straight segments.

  ![a section of the dragon at z = 9.62 seen from above, with the facets within 0.5 of a point of the section removed so that one loop is open, its two ends marked in ember; the background cells are coloured by the generalized winding number, teal 1 inside the five closed regions, ember 0 outside, and a cream fan of fractional values spreading out of the gap of the open loop](../../examples/images/ref-2d-geometry-winding.svg){.fossil-figure}

  On the open section of the dragon above, 645 of the 11400 cells have 0.05 < w < 0.95: all of them in the fan
  around the gap of 0.957 between the two open ends. Elsewhere w is 0 or 1, as for closed loops.
- **Slicing.**
  - A vertex on the plane counts as above it. A facet is therefore crossed
    along exactly two edges or not at all, and a facet lying in the plane is
    not crossed.
  - Crossing points are built from the edge endpoints in lexicographic
    order, so the two facets sharing an edge produce the same bits. Loops
    close by exact equality, without tolerances.
  - Each segment's direction comes from the facet's vertex order, so
    outward-oriented facets put the solid on the left.
- **Device.** The traversals are free `_oac` kernels (`!$acc routine seq`) on
  the flat arrays, and the CPU queries call them too. `distance_many_device`
  follows the [GPU rules of #20](/guide/advanced/device-distance).

## Results

| Case | Result |
|---|---|
| near-collinear `orient2d` grid (2048 points) | naive formula wrong on 482, `orient2d` on 0 |
| 2000 random points (NACA 0012, cusp, 64-gon) | BVH distance = brute force bit for bit; sign = exact winding number |
| fast winding number, NACA 0012 | inside decisions exact; max error 1.4e-2, 2.0e-3, 2.3e-4 at β = 2, 4, 8 |
| slice of the unit cube at z = 0.4, 0.5, 0.25, 1 | the exact unit square after merging, area exactly 1 |
| slice of an icosphere at z = 0.3, levels 2-4 | area error 6.9e-2, 1.8e-2, 4.7e-3 (O(h²)) |
| `distance_many` / device path vs `compute_distance` | bit for bit (4000 points, both signs) |
| 2-D field, circle, 32², 64², 128² | band bit-exact; far-field error 2.8e-2, 2.0e-2, 1.3e-2 (O(h \|log h\|)) |

## Limitations

- **GPU not validated yet.** The device path is written and runs bit-exact on
  the host, through the same kernel, in the `tests-gnu` build. The GPU check
  (`fobis build --mode tests-nvidia-curve --varset local_nvf`, then
  `fossil_test_curve_bench`) is pending.
- **Coarse fast winding number near thin bodies.** It is about 1e-2 at the
  default β = 2. The inside decision has a margin of 0.5; raise β where the
  value itself matters.
- **No curved segments.** Curves are polylines: arcs and splines must be
  sampled first.
- **Slicing with planes along axes only.** An arbitrary plane needs a
  rotation of the surface first.
- **A coarse join rule in repair.** It joins open chains greedily by nearest
  end within the tolerance. Ambiguous gaps (several ends within the
  tolerance) take the nearest.
- **One tolerance for joins and merges.** `repair` also merges consecutive
  vertices closer than `tolerance`. A tolerance large enough to bridge a
  real gap (a missing facet, as in the figures above) also drops every
  segment shorter than the gap and coarsens the curve. It is meant for
  round-off gaps, as in the example; a real gap needs the open-loop queries
  (fast winding number sign) or a repaired surface.

## API

See [`curve_object`](/guide/api/curve-object) for every method and argument.
