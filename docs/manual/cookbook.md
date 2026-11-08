---
title: Cookbook
---

# Cookbook

Short answers to "how do I ...?". Each recipe is a complete program with its real output; the
[reference](/guide/api/) has the details. Most programs read the test bodies of the repository, `dragon.stl` (binary,
6588 facets, with holes), `cube-01.stl` (the unit cube, binary) and `naca0012-ascii.stl` (a wing section, ASCII),
which are in [`src/tests`](https://github.com/szaghi/FOSSIL/tree/master/src/tests); the others build their geometry in
code. The [tutorial](./index) builds the same calls up step by step.

[[toc]]

## Files

### Read an STL file, ASCII or binary

<<< @/examples/snippets/ck_read_any.f90

<<< @/examples/output/ck_read_any.ansi{ansi}

<div class="fossil-pair">

![the dragon read from the binary dragon.stl: 6588 facets](../examples/images/ck-read-any-binary.svg)

![the NACA 0012 wing section read from the ASCII naca0012-ascii.stl: 188 facets](../examples/images/ck-read-any-ascii.svg)

</div>

`guess_format=.true.` tells the formats apart by the file size: a binary STL of `n` facets is `84 + 50 n` bytes. See
[`surface_stl_object`](/guide/api/surface-stl-object).

### Reject a file with NaN or Inf coordinates

<<< @/examples/snippets/ck_reject_nan.f90

<<< @/examples/output/ck_reject_nan.ansi{ansi}

With `status=` present the load returns `STATUS_INVALID_INPUT` (5) instead of stopping, and the surface stays empty.
The codes are listed in [constants](/guide/api/constants#status-codes).

### Save as ASCII or binary

<<< @/examples/snippets/ck_save_formats.f90

<<< @/examples/output/ck_save_formats.ansi{ansi}

`is_ascii` defaults to `.true.`. The binary file of the 48-facet cube is `84 + 50 x 48 = 2484` bytes.

### Wall, inlet and outlet from one STL

<<< @/examples/snippets/ck_patches.f90

<<< @/examples/output/ck_patches.ansi{ansi}

![the box of the recipe with its three patches, the wall translucent, the inlet (z = 0) teal and the outlet (z = 1) ember; the three probes joined to their closest points and labelled inlet, outlet and wall](../examples/images/ck-patches.svg){.fossil-figure}

Each `solid <name>` block of an ASCII STL is a named patch, kept on every facet; `load_from_files` reads one file per
patch. The query outputs `patches=` give the boundary condition at each point: the patch of its closest facet. The
patch travels through repair, booleans and remeshing. See [named patches](/guide/advanced/patches).

### Load several files as one surface

<<< @/examples/snippets/ck_load_files.f90

<<< @/examples/output/ck_load_files.ansi{ansi}

![cube-01.stl and naca0012-ascii.stl loaded by one call into one surface: the cube translucent (patch 1, 48 facets), the wing in ember (patch 2, 188 facets) passing through it](../examples/images/ck-load-files.svg){.fossil-figure}

`load_from_files` reads the files in the order given, each in its own format (here one binary, one ASCII), into one
surface: 236 facets in two patches, one per file, and one bounding box over both. A patch takes the `solid` name of
its file, or the file name when that is empty; a multi-solid ASCII file brings all its blocks. The facets keep their
patch through repair and the queries (`facet_patch`, `patches=`). See
[reading and writing patches](/guide/advanced/patches#reading-and-writing).

### Load only a region of a large file

<<< @/examples/snippets/ck_load_clipped.f90

<<< @/examples/output/ck_load_clipped.ansi{ansi}

![the clip box over the dragon: the 1869 facets the clipped load keeps, whole and reaching out of the box, and the rest of the dragon as a ghost](../examples/images/ck-load_clipped.svg){.fossil-figure}

A facet is kept whole when it overlaps the box, so the loaded part reaches below the box (z 13.6403 for a box
starting at 15). To select a region of a surface already in memory, use `clip`
([Split a surface by a box](#split-a-surface-by-a-box)); it keeps whole facets too.

## Transform

### Rotate a body

<<< @/examples/snippets/ck_rotate.f90

<<< @/examples/output/ck_rotate.ansi{ansi}

![the NACA 0012 wing seen along its span before (translucent) and after (stone) a 10-degree nose-up turn about the quarter chord, marked as the centre](../examples/images/ck-rotate.svg){.fossil-figure}

`rotate(axis, angle, center)` turns every facet by `angle` radians, right-handed about `axis`, through `center`
(default: the origin). A negative angle about +z turns the wing clockwise, nose up: its lowest point, now the
trailing edge, drops from y = −0.0598 to y = −0.1315. The facets move, nothing else: call `analyze` before the next
query to rebuild the box, the normals and the tree. `translate`, `mirror` and `resize` follow the same pattern. See
[rotate](/guide/usage#rotate).

### Split a surface by a box

<<< @/examples/snippets/ck_clip.f90

<<< @/examples/output/ck_clip.ansi{ansi}

![the dragon split by a box over its front half: the 3811 facets clip keeps in stone, the 2777 of the remainder translucent, and the box](../examples/images/ck-clip.svg){.fossil-figure}

`clip` keeps the facets that overlap the box and, with `remainder=`, returns the others as a second surface; both are
analyzed. The facets are kept whole, not cut, so the parts reach out of the box, and neither part is watertight: the
cut leaves open boundaries. Together they hold the 6588 facets of the dragon. See [clip](/guide/usage#clip).

## Repair

### Close the holes of a surface

<<< @/examples/snippets/ck_watertight.f90

<<< @/examples/output/ck_watertight.ansi{ansi}

<div class="fossil-pair">

![the cube loaded without its top face, its boundary loop in ember](../examples/images/ck-watertight-open.svg)

![the same box after sanitize with fill_holes: six new facets in ember close the top](../examples/images/ck-watertight-closed.svg)

</div>

`fill_holes=.true.` closes the boundary loops that pass two limits: `max_hole_edges` (default 200) and
`max_hole_size`, the loop's bounding-box diagonal over the surface's (default 0.1). The missing face is a large hole,
so the recipe raises the size limit to 1. `sanitize` reports each defect class on standard error. See
[hole filling](/guide/advanced/hole-filling).

### Check that a surface is valid

<<< @/examples/snippets/ck_validity.f90

<<< @/examples/output/ck_validity.ansi{ansi}

![the NACA 0012 wing section with its 72 boundary edges in ember, along the seams where the two caps meet the side](../examples/images/ck-validity.svg){.fossil-figure}

`is_watertight` asks for two facets on every edge, `is_manifold` only for no edge with three or more, and `is_volume`
for a watertight surface with a positive volume. The wing has edges with a single facet and none with three:
manifold, not watertight. See
[getters and predicates](/guide/api/surface-stl-object#getters-and-predicates).

### Fix an inside-out surface

<<< @/examples/snippets/ck_flip.f90

<<< @/examples/output/ck_flip.ansi{ansi}

<div class="fossil-pair">

![the inside-out cube, a translucent ghost, its 48 facet normals in ember pointing inward](../examples/images/ck-flip-before.svg)

![after sanitize_normals the 48 normals point outward](../examples/images/ck-flip-after.svg)

</div>

`sanitize_normals` makes the windings consistent and points them outward; `analyze` then recomputes the volume.
`sanitize` runs both. See the [repair pipeline](/guide/api/surface-stl-object#repair-pipeline).

### Find self-intersections and remove them

<<< @/examples/snippets/ck_self_intersections.f90

<<< @/examples/output/ck_self_intersections.ansi{ansi}

<div class="fossil-pair">

![two unit cubes merged into one surface, as a translucent ghost, and in ember the closed curve along which their facets cross](../examples/images/ck-self-intersections-before.svg)

![after resolve_self_intersections: 136 facets, the walls inside the other cube removed](../examples/images/ck-self-intersections-after.svg)

</div>

`find_self_intersections` lists every pair of facets that cross, with the segment they share: here 32 pairs along one
closed curve where the walls of the two cubes meet. `resolve_self_intersections` removes them: it keeps the boundary
of the region the surface winds around at least once, watertight, no crossing pair left, volume
1.7270 = 2 − 0.7 × 0.65 × 0.6, the same as the `boolean(BOOL_UNION)` of the two cubes as separate surfaces. See
[self-intersection](/guide/advanced/self-intersection).

### Resolve the self-intersections of one surface

<<< @/examples/snippets/ck_resolve_self.f90

<<< @/examples/output/ck_resolve_self.ansi{ansi}

<div class="fossil-pair">

![the unit cube with its corner pushed through its own bottom face, translucent, the crossing segments in ember](../examples/images/ck-resolve-self-folded.svg)

![what resolve_self_intersections returns: the dent cut at the bottom face, the pocket below it gone, watertight](../examples/images/ck-resolve-self-resolved.svg)

</div>

One shell can pass through itself: here the unit cube with its corner (1, 1, 1) pushed to (0.3, 0.35, −0.3), through
its own bottom face. It is watertight and one shell, so only `find_self_intersections` sees the defect, and there are
no separate bodies to unite. `resolve_self_intersections` keeps the faces with winding number 0 in front and 1 behind:
the pocket below the bottom face, wound −1, goes, and the dent opens through the bottom face. Check the status: on
`BOOL_STATUS_NOT_CLOSED` (an input that is not closed, as `dragon.stl`) or `BOOL_STATUS_CDT_FAILED` (a cut along a
coplanar facet's edge crosses another cut) the surface is left unchanged; [`alpha_wrap`](/guide/advanced/alpha-wrap) is the fallback.
See [self-intersection](/guide/advanced/self-intersection#known-limitations).

### Wrap a triangle soup into a closed surface

<<< @/examples/snippets/ck_alpha_wrap.f90

<<< @/examples/output/ck_alpha_wrap.ansi{ansi}

<div class="fossil-pair">

![the dragon as loaded, 6588 facets, its boundary edges around the holes in ember](../examples/images/ck-alpha-wrap-input.svg)

![the alpha wrap of the dragon from the same camera: 1840 facets, closed, details smaller than the ball smoothed over](../examples/images/ck-alpha-wrap-wrap.svg)

</div>

`alpha_wrap` returns a new surface, watertight and manifold by construction whatever the input: holes, duplicates and
non-manifold edges included. `alpha` is the diameter of the ball that carves the wrap from outside: gaps narrower than
about alpha + 2 offset are closed. Here alpha is 1/20 of the box diagonal and the wrap has 1840 facets; a smaller
alpha follows the input more closely and costs more. See [alpha wrap](/guide/advanced/alpha-wrap).

## Queries

### Tell whether a point is inside

<<< @/examples/snippets/ck_point_inside.f90

<<< @/examples/output/ck_point_inside.ansi{ansi}

![the unit cube as a ghost, p inside it in teal and q outside it in ember, as all four sign algorithms classify them](../examples/images/ck-point_inside.svg){.fossil-figure}

On a watertight surface the four sign algorithms agree. The default is `SIGN_PSEUDO_NORMAL`, the cheapest; see
[sign algorithms](/guide/api/constants#sign-algorithms).

### Find the closest point and its normal

<<< @/examples/snippets/ck_closest_point.f90

<<< @/examples/output/ck_closest_point.ansi{ansi}

![the query point joined to its closest point on the cube's top edge at x = 1, the returned facet in ember and the edge pseudo-normal there](../examples/images/ck-closest_point.svg){.fossil-figure}

The closest point is on the edge between the faces x = 1 and z = 1, so the normal is the edge pseudo-normal, the
angle-weighted mean of the two face normals. See
[distance and inside queries](/guide/api/surface-stl-object#distance-and-inside-queries).

### Tell whether the closest point is on a face, an edge or a vertex

<<< @/examples/snippets/ck_query_point.f90

<<< @/examples/output/ck_query_point.ansi{ansi}

![the unit cube and three query points joined to their feet: one above the top face (face), one off the edge x = 1, z = 1 (edge), one off the corner (1, 1, 1) (vertex)](../examples/images/ck-query-point.svg){.fossil-figure}

`query_point` returns, with the distance and the closest point, the feature the closest point lies on:
`FEATURE_FACE` (1), `FEATURE_EDGE` (2) or `FEATURE_VERTEX` (3), and the pseudo-normal of that feature. Off the
top face the normal is the face's, (0, 0, 1); off an edge it is the mean of the two face normals, (0.707, 0, 0.707);
off a corner it is the angle-weighted mean of the three, (0.577, 0.577, 0.577). An immersed-boundary method uses the
feature to treat ghost points near edges and corners. `query_many` does the same for an array of points. See
[immersed-boundary queries](/guide/advanced/ib-queries#definitions).

### Get a reliable sign on a mesh with holes

<<< @/examples/snippets/ck_sign_holes.f90

<<< @/examples/output/ck_sign_holes.ansi{ansi}

![the cube without its top face and the 588 grid points the pseudo-normal sign misclassifies: all above the opening, all put inside](../examples/images/ck-sign_holes.svg){.fossil-figure}

The pseudo-normal sign is local: near an open boundary it misclassifies 588 of the 8000 points. The fast winding
number counts how much of the surface surrounds the point and misclassifies none. See
[winding number](/guide/advanced/winding-number).

### Measure how much of a surface surrounds a point

<<< @/examples/snippets/ck_winding_number.f90

<<< @/examples/output/ck_winding_number.ansi{ansi}

![the winding number of the open box on the plane x = 0.5, seen side-on: teal (w near 1) inside the box, ember (w near 0) outside, a smooth plume over the opening, the four recipe points labelled with their w](../examples/images/ck-winding-number.svg){.fossil-figure}

The generalized winding number w is the solid angle the surface subtends at the point over 4π: 1 inside a closed
surface, 0 outside, and smooth on an open one. At the centre of the box without a lid w = 0.83333333 = 5/6, the five
faces' share. Across the opening it passes through 1/2 (0.5886 just below, 0.4114 just above), so `w > 0.5` stays a
sensible inside test. `beta=-1` gives the exact O(N) sum; the default Barnes–Hut evaluation agrees with it to the
digits shown. See [winding number](/guide/advanced/winding-number).

### Cast a ray

<<< @/examples/snippets/ck_ray_hit.f90

<<< @/examples/output/ck_ray_hit.ansi{ansi}

![a ray cast along +y from (0, -20, 12) to its first hit on the dragon, the facet hit in ember](../examples/images/ck-ray_hit.svg){.fossil-figure}

`t` is the parameter along the ray, `point = origin + t direction`; the direction need not be a unit vector.
`intersect_ray_all` returns every crossing, once each even through a shared edge or vertex (so the parity of the
count on a closed surface tells inside from outside), and `intersect_ray_any` stops at the first one. A ray cast from
a facet skips it with `exclude=`. See [ray queries](/guide/advanced/ray-queries).

### Choose the acceleration tree

<<< @/examples/snippets/ck_tree_kind.f90

<<< @/examples/output/ck_tree_kind.ansi{ansi}

Every distance query descends a tree of boxes over the facets. `aabb_tree_kind=` on `load_from_file`, `adopt_facets`
or `analyze` picks it: `AABB_TREE_SAH_BVH` (the default, a binary tree built by the surface-area heuristic) or
`AABB_TREE_OCTREE` (the earlier eight-way space partition, kept for comparison). Both find the true closest facet, so
the distances agree bit for bit: the largest difference over the 1000 points is 0. They differ in speed only, and
the SAH tree is the faster on every measured mesh. See
[choosing the AABB tree kind](/guide/usage#choosing-the-aabb-tree-kind).

### Read the acceleration tree

<<< @/examples/snippets/ck_tree_boxes.f90

<<< @/examples/output/ck_tree_boxes.ansi{ansi}

`tree_boxes` returns the boxes of the surface's own tree, breadth first from the root, with the optional depth, leaf
flags, parents and the facets of each leaf (`facet_ids(facet_start(i):facet_start(i+1)-1)`). Every facet is in
exactly one leaf. Above, the leaf boxes that touch a region give the facets a query there would test; the same list
serves to tag the cells of a block or to see where a query spends its time. `max_depth=` stops at a level,
`curve%tree_boxes` reads the segment BVH of a curve. See [reading the tree](/guide/usage#reading-the-tree).

### Run a batch of queries on the GPU

<<< @/examples/snippets/ck_device.f90

<<< @/examples/output/ck_device.ansi{ansi}

`enter_device` copies the BVH arrays to the device once, `query_many_device` takes the same arguments as
`query_many`, and `exit_device` releases the copies. The device kernels are the CPU ones (`_oac` routines), so the
results are bit-identical: 0 mismatches over 8000 points, in distances, facets and foot points. This page is built
with gfortran without OpenACC, so here the device routines run on the host; build with the `tests-nvidia` mode for the
GPU (see the [GPU chapter of the tutorial](./tutorial/16-gpu)). See [device distance](/guide/advanced/device-distance).

## Grids

### Compute the signed distance on a whole grid

<<< @/examples/snippets/ck_grid_distance.f90

<<< @/examples/output/ck_grid_distance.ansi{ansi}

![the signed distance of cube-01 on the grid plane z = 0.5, teal inside and ember outside, with the band nodes framed and the cube's section outlined](../examples/images/ck-grid-distance.svg){.fossil-figure}

The cells near the surface (the band) get exact distances; fast sweeping fills the others with a first-order
solution, -0.4750 at the centre where the exact value is -0.5. For exact values at a few points use `distance_many`.
See [distance field](/guide/advanced/distance-field).

### Check a distance field by the eikonal diagnostic

<<< @/examples/snippets/ck_eikonal.f90

<<< @/examples/output/ck_eikonal.ansi{ansi}

![q = |1 − |∇φ|| of the cube's distance field on the grid plane z = 0.5, cream where the field is a distance and darker where it is not: the flagged nodes form an X along the diagonals inside the cube (its medial axis) and four small clusters off its corners outside](../examples/images/ck-eikonal.svg){.fossil-figure}

`gradient_diagnostic` measures how far the field is from a distance, q = |1 − |∇φ||, by central differences at every
node, and flags the nodes where q > 0.05: 5983 of 68921 here. Inside the cube they lie on the medial axis, where φ has
a kink and q reaches 1.0000: not an error, mask those nodes. The 3352 flagged outside, with q at most 0.1480, show the
first-order error of fast sweeping off the corners. See
[wall distance and field diagnostic](/guide/advanced/ib-queries#wall-distance-and-field-diagnostic).

### Compute the wall distance up to a cutoff

<<< @/examples/snippets/ck_wall_distance.f90

<<< @/examples/output/ck_wall_distance.ansi{ansi}

![the wall distance to cube-01 on the plane z = 0.5 with cutoff 0.25, a flat plateau wherever no wall is closer, and the recipe's first three points; the fourth, (3, 3, 3), is off the frame](../examples/images/ck-wall-distance.svg){.fossil-figure}

The distance is unsigned. A point with no wall closer than the cutoff gets the cutoff and body 0, the centre of the
cube included; the traversal stops early for it. See
[wall distance](/guide/advanced/ib-queries#wall-distance-and-field-diagnostic).

### Tag the cells near the wall

<<< @/examples/snippets/ck_tag_cells.f90

<<< @/examples/output/ck_tag_cells.ansi{ansi}

![one layer of the recipe's cells around cube-01, seen from above: cut cells on the wall, a ring of near cells outside them, inside cells within, outside cells beyond](../examples/images/ck-tag-cells.svg){.fossil-figure}

A box is cut when it meets the surface, near when it is within `radius` of it, otherwise inside or outside by the
sign at its centre (fast winding number by default). The 152 cut cells are the shell between the 6^3 and 4^3 blocks
that contain the cube's faces. See [AMR tagging](/guide/advanced/amr-tagging).

### Tag refinement bands on a block

<<< @/examples/snippets/ck_tag_block.f90

<<< @/examples/output/ck_tag_block.ansi{ansi}

![one cell layer through the middle of the cube, each of its 20 x 20 cells labelled with its tag: 0 on the ring of cut cells along the cube's section, −1 to −4 inward in teal, 1 to 3 outward in ember, and 4 (band + 1) beyond](../examples/images/ck-tag-block.svg){.fossil-figure}

`tag_block` gives every cell of a block its Chebyshev cell distance to the nearest cut cell: 0 on the 488 cut cells,
negative inside, positive outside, capped at band + 1 = 4. Each AMR level refines where |tag| ≤ its band. For a list
of boxes rather than a block, use `box_classify_many` ([above](#tag-the-cells-near-the-wall)). See
[AMR tagging](/guide/advanced/amr-tagging#definitions).

### Refine the cells along sharp edges

<<< @/examples/snippets/ck_feature_boxes.f90

<<< @/examples/output/ck_feature_boxes.ansi{ansi}

![the unit cube, translucent, and the 116 cells of side 0.1 that touch one of its 12 sharp edges, outlined: a frame of cells along every edge](../examples/images/ck-feature-boxes.svg){.fossil-figure}

`build_features` collects the feature edges (sharp above `angle`, here 30 degrees; boundary and non-manifold edges
too, by default) into a tree of their own; `box_feature_many` then tells, for each box, whether a feature edge
touches it. The grid is shifted by 0.03 so that no edge lies on a cell face: each of the 12 edges crosses 11 cells,
and the 8 corners are shared by three edges, so 132 − 16 = 116 cells. The test is conservative (the box is inflated
by 1e-12 of the problem size), so it never misses a contact. `tag_block(feature_band=, feature_tags=)` turns the same
test into tag bands. See [refining along feature edges](/guide/advanced/amr-tagging#refining-along-feature-edges).

### Refine the cells across thin walls

<<< @/examples/snippets/ck_thin_boxes.f90

<<< @/examples/output/ck_thin_boxes.ansi{ansi}

![one layer of 40 x 8 cells of side chord/32 over the NACA 0012 wing seen along the span, the section in ember; the 12 cells box_thin_many flags, dark, at the trailing edge and at the nose](../examples/images/ck-thin-boxes.svg){.fossil-figure}

`thickness` measures each facet; `build_thin` stores the local size min(thickness, gap) of every facet in the tree.
`box_thin_many` flags the boxes that meet a facet thinner than `factor` times the box size: 12 of the 320 cells, from
x = −0.2687 to 0.7625. Most sit along the trailing edge, where the section closes to its blunt base; two sit at the
nose, where the medial ball is bounded by the nose's curvature. A cell of side h cannot resolve those walls; refine
them. See [refining thin walls and gaps](/guide/advanced/amr-tagging#refining-thin-walls-and-gaps).

## Immersed boundaries

### Image points and probes for a ghost-cell method

<<< @/examples/snippets/ck_image_points.f90

<<< @/examples/output/ck_image_points.ansi{ansi}

![an ellipse over its grid of spacing 0.05, zoomed on its left end: the ghost nodes inside the wall in teal, each joined through its foot point on the wall to its mirror image outside in ember, and a probe point 2h out along the normal from each foot](../examples/images/ck-image-points.svg){.fossil-figure}

The ghost points are the grid nodes inside the body within one spacing of the wall, 31 here. `image_points` mirrors
each ghost through its foot point; `probe_many` places a point at a given length along the normal from the foot. Both
say whether the segment from the foot crosses the wall anywhere else (`visible`): on this convex ellipse all 31 are
visible, next to a thin wall or a narrow gap some are not. `surface_stl_object` has the same calls. See
[IB queries](/guide/advanced/ib-queries#definitions).

### Place Lagrangian markers for an IB method

<<< @/examples/snippets/ck_markers.f90

<<< @/examples/output/ck_markers.ansi{ansi}

![the unit cube's edges and the markers surface_markers places on its three faces in view for a grid of spacing 0.1, evenly spread over each face, each with its unit normal](../examples/images/ck-markers.svg){.fossil-figure}

`surface_markers` remeshes a copy of the surface to an edge length of sqrt(3) dx and returns the facet centroids,
areas and normals. The 256 areas add up to 6.000000, the cube's area. Direct forcing wants the markers dx apart:
0.6055 of the spacings lie within [0.8, 1.2] dx (the locked edges of the cube keep the spacing of their splits, so
the markers next to them are farther apart; on a smooth body the fraction is higher). `curve_markers` is the 2-D twin, exact along the arc length. See
[grid-adapted markers](/guide/advanced/ib-queries#grid-adapted-markers).

### Place markers on a curve

<<< @/examples/snippets/ck_curve_markers.f90

<<< @/examples/output/ck_curve_markers.ansi{ansi}

![the NACA 0012 section with the 41 markers curve_markers places for dx = 0.05, evenly spaced along the contour, each with its outward normal](../examples/images/ck-curve-markers.svg){.fossil-figure}

`curve_markers` splits each loop into n = nint(L / dx) pieces of equal arc length and puts a marker at the middle of
each, on the polyline itself, with the length of its piece and the outward normal of its segment. The section,
2.040981 long, gets 41 markers of length 0.049780; 97.56 % of the chord spacings between neighbouring markers are
within [0.8, 1.2] dx, the others at the nose, where the contour turns fastest. See
[grid-adapted markers](/guide/advanced/ib-queries#grid-adapted-markers).

## Cut cells

### Find where a grid line crosses a curve

<<< @/examples/snippets/ck_line_crossings.f90

<<< @/examples/output/ck_line_crossings.ansi{ansi}

![the diamond and the five grid lines v = 0, 0.25, ..., 1, each labelled with its number of crossings: two on v = 0 at the bottom vertex (entering and leaving at the same point), none on v = 1 at the top vertex; entering crossings in teal, leaving ones in ember](../examples/images/ck-line-crossings.svg){.fossil-figure}

`line_crossings` returns the crossings of one grid line with the curve, sorted, each with its segment and whether the
line enters the solid there. The side test is exact, with a symbolic tie-break for points on the line: a grid point
belongs to the low side of the line. The diamond's bottom vertex lies on v = 0, so that line enters and leaves at
u = 0.5 (two crossings, parity kept); its top vertex lies on v = 1, which then misses it (no crossing). Every cell
edge sees the same decision, so fractions and inside flags never disagree. 3-D: `surface%line_crossings`. See
[cut-cell geometry](/guide/advanced/cut-cells).

### Find where the wall cuts the links between cell centres

<<< @/examples/snippets/ck_wall_fractions.f90

<<< @/examples/output/ck_wall_fractions.ansi{ansi}

![a disc of radius 0.3 over 16 x 16 cell centres, the 76 inside it in teal; every link between centres that the wall cuts drawn from both ends to the wall point its fraction gives, yellow](../examples/images/ck-wall-fractions.svg){.fossil-figure}

Ghost-cell and cut-link schemes need, along each grid line, where the wall falls between two cell centres.
`cell_wall_fractions` gives it for every centre and direction (−u, +u, −v, +v; −z, +z too in 3-D): the fraction of
the distance to the next centre at which the line meets the wall, or −1 if it does not. Centre (4, 8) is inside, and
the wall lies 0.2672 of the way to centre (3, 8); seen from (3, 8), outside, it lies 0.7328 of the way back: the two
add up to 1. `inside` is the exact parity; `bad_lines` counts the lines whose crossings do not alternate (0 on a
closed curve). See [cut-cell geometry](/guide/advanced/cut-cells).

### Cut a 3-D grid by a body

<<< @/examples/snippets/ck_cut_cells_3d.f90

<<< @/examples/output/ck_cut_cells_3d.ansi{ansi}

![the turned unit cube, translucent, and one layer of its 16^3 cells, 0.5 < z < 0.625, coloured by fluid fraction: sepia inside the cube, cream outside, the cut cells in between](../examples/images/ck-cut-cells-3d.svg){.fossil-figure}

`cut_cells_3d` cuts the cells of a Cartesian block by a closed surface: the fluid fraction of every cell (`alpha`) and
face (`face_x/y/z`), the fluid centroids, the wetted area, centroid and normal of each cell, and its fluid pieces. The
cube is turned off the grid planes: of the 4096 cells, 554 are cut and 274 solid. The checks need no reference: the
fluid volumes add up to the block's minus the cube's, 7, and the wetted areas to the cube's area, 6, both to the 12
digits printed. Optional outputs: `links=`, `edges=`, `moments_order=`, `quadrature=`. See
[cut-cell geometry](/guide/advanced/cut-cells).

### Moments and quadrature of cut cells

<<< @/examples/snippets/ck_cut_moments.f90

<<< @/examples/output/ck_cut_moments.ansi{ansi}

![cells near the left edge of a disc on a 16 x 16 grid, coloured by fluid fraction, with the disc's polygon and the quadrature points of the 3 x 3 cells around cell (4, 8): Gauss points in the fluid cells, points crowded along the wall in the cut cells, teal for positive and ember for negative weights](../examples/images/ck-cut-moments.svg){.fossil-figure}

`moments_order=` returns the moments of every fluid piece in its cell's scaled coordinates, in the order of
`mono_exponents` (index 4 is X²). `quadrature=.true.` returns a rule per piece, exact for polynomials: the weights add
up to the piece areas within 3.69E-17, and ∫x² over the fluid comes out at 0.2565237066412, 1/3 minus the disc's
moment. The rule uses signed triangles from the cell centre, so a point can fall in the solid with a negative weight.
3-D: `cut_cells_3d(..., moments_order=, quadrature=.true.)`. See
[moments and quadrature](/guide/advanced/cut-solver#moments-and-quadrature).

### Neighbourhoods for a small-cell cure

<<< @/examples/snippets/ck_cut_solver.f90

<<< @/examples/output/ck_cut_solver.ansi{ansi}

![the disc on 16 x 16 cells, the cut cells coloured by their area fraction; the smallest piece, the cells of its state-redistribution neighbourhood outlined, its merge target marked](../examples/images/ck-cut-solver.svg){.fossil-figure}

`cut_links_2d` (3-D: `cut_cells_3d(..., links=)`) links the pieces of neighbouring cells; the solver support builds
on it. The smallest piece holds 4.5 % of its cell. Its state-redistribution neighbourhood grows along its wall
normal (0.58, −0.81) until it holds half a cell: the diagonal piece (+1, −1) comes first, reached through the links
inside the 3 × 3 tile. Its single merge target is the linked face neighbour that best faces the normal. The
weighted factor η^nb = 0.45 is the share it hands to its neighbour. The moments are those of the piece in its
cell's scaled coordinates. See [Cut cells: solver support](/guide/advanced/cut-solver).

### Flux-redistribution weights of a small cell

<<< @/examples/snippets/ck_frd.f90

<<< @/examples/output/ck_frd.ansi{ansi}

![the disc's cut cells near its smallest fluid piece, coloured by fluid fraction; the six cells of the piece's FRD neighbourhood outlined, each labelled with its weight, the small piece itself in ember](../examples/images/ck-frd.svg){.fossil-figure}

Flux redistribution updates a small piece with a stable fraction of its flux and spreads the rest over its
neighbours. `frd_neighbourhoods` gives each piece its neighbours, the pieces reached by monotone paths of links
inside the 3 × 3 block of cells around it, and the weights κⱼ / Σκ, κ the fluid fraction. The smallest piece, 4.52 %
of its cell, gets six members (itself first) whose weights add up to 1; the full fluid cell (12, 4) takes the largest
share, 0.3189. With `z=` the same call works on 3-D cut cells. See
[neighbourhoods for flux redistribution](/guide/advanced/cut-solver#neighbourhoods-for-flux-redistribution).

### Coarsen cut cells

<<< @/examples/snippets/ck_cut_coarsen.f90

<<< @/examples/output/ck_cut_coarsen.ansi{ansi}

<div class="fossil-pair">

![the disc's cut cells on a 32 x 32 grid, coloured by fluid fraction](../examples/images/ck-cut-coarsen-fine.svg)

![the 16 x 16 cut cells cut_coarsen_2d builds from them, the same as cutting the coarse grid directly](../examples/images/ck-cut-coarsen-coarse.svg)

</div>

`cut_coarsen_2d` builds the coarse level of an AMR hierarchy from the fine cut cells and their links, without cutting
the geometry again: the same 204 pieces and 720 links as cutting the 16 x 16 grid directly, the fractions equal to
1.42E-14 and the apertures exactly. A coarse piece is a set of fine pieces linked inside the coarse cell;
`fine_to_coarse` maps each fine piece to it. 3-D: `cut_coarsen_3d`. See
[coarsening](/guide/advanced/cut-solver#coarsening).

### Coarsen 3-D cut cells

<<< @/examples/snippets/ck_cut_coarsen_3d.f90

<<< @/examples/output/ck_cut_coarsen_3d.ansi{ansi}

<div class="fossil-pair">

![the fluid fraction of the turned cube's cut cells on the fine 16^3 grid, layer 0.5 < z < 0.625, seen from above](../examples/images/ck-cut-coarsen-3d-fine.svg)

![the 8^3 cut cells cut_coarsen_3d builds from them, layer 0.5 < z < 0.75](../examples/images/ck-cut-coarsen-3d-coarse.svg)

</div>

`cut_coarsen_3d` builds the cut cells of the grid with every `ratio`-th edge from the fine ones and their links, with
no new geometry: a coarse piece is a set of fine pieces linked inside the coarse cell, its measures the sums. The
3836 fine pieces become 500 coarse ones, and the coarse cells equal the cells cut directly on the coarse grid to
round-off (6.66e-16 in α, 7.63e-17 in wetted area). `fine_to_coarse` maps every fine piece to its coarse one, for
the restriction of a multigrid or AMR solver. See [coarsening](/guide/advanced/cut-solver#coarsening).

## Structured and overset grids

### Locate a point in a curvilinear block

<<< @/examples/snippets/ck_block_locate.f90

<<< @/examples/output/ck_block_locate.ansi{ansi}

![an O-grid of 32 x 8 cells around the unit cylinder, the cell holding a point next to the seam, and the 16 vertices of its width-4 stencil on both sides of the seam](../examples/images/ck-block-locate.svg){.fossil-figure}

`detect_topology_2d` finds the seam the grid writes exactly (faces 1 and 2 abut each other); `validity` certifies
every cell, exactly in 2-D. `locate_many` returns the cell and its local coordinates: the point lies in the first cell
after the seam. The width-4 stencil reaches across the seam through the record (vertex columns 32, 1, 2, 3), and its
weights use the order-matched ξ, so the interpolation is fourth order on the curved grid. 3-D blocks
(`curvilinear_block_t`) work the same way. See [Curvilinear blocks](/guide/advanced/curvilinear-blocks).

### Find which block holds a point

<<< @/examples/snippets/ck_locate_blocks.f90

<<< @/examples/output/ck_locate_blocks.ansi{ansi}

![the two blocks seen from above, the shell's equator (circles r = 1 and 1.5 and its radial lines) and the box's outline; the four points labelled with the block that holds them: 2 inside the sphere, 1 in the shell, 2 in the box's corner, 0 outside both](../examples/images/ck-locate-blocks.svg){.fossil-figure}

`locate_in_blocks` locates every point in a set of blocks, each prepared with `build_search`: the first block that
contains the point wins, so the order of the blocks is a priority. The second point lies in both the shell and the
box and gets the shell (block 1). The first lies in the box only, inside the sphere where the shell has no cells. The
fourth lies in no block: block 0 and `LOC_OUTSIDE` (2). Each found point comes with its cell and local coordinates
ξ, ready for interpolation. See [point location](/guide/advanced/curvilinear-blocks#point-location).

### Locate many points on the GPU

<<< @/examples/snippets/ck_locate_device.f90

<<< @/examples/output/ck_locate_device.ansi{ansi}

`enter_device` maps a block's search arrays to the device once (after `build_search`); `locate_many_device` takes the
arguments of `locate_many`; `exit_device` releases the arrays. The kernels are the CPU ones, so cells, statuses and
local coordinates agree bit for bit: 0 mismatches over 4000 points. Without `enter_device` the call maps the block
itself, around the one call. This page is built without OpenACC, so the device routine runs on the host here. See
[point location](/guide/advanced/curvilinear-blocks#point-location) and the
[GPU chapter of the tutorial](./tutorial/16-gpu).

### Check a block's validity and quality

<<< @/examples/snippets/ck_block_quality.f90

<<< @/examples/output/ck_block_quality.ansi{ansi}

![a 16 x 8 channel block over a bump, cells shaded by equiangle skewness, and the two cells folded by the dragged vertex (marked) outlined in ember](../examples/images/ck-block-quality.svg){.fossil-figure}

`validity` decides the sign of det J in every cell exactly (in 2-D, by orient2d at the corners). Dragging one vertex
past its neighbour folds the two cells next to it: 126 valid, 2 invalid. `quality` measures the rest per cell:
skewness up to 0.930, orthogonality down to −0.978 and scaled min J down to −9.930, the last two negative only in the
folded cells. Skewness alone does not detect a fold: the folded cell (5, 4) has 0.799. 3-D blocks
(`curvilinear_block_t`) also return the warp of each face. See
[cell validity](/guide/advanced/curvilinear-blocks#cell-validity) and [quality metrics](/guide/advanced/curvilinear-blocks#quality-metrics).

### Find how blocks connect

<<< @/examples/snippets/ck_topology.f90

<<< @/examples/output/ck_topology.ansi{ansi}

![two 8 x 8 blocks side by side, the second indexed with i along −y: the abutting interface between them in ember, the periodic outer faces in yellow, each block's i and j directions in teal](../examples/images/ck-topology.svg){.fossil-figure}

`detect_topology_2d` compares the faces of every block with every other, its own included, and appends a record
wherever the vertices coincide exactly (`tolerance=` relaxes this). With `translations=` it also matches faces that
coincide after a period, here (0, 1) and (2, 0). The two blocks get eight records: the shared interface, whose
transform (2, −1) maps block 1's i to block 2's j and its j to block 2's −i, and three periodic pairs. See
[topology records](/guide/advanced/curvilinear-blocks#topology-records).

### Find the seam and the poles of a 3-D block

<<< @/examples/snippets/ck_topology_3d.f90

<<< @/examples/output/ck_topology_3d.ansi{ansi}

![the spherical shell, its inner sphere in stone and its outer sphere translucent; the seam where the faces i = 1 and i = 25 meet in ember and the two collapsed pole lines in yellow](../examples/images/ck-topology-3d.svg){.fossil-figure}

A grid around a sphere closes on itself: the azimuth seam (i = 1 and i = 25 hold the same vertices) and the two poles
(every vertex of the lines j = 1 and j = 13 coincides along i). `detect_topology` finds both by exact coincidence:
faces 1 and 2 abut each other with the identity transform, and faces 3 and 4 are singular, collapsed along direction
1 (i). The records feed the stencils, which reach across the seam, and the validity check, which reports the pole
cells as singular rather than invalid. See [topology records](/guide/advanced/curvilinear-blocks#topology-records).

### Assemble an overset grid

<<< @/examples/snippets/ck_overset.f90

<<< @/examples/output/ck_overset.ansi{ansi}

![a cylinder O-grid over a Cartesian background: the background vertices inside the cylinder are holes, the ring around them receives from the O-grid, the O-grid's outer ring receives from the background](../examples/images/ck-overset.svg){.fossil-figure}

The O-grid's wall becomes the cutter, which blanks the background vertices inside the cylinder (the O-grid itself is
not cut: `cuts`). Every vertex is located in the other block, and `classify_2d` picks the donors: the ring of
background vertices next to the holes and the O-grid's `overlap` boundary receive, from stencils of width 3 that hold
no hole and no receptor. No vertex is left an orphan; `build_connectivity` lists the 89 receptors with their donor
cells and local coordinates, as a flow solver reads them. FOSSIL does not choose the assembly policy: the rule, the
stencil width and the fringe depth are the caller's. See [Overset assembly](/guide/advanced/overset).

### Assemble a 3-D overset grid

<<< @/examples/snippets/ck_overset_3d.f90

<<< @/examples/output/ck_overset_3d.ansi{ansi}

![the plane y = 0 of the assembly seen along y: the shell's grid lines, the background holes inside the sphere in teal, the ring of background receptors around them and the shell's outer ring of receptors in ember, the field vertices faint](../examples/images/ck-overset-3d.svg){.fossil-figure}

The 3-D twin of [Assemble an overset grid](#assemble-an-overset-grid): a spherical shell with its wall on k = 1 and
its outer face (`overlap`) receiving from a Cartesian background. `donor_candidates` locates every vertex in the other
block, with the quality of each candidate; `classify` picks donors by the explicit rule with stencils of width 2:
147 background holes, 314 background receptors around them, and the 325 vertices of the shell's outer face, with no
orphan. `fringe_levels` marks the rings: INNER 1 around the holes, OUTER 1 on the overlap face. See
[donor candidates](/guide/advanced/overset#donor-candidates) and
[selection](/guide/advanced/overset#selection-validity-and-the-fixed-point).

### Find donors next to a curved wall that two blocks discretize differently

<<< @/examples/snippets/ck_wall_projection.f90

<<< @/examples/output/ck_wall_projection.ansi{ansi}

![the plane z = 0.5 near the wall, zoomed: the fine block's grid lines faint, the coarse block's in stone; the fine wall vertices that lie outside the coarse chords in ember, each joined to its shifted search point on the coarse wall in teal](../examples/images/ck-wall-projection.svg){.fossil-figure}

Two blocks that share a curved wall each discretize it by their own chords. On a concave wall the fine block's wall
vertices between two coarse vertices lie outside the coarse block, so they find no donor: 12 of the 13 in the
overlap. `wall_projection` computes, for each wall vertex of block a, the shift onto block b's wall (nudged into the
fluid), decaying to zero over `max_distance` along the grid line; `donor_candidates(shifts=)` searches at the shifted
points. With it all 13 find a donor. See [wall projection](/guide/advanced/overset#wall-projection-for-viscous-walls).

### Find the points a moving body uncovers

<<< @/examples/snippets/ck_exposed_points.f90

<<< @/examples/output/ck_exposed_points.ansi{ansi}

![the background around the moving cylinder: its wall at step n faint and at step n + 1, 0.6 h along x, in stone; the background holes at n + 1 in teal and the 8 exposed points, behind the cylinder, in ember](../examples/images/ck-exposed-points.svg){.fossil-figure}

When a body moves, background vertices that were holes or receptors become field points and need a value.
`exposed_points` compares two classifications and flags them; `violations` counts the holes that became field points
in one step: a hole has no old value to start from, so a step must not do that (Meakin and Suhs). The O-grid moves through
`poses=` without touching its vertices. A step of 0.6 h exposes 8 points, none of them an old hole; a step of 3 h
exposes 45, of which 18 were holes: the step is too long for this fringe. See
[moving grids](/guide/advanced/overset#moving-grids).

### Cut holes with a hole map or X-rays

<<< @/examples/snippets/ck_hole_map.f90

<<< @/examples/output/ck_hole_map.ansi{ansi}

![the nose of the NACA 0012 section over a Cartesian background: the hole-map holes in teal, the field vertices faint, and the background vertices next to the wall that the X-rays misclassify in ember](../examples/images/ck-hole-map.svg){.fossil-figure}

The cutter is a closed curve, the airfoil section, and the background a plain block. `hole_flags_2d` blanks the
vertices whose signed distance (by the winding number) is negative: 522 holes. `hole_map_2d` paints a 16 × 16 map
from its cut elements and settles the vertices in cut elements by the exact crossing parity: the same 522 holes, none
unresolved. `xray_flags_2d` takes the parity on the nearest of 33 rays and gets 13 vertices next to the wall wrong,
all within the ray spacing. See [hole cutting](/guide/advanced/overset#hole-cutting).

### Cut holes around a 3-D body

<<< @/examples/snippets/ck_hole_cutting_3d.f90

<<< @/examples/output/ck_hole_cutting_3d.ansi{ansi}

![the wall cutter of the shell, a translucent sphere, the 595 background vertices the hole map blanks inside it in teal, and the 68 the X-rays get wrong in ember, all near the wall](../examples/images/ck-hole-cutting-3d.svg){.fossil-figure}

`wall_cutter` turns the wall faces of the blocks (records named `wall`) into a closed surface: 528 facets,
watertight. `hole_flags` blanks every vertex whose signed distance to it is negative: 595 background vertices (the
shell is not cut, `cuts=`). `hole_map` paints a 12³ Cartesian map from its cut elements and decides the vertices in
cut elements by the exact crossing parity: the same 595, none unresolved. `xray_flags` takes the parity on the
nearest of 13 × 13 rays along z and gets 68 vertices near the wall wrong. See
[hole cutting](/guide/advanced/overset#hole-cutting).

### Write the overset connectivity for a flow solver

<<< @/examples/snippets/ck_save_connectivity.f90

<<< @/examples/output/ck_save_connectivity.ansi{ansi}

The assembly is the one of [Assemble an overset grid](#assemble-an-overset-grid). `save_connectivity` writes it as
plain text. The header gives each block's field, hole, receptor and orphan counts. Each of the 89 following lines is
one receptor: receptor block and vertex (linear index, i fastest), donor block, donor cell corner (i, j, k), local
coordinates ξ, stencil widths, a low-order flag and the donor quality. See [connectivity output](/guide/advanced/overset#connectivity-output).

### Put a block face on a body with sharp edges

<<< @/examples/snippets/ck_project_face.f90

<<< @/examples/output/ck_project_face.ansi{ansi}

![the unit cube and a 9 x 9 face grid near its face x = 1: before the projection, off the face and wider than it (faint lines); after it, its corners on the cube's corners, its edge lines on the cube's edges and its interior on the face](../examples/images/ck-project-face.svg){.fossil-figure}

`projection_target` collects the body's patches, its feature edges and patch boundaries (as curves between corners)
and its corners, once. `project_block_face` tags each vertex of the face by its distance to them (`tolerance`): the
face's corners go to the cube's corners, its edge lines to the cube's edges, its interior to the face. A plain closest
point would pull grid lines across the edges. See [body-fitted grid support](/guide/advanced/body-fitted).

### Snap points onto a body with sharp edges

<<< @/examples/snippets/ck_project_points.f90

<<< @/examples/output/ck_project_points.ansi{ansi}

![the unit cube, translucent, and the five points of the recipe joined to their feet, each foot labelled with its constraint: face, edge (three) and corner](../examples/images/ck-project-points.svg){.fossil-figure}

`project_block_face` does this for a whole block face; `classify_points` and `project_points` do it for any points.
`classify_points` gives each point a constraint within the tolerance: a corner, else a constraint curve (a sharp edge
or a patch boundary), else the patch of its closest facet. `project_points` moves each point onto its constraint.
The points within 0.12 of the edge x = 1, z = 1 land on it, the third from inside the cube; the last one lands on
the corner. Set `kinds` and `ids` by hand to force a constraint. See
[projection and classification](/guide/advanced/body-fitted#projection-and-classification).

### Extrude an O-grid and check its wall spacing

<<< @/examples/snippets/ck_wall_spacing.f90

<<< @/examples/output/ck_wall_spacing.ansi{ansi}

![an O-grid marched from a NACA0012 section, every second grid line, zoomed on the nose: the lines leave the wall at right angles and spread away from it](../examples/images/ck-wall-spacing.svg){.fossil-figure}

`hyperbolic_grid_2d` marches layers normal to the section, with heights from `geometric_spacing`. The first one is set
for y₁⁺ = 1 with the given u_τ. `block_report_2d` measures what the grid gives: y₁⁺ at every wall vertex, from the
true distance to the wall, the growth ratio along the grid lines, and WP2's quality. See
[body-fitted grid support](/guide/advanced/body-fitted).

### Save the wall-spacing report of a grid

<<< @/examples/snippets/ck_save_report.f90

<<< @/examples/output/ck_save_report.ansi{ansi}

The grid and the report are those of
[Extrude an O-grid and check its wall spacing](#extrude-an-o-grid-and-check-its-wall-spacing). `save_block_report`
writes the totals as `#` lines, then one line per wall point: indices, face, the first spacing along the grid line and
to the wall, y₁⁺ and the growth ratio. The 130 wall points reach y₁⁺ between 1.000 and 1.478, the largest at the
trailing-edge vertex (1, 1). The fraction above the target is 1.0000 because y₁⁺ comes out at 1.00000001, a round-off
over the target of 1. See [the wall-spacing report](/guide/advanced/body-fitted#the-wall-spacing-report).

### Check the first-cell height of a 3-D wall block

<<< @/examples/snippets/ck_block_report_3d.f90

<<< @/examples/output/ck_block_report_3d.ansi{ansi}

![y1+ at the wall points of the mid-span line, as bars along the plate from x = 0.05 to 1.05, and the target y1+ = 1 as a line: above it near the leading edge, below it downstream](../examples/images/ck-block-report-3d.svg){.fossil-figure}

The first height is set for y₁⁺ = 1 at mid-plate, with u_τ from the turbulent flat-plate estimate (`utau_flat_plate`).
`block_report` measures y₁⁺ at every wall point of the 3-D block: with `u_inf=` and `x_le=` it estimates u_τ at each
point from its distance to the leading edge (`utau_estimated` T). u_τ falls along the plate, so y₁⁺ runs from 1.2589
near the leading edge to 0.9285 at the end; 46.34 % of the 205 wall points are above the target. Pass `u_tau=` or
`u_tau_wall=` when a solution gives it. See [the wall-spacing report](/guide/advanced/body-fitted#the-wall-spacing-report).

## Curves

### Cut a 2-D section of an STL

<<< @/examples/snippets/ck_slice.f90

<<< @/examples/output/ck_slice.ansi{ansi}

![the NACA 0012 wing cut by the plane z = 0.5, its airfoil section one closed loop](../examples/images/ck-slice.svg){.fossil-figure}

`slice_surface(surface, normal_axis, level, curve)` cuts by an axis-aligned plane and returns a `curve_object`; the
section of the NACA 0012 wing is one closed loop of 60 vertices with area 0.08165. See
[2-D geometry](/guide/advanced/2d-geometry).

### Compute the signed distance to a polygon

<<< @/examples/snippets/ck_polygon_distance.f90

<<< @/examples/output/ck_polygon_distance.ansi{ansi}

![the signed distance to the unit square, negative inside, with the recipe's three query points joined to their closest points](../examples/images/ck-polygon-distance.svg){.fossil-figure}

Outer loops run counter-clockwise and holes clockwise; the distance is negative inside. See
[`curve_object`](/guide/api/curve-object).

### Second moments of area of a section

<<< @/examples/snippets/ck_mass_properties_2d.f90

<<< @/examples/output/ck_mass_properties_2d.ansi{ansi}

![the hollow 0.1 x 0.2 m section seen from above, its clockwise hole outlined in ember, its centroid at the origin and its radii of gyration along u and v, labelled with the integrals of u² and v²](../examples/images/ck-mass-properties-2d.svg){.fossil-figure}

The hole is wound clockwise, so it subtracts. `second_moments(2)`, the integral of v², is the second moment of area
for bending about the u axis, here (b h³ − b′h′³)/12. See [mass properties](/guide/advanced/mass-properties#2-d-curves).

### Find the corners of a section

<<< @/examples/snippets/ck_curve_corners.f90

<<< @/examples/output/ck_curve_corners.ansi{ansi}

![the trailing edge of the NACA 0012 section, zoomed: the two corners feature_corners finds at the ends of the blunt base in ember, the other vertices faint](../examples/images/ck-curve-corners.svg){.fossil-figure}

`feature_corners` returns the vertices where a closed loop turns by more than `angle` (default 30 degrees) and both
ends of every open loop: the 2-D feature edges. Of the 60 vertices of the section, two qualify: the ends of the
blunt trailing-edge base, 0.00252 thick. The rounded nose turns in many small steps and has none. Grid generators keep
corners as fixed points; `tag_block(feature_band=)` refines around them. See
[2-D corners](/guide/advanced/feature-edges#2-d-corners).

## Geometry processing

### Measure how far a repair moved the surface

<<< @/examples/snippets/ck_deviation.f90

<<< @/examples/output/ck_deviation.ansi{ansi}

![the dragon after sanitize(fill_holes): the facets it kept in stone and the 60 patch facets in ember, over the holes of the original](../examples/images/ck-deviation.svg){.fossil-figure}

The two directions differ: the hole patches are far from the original (they span the holes), while the original lies on
the repaired surface except where `sanitize` removed facets. `lower <= h <= upper` brackets the one-sided Hausdorff
distance (certified, to 1e-3 of the bounding-box diagonal by default); `mean` and `rms` are estimates. See
[surface deviation](/guide/advanced/surface-deviation).

### Compare two versions of a part

<<< @/examples/snippets/ck_hausdorff.f90

<<< @/examples/output/ck_hausdorff.ansi{ansi}

![the wing (translucent) and the same wing at 2 degrees more incidence (stone) seen along the span and zoomed on the trailing edge, where the two are farthest apart: Hausdorff distance 0.0258989](../examples/images/ck-hausdorff.svg){.fossil-figure}

The symmetric Hausdorff distance is the largest distance from a point of either surface to the other: the one number
for "how different are these two parts". `hausdorff` brackets it between a sampled lower bound and a certified upper
bound, here [0.0258979, 0.0258989] for a tolerance of 1e-6. It is a little less than the 0.0261786 the trailing edge
moves: the nearest point of the other wing is not its old trailing edge but its lower surface, closer by. See
[surface deviation](/guide/advanced/surface-deviation).

### Find the sharp edges of a part

<<< @/examples/snippets/ck_feature_edges.f90

<<< @/examples/output/ck_feature_edges.ansi{ansi}

![the dragon's head with the sharp edges above 60 degrees in gold and the non-manifold edges in ember](../examples/images/ck-feature-edges.svg){.fossil-figure}

`feature_edges` reports sharp edges (30 degrees by default), boundary and non-manifold edges, with their corners and
polylines as vertex-pool ids. See [feature edges](/guide/advanced/feature-edges).

### Find walls too thin for my grid

<<< @/examples/snippets/ck_thin_walls.f90

<<< @/examples/output/ck_thin_walls.ansi{ansi}

![the NACA 0012 wing with the 76 facets thinner than h = chord / 32 in ember, most of them along the trailing edge](../examples/images/ck-thin-walls.svg){.fossil-figure}

`thickness(method=THICKNESS_MEDIAL)` gives each facet the diameter of the largest ball tangent to it inside the body
(and outside, for `gap`). `thin_regions` groups the facets thinner than `factor * h` into connected regions. The
thickness is computed once; each cell size is one cheap call. On this NACA 0012 section, the thin facets shrink towards
the trailing edge as h decreases. For AMR, `build_thin` and `tag_block(thin_factor=, thin_band=, thin_tags=)` turn the
same values into refinement tags. See [thickness and gaps](/guide/advanced/thickness).

### Decimate a mesh

<<< @/examples/snippets/ck_decimate.f90

<<< @/examples/output/ck_decimate.ansi{ansi}

![the dragon decimated again and again to half its facets, 6588, 3294, 1646 and 822, then asked for 411 and stopped at 572 with NO_PROGRESS, from the same camera](../examples/images/anim-decimate.svg){.fossil-figure}

Quadric edge collapses remove two facets each, so the result stops at or just below the target (1646 for 1647). No
collapse may change the topology: the dragon is not a manifold (61 non-manifold edges, 83 shells), and halving it a
fourth time stops at 572 facets with `DEC_STATUS_NO_PROGRESS`, a valid surface larger than asked. Check the status.
See [decimation](/guide/advanced/decimate).

### Unite two bodies

<<< @/examples/snippets/ck_union.f90

<<< @/examples/output/ck_union.ansi{ansi}

![the two overlapping unit cubes, then their union, intersection, difference and symmetric difference with each computed volume, the void of the symmetric difference in ember](../examples/images/anim-boolean.svg){.fossil-figure}

Two unit cubes overlapping on a cube of side 0.5 give 2 - 0.125 = 1.875. `BOOL_INTERSECT`, `BOOL_DIFFERENCE` and
`BOOL_SYMDIFF` work the same way. See [booleans](/guide/advanced/booleans).

### Remesh to a uniform edge length

<<< @/examples/snippets/ck_remesh.f90

<<< @/examples/output/ck_remesh.ansi{ansi}

<div class="fossil-pair">

![the unit cube as loaded: 48 facets, edges from 0.5 to 0.71](../examples/images/ck-remesh-before.svg)

![the same cube after isotropic_remesh with target length 0.1 and features preserved: 1336 facets, sharp edges kept](../examples/images/ck-remesh-after.svg)

</div>

`preserve_features=.true.` keeps the twelve sharp edges of the cube, so the volume stays 1. After 10 iterations the
edges are between 0.0639 and 0.1409 for a target of 0.1: the band (0.08, 0.133) except along the locked edges. See [isotropic remeshing](/guide/advanced/isotropic-remesh).

### Measure the curvature at a point

<<< @/examples/snippets/ck_curvature.f90

<<< @/examples/output/ck_curvature.ansi{ansi}

![the marching-cubes unit sphere coloured by the mean-curvature error H - 1 at each facet's foot (teal below 1, ember above, saturating at 0.05), dots where it exceeds 0.1 (none since issues #48 and #51; 25 slivers before), and the recipe's query point above with its foot](../examples/images/ck-curvature.svg){.fossil-figure}

`curvatures=` returns the mean and Gaussian curvature at the foot point, by a cubic fit over the two rings of
vertices around it: on this unit sphere of 6632 facets, H = 0.9944 and K = 0.9892. It runs on the CPU only. See
[curvature at the foot point](/guide/advanced/ib-queries#curvature-at-the-foot-point).

### Curvature at every vertex

<<< @/examples/snippets/ck_vertex_curvature.f90

<<< @/examples/output/ck_vertex_curvature.ansi{ansi}

![the marching-cubes unit sphere twice, coloured by H − 1 per facet: left the cotangent-Laplacian H, mottled teal and ember; right the cubic-jet H, nearly uniform](../examples/images/ck-vertex-curvature.svg){.fossil-figure}

Three per-vertex estimates on a unit sphere (exact H = K = 1) meshed by marching cubes, on a grid with 30 nodes on the
sphere: the extraction has no degenerate facet there (issue #87), and the three return `CURV_STATUS_OK`. `gaussian_curvature` (angle
defect) and `mean_curvature` (cotangent Laplacian) are integral measures divided by a vertex area: on an irregular
mesh they do not converge pointwise, and here they range from 0.1956 to 2.5433 (K) and 0.4169 to 2.6044 (H).
`jet_curvature` fits a cubic over the 2-ring and stays within 0.8243 to 1.1003 (K) and 0.9088 to 1.0492 (H), with an
error O(h²). Use the jet for pointwise values; the angle defect keeps its global meaning, Σ K A = 4π on a closed
genus-0 surface (discrete Gauss–Bonnet). See
[per-vertex curvature](/guide/advanced/curvature).

### Assemble the cotangent Laplacian of a surface

<<< @/examples/snippets/ck_laplacian.f90

<<< @/examples/output/ck_laplacian.ansi{ansi}

`cotangent_laplacian` builds the cotangent Laplacian `L` and the barycentric mass matrix `M` over the surface's
unique vertices, as CSR matrices. Constants are in the kernel of `L` (here exactly: max |L 1| = 0), and the masses add
up to the surface area. FOSSIL provides the matrices; the solve (implicit smoothing, heat-method geodesics) is the
caller's. See [cotangent Laplacian](/guide/advanced/cotangent-laplacian).

### Smooth a noisy surface

<<< @/examples/snippets/ck_smooth.f90

<<< @/examples/output/ck_smooth.ansi{ansi}

<div class="fossil-pair">

![a marching-cubes unit sphere with its vertices moved radially by up to 3 %: a rough, faceted surface](../examples/images/ck-smooth-before.svg)

![the same sphere after 10 Taubin pairs: the roughness gone, the shape kept](../examples/images/ck-smooth-after.svg)

</div>

Taubin smoothing alternates a shrinking step (λ) and an inflating one (μ), so it removes noise without shrinking the
body. Ten pairs bring the rms radial error from 0.02077 to 0.00762 and the area from 14.428 back towards the sphere's
12.566 (12.610); the volume moves by 0.3 %. The grid has nodes on the sphere, and the extraction has no degenerate
facet (issue #87), so `smooth` returns `SMOOTH_STATUS_OK`. See [smoothing](/guide/advanced/smoothing).

### Extract an isosurface from a sampled field

<<< @/examples/snippets/ck_isosurface.f90

<<< @/examples/output/ck_isosurface.ansi{ansi}

![the torus extract_isosurface draws from 31 x 31 x 11 samples of its signed distance, 4208 facets, inside the sampling box](../examples/images/ck-isosurface.svg){.fossil-figure}

`extract_isosurface` runs marching cubes on values sampled at the nodes of a uniform grid between `bmin` and `bmax`
and returns facets for `adopt_facets`; the result is watertight. At spacing 0.1 for a tube radius 0.4 the volume is
1.25 % short of 2π²Rr² (3.11864 for 3.15827): the facets cut the curved surface as chords. Any field works: a
distance, a density, a level set. See [marching cubes](/guide/advanced/marching-cubes).

### Segment a part by its shape diameter

<<< @/examples/snippets/ck_segment.f90

<<< @/examples/output/ck_segment.ansi{ansi}

![the dumbbell coloured by segment_sdf's two clusters: the thin bar (480 facets) in ember, the two balls (4080 facets) in the body colour](../examples/images/ck-segment.svg){.fossil-figure}

`segment_sdf` measures at each facet the shape diameter (the median length of a cone of rays cast inward) and clusters
the values; label 1 is the thinnest cluster. The bar of diameter 0.4 gets a mean of 0.4551, raised by its junctions
with the balls; the balls of diameter 1.2 get 0.9027, since the cone's off-axis rays are chords shorter than the
diameter. Facets whose rays find no opposite wall get label 0 (none here). See
[SDF segmentation](/guide/advanced/sdf-segmentation).

## Moving bodies

### Mass and inertia of a body

<<< @/examples/snippets/ck_mass_properties.f90

<<< @/examples/output/ck_mass_properties.ansi{ansi}

![the 0.1 x 0.2 x 0.3 m block, translucent, its centre of mass, its three principal axes along the edges labelled 0.19625, 0.39250 and 0.51025 kg m², and the corner (0, 0, 0) that the second call takes as about](../examples/images/ck-mass-properties.svg){.fossil-figure}

The body must be closed (`status` says so otherwise). The tensor is about the centre of mass unless `about=` names another
point; the principal moments are sorted ascending, here the block's three diagonal entries. See
[mass properties](/guide/advanced/mass-properties).

### Inertia of a turned body

<<< @/examples/snippets/ck_inertia_world.f90

<<< @/examples/output/ck_inertia_world.ansi{ansi}

![the 0.1 x 0.2 x 0.3 m block before (translucent) and after (stone) a 30-degree turn about z, seen from above, with the principal axes of the turned block from its centre of mass, labelled 0.51025, 0.39250 and 0.19625 (z) kg m²](../examples/images/ck-inertia-world.svg){.fossil-figure}

A rigid body's tensor about its centre of mass moves with it: in the world frame it is R I Rᵀ. `inertia_to_world`
applies a pose's rotation to the body-frame tensor, with no new integration; it matches the tensor integrated on the
turned surface to 4.163e-17. `principal_frame` packs the principal axes and the centre of mass as a pose: its
rotation diagonalizes the tensor (R diag(moments) Rᵀ equals it to 1.110e-16), and its axis of the largest moment is
the block's x axis turned by 30 degrees, (−0.8660, −0.5000) up to sign. See
[frames and poses](/guide/advanced/mass-properties#frames-and-poses).

### How close are my two bodies?

<<< @/examples/snippets/ck_body_distance.f90

<<< @/examples/output/ck_body_distance.ansi{ansi}

![the unit cube and a second unit cube turned 45 degrees about the vertical axis, a corner pointing at the first cube's face x = 1: the closest pair, on the face (yellow) and at the corner (teal), joined by a segment](../examples/images/ck-body-distance.svg){.fossil-figure}

`distance_to` returns the least distance over all pairs of facets, the two facets and the points that realize it, and
`intersecting` when the bodies touch or cross (0, decided exactly). `pose_other=` places the second body without
moving it or rebuilding its tree, here the corner of the turned cube at x = 2 − sqrt(2)/2. For many bodies,
`proximity_pairs(cutoff=)` lists the pairs closer than the cutoff. See [body–body distance](/guide/advanced/body-distance).

### When do my moving bodies first touch?

<<< @/examples/snippets/ck_ccd.f90

<<< @/examples/output/ck_ccd.ansi{ansi}

![the cube at rest, the moving cube at the start of its step (ghost) and at the time of first contact (solid), turned about its own corner, its leading edge touching the resting cube: the witness points on the two cubes](../examples/images/ck-ccd.svg){.fossil-figure}

`time_of_impact` moves each body from its pose at s = 0 to its pose at s = 1 (the motion of `pose_interpolate`: the
rotation about the body's own origin, angle and translation linear in s) and returns the first s at which they come
within the tolerance (default 1e-6 of the larger box diagonal), never later than the true contact: a plate 1e-3 thick
crossed within one step, apart at both ends, is still caught. Bodies sliding along each other just above the
tolerance converge with the default arguments too: the loop switches to the directional bound when its plain steps
stall (`directional=.false.` / `.true.` force either bound). `first_impact` gives the earliest contact among several
bodies. See [continuous collision detection](/guide/advanced/collision#the-default-adaptive-issue-75).

### Which of my moving bodies touch first?

<<< @/examples/snippets/ck_first_impact.f90

<<< @/examples/output/ck_first_impact.ansi{ansi}

![four unit cubes on a row seen from above at the first contact: cube 2 sliding up out of the row, cube 4 turning towards cube 3 and touching it with a corner; the pair 3, 4 in ember, the moving cubes at the start of the step translucent](../examples/images/ck-first-impact.svg){.fossil-figure}

`first_impact` finds the earliest contact among many bodies over a step, each moving from `poses0(i)` to
`poses1(i)`. A broad phase on boxes grown by each body's motion bound drops the pairs that cannot meet, and the rest
run `time_of_impact` in the order of their lower bounds, stopping when the next bound is later than the best contact.
Here 1 of the 6 pairs runs it: cubes 3 and 4 touch at s = 0.440563. The result equals the loop over all pairs, ties
going to the smallest pair. `first_impact_2d` does the same for curves; the squares give the same s. See
[several bodies](/guide/advanced/collision#several-bodies).

### Query a rotated body without rebuilding it

<<< @/examples/snippets/ck_pose.f90

<<< @/examples/output/ck_pose.ansi{ansi}

![the cube rotated 45 degrees by the pose, seen from above over the translucent cube at rest: the query point sits just inside the rotated cube's corner, 0.2 outside the cube at rest](../examples/images/ck-pose.svg){.fossil-figure}

With `pose=` the world-frame point is mapped into the body frame, and the tree built at load time serves every pose.
`pose_from_axis_angle` rotates about the origin, so the translation sets the pivot. A corner of the rotated cube reaches x = 0.5 +
sqrt(2)/2, so the point is inside it, 0.0050 from a face. See [moving bodies](/guide/advanced/moving-bodies).

### Swept volume of a moving body on a grid

<<< @/examples/snippets/ck_swept_cells.f90

<<< @/examples/output/ck_swept_cells.ansi{ansi}

![a 0.3 x 0.1 plate at its start pose (blue outline), aligned with the grid, and at its end pose (yellow outline), turned 30 degrees and moved; the cells coloured by the swept area: ember where the plate left and fluid was gained, teal where it arrived and fluid was lost](../examples/images/ck-swept-cells.svg){.fossil-figure}

`swept_cells_2d` moves the body from `pose0` to `pose1` by screw motion and returns, per cell, the fluid area change,
the time-averaged face apertures and wall, and the cells cleared or covered, integrated exactly for translation and
rotation alike. The body keeps its area, so the 0.01827 gained equals the 0.01827 lost (sum −2.33E-17); the volume
balance residual, at most 3.28E-17, checks that the faces and the swept areas agree. 3-D:
`swept_cells_3d(surface, ...)`. See [swept geometry](/guide/advanced/moving-bodies#swept-geometry).
### Swept volume of a moving body on a 3-D grid

<<< @/examples/snippets/ck_swept_cells_3d.f90

<<< @/examples/output/ck_swept_cells_3d.ansi{ansi}

![a 0.3 x 0.1 x 0.2 box seen from above at its start pose (translucent) and its end pose (stone), turned 20 degrees about z and moved, over one layer of cells coloured by swept volume: ember where it left and fluid was gained, teal where it arrived and fluid was lost](../examples/images/ck-swept-cells-3d.svg){.fossil-figure}

The 3-D twin of the 2-D recipe above: `swept_cells_3d` gives the change of fluid volume of every cell over the step
(from the cut cells at both poses, exact), the time-averaged face apertures (exactly integrated here: `exact` T) and
the wall term. The box is a closed body, so the fluid gained, 0.004102, equals
the fluid lost and the swept volumes add up to round-off (−8.23e-18); the volume balance of the rigid field holds to
6.88e-18 in every cell. 2 cells are cleared (no fluid at the start, fluid at the end). See
[swept geometry](/guide/advanced/moving-bodies#swept-geometry).

