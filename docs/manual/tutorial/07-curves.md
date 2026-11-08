# 7. Curves in 2-D

A 2-D run of the same case needs the body as a curve in a plane, not as triangles. `slice_surface` cuts the dragon
by a plane, and the section is a `curve_object` with the same queries as a surface, in two dimensions. This chapter
slices, queries the section and builds its distance field.

## Slicing the surface

<<< @/examples/snippets/tut07_a_slice-slice.f90

- `slice_surface(surface, normal_axis, level, curve)` cuts by the plane normal to axis `normal_axis` (1, 2, 3) at
  coordinate `level`. The curve lives in that plane with in-plane coordinates (u, v), the other two axes in cyclic
  order: for `normal_axis = 3` they are (x, y). The crossings are exact constructions, with no tolerance: a vertex on
  the plane counts as above it.
- `merge_collinear=.true.` drops vertices between collinear segments. `open_loops` counts the chains that did not
  close, as happens on a surface with holes; such a curve can be closed by `repair(tolerance=)`.
- The solid is on the left of every segment: outer loops run counter-clockwise, holes clockwise. `check_validity()`
  returns `CURVE_VALID` or the first defect found (`CURVE_OPEN_LOOP`, `CURVE_WRONG_ORIENTATION`, ...).
- `get_area` and `get_perimeter` are the 2-D metrics. `save_into_file` writes the ASCII loop format, one `u v` per
  line, read back by `load_from_file`. Chapters 10 and 11 read `dragon-section.dat`.
- A curve can also be built from arrays: `adopt_loops(u=, v=, loop_sizes=)`.

<<< @/examples/output/tut07.ansi{ansi}

The section has two closed loops, 152 vertices and an area of 0.072609.

![the placed dragon cut by the plane z = 0, drawn as a rectangle: the part below the plane solid, the part above translucent, the section curve in ember along the cut](../../examples/images/tut07-section.svg){.fossil-figure}

Each level gives another curve. The second program slices the dragon by 13 planes, from z = -0.30 to 0.30:

<<< @/examples/output/tut07-planes.ansi{ansi}

![a horizontal plane sweeping up through the translucent dragon from z = -0.30 to 0.30 in steps of 0.05, one frame per plane, the section curves in ember on it: 8 loops near the bottom, 2 at z = 0, two small ones at the top; the caption gives the loops and the area](../../examples/images/anim-tut07-planes.svg){.fossil-figure}

No section is open: `sanitize` left no boundary edge (chapter 2). The areas, summed times the spacing 0.05, give
0.053, close to the volume 0.054283 of chapter 3.

## A fold of the section

`check_validity()` does not look for crossings by default; `check_validity(check_self_intersections=.true.)` does
(`CURVE_SELF_INTERSECTING`), and `self_intersections(pairs, n)` returns the pairs of segments that cross, touch or
overlap, decided exactly. The section at z = 0 is valid and still crosses itself once, at a fold of the non-manifold
dragon: segment 12 doubles back along segment 11, and 13 crosses 11.

![the section z = 0 seen from above, zoomed on its fold: segment 11 and segment 13 in ember crossing, segment 12 doubling back between them, the rest of the loop in gold with its vertices](../../examples/images/tut07-crossing.svg){.fossil-figure}

The distance queries do not mind the fold; the cut cells and the swept geometry do ([chapter 11](./11-moving-bodies)
shows how the swept-geometry residual detects it).

## Distance to a curve

<<< @/examples/snippets/tut07_c_query-query.f90

- Points are `p = [u, v]`. `compute_distance` returns the signed distance, negative inside. Unlike the surface
  queries, curve distances are signed and not squared by default (`is_signed=.false.` for the unsigned one). It
  comes with the closest point, `closest_segment`, and `closest_vertex` (0 when the closest point is inside a
  segment).
- `distance(p)` and `is_point_inside(p)` are the short forms. `distance_many(points=, distances=)` is the batch, with
  the points as columns of a `(2, n)` array.

<<< @/examples/output/tut07-query.ansi{ansi}

The origin is inside, 0.015933 from the wall, its closest point inside segment 20. Three of the twelve points of the
ellipse are inside too.

![the section z = 0 seen from above, x pointing down and y to the right: the origin and twelve points on an ellipse around it, in teal inside and ember outside, each joined by a gold segment to its closest point on the two loops](../../examples/images/tut07-query.svg){.fossil-figure}

## A 2-D field

<<< @/examples/snippets/tut07_d_field-field.f90

- `initialize(x=, y=)` without `z` makes a 2-D field; `compute(curve=)` runs the same band, sign and sweep stages as
  in chapter 6.

<<< @/examples/output/tut07-field.ansi{ansi}

The 2-D field has 327 nodes inside on a spacing of 0.015, an area of 0.0736, against 0.0726 for the curve.

![the 2-D signed distance of the section seen from above, x pointing down and y to the right: the nodes of the 81 by 81 field around the curve coloured by tanh(phi/0.05), teal inside, cream on the wall, ember outside, with the two loops of the curve in yellow (the larger one runs through the middle and right lobes)](../../examples/images/tut07-field.svg){.fossil-figure}

The colour saturates beyond about 0.1 from the wall; the field itself runs to the edges of the domain.

::: tip What you learned
`slice_surface` turns a surface into a `curve_object`; curves have the signed distance, the inside test and the
distance field of surfaces, in (u, v) coordinates, and `self_intersections` finds where a section folds. Reference:
[2-D geometry (curves)](/guide/advanced/2d-geometry), [curve_object](/guide/api/curve-object).
:::

Next: [8. Immersed-boundary queries](./08-ib-queries).
