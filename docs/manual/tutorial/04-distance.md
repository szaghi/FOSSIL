# 4. Signed distance

With the dragon placed in the domain, the first question a grid asks is how far each point is from the wall, and on
which side. This chapter queries single points.

<<< @/examples/snippets/tut04_distance-distance.f90

- `distance(point=)` returns the distance from a point to the closest facet, found by a best-first traversal of the
  AABB tree. The defaults are unsigned and **squared**: pass `is_square_root=.true.` for the distance itself and
  `is_signed=.true.` for a sign, negative inside.
- `sign_algorithm=` chooses how the sign is decided:
  - `SIGN_PSEUDO_NORMAL` (the default): the angle-weighted pseudo-normal of the closest face, edge or vertex, read in
    the same traversal. Exact on a closed, outward-oriented surface.
  - `SIGN_FAST_WINDING_NUMBER`: inside when the winding number of the surface around the point exceeds 0.5, evaluated
    hierarchically on the tree. It tolerates holes and gaps.
  - `SIGN_RAY_INTERSECTIONS` and `SIGN_SOLID_ANGLE`: ray parity and the solid-angle sum, O(N) per query.
- `is_point_inside(point=, sign_algorithm=)` returns the inside test alone, with the same algorithms and the same
  default.

<<< @/examples/snippets/tut04_distance-rich.f90

- `compute_distance` is the subroutine form with the details of the closest feature: `facet_index`, the local
  `edge_index` (1 to 3, 0 when not on an edge) and `vertex_index` (1 to 3, 0 when not at a vertex), the
  `closest_point` and the pseudo-normal `normal` of that feature. Chapter 8 gets the same for many points at once.

## Running it

<<< @/examples/output/tut04.ansi{ansi}

The five points on the y axis go from outside (0.150398) through the body to outside again (0.114892). The two signs
agree at these points. The point (0, 0.6, 0) is closest to vertex 3 of its facet, so its normal is the vertex
pseudo-normal, not the facet normal. On this open dragon the two sign algorithms do not agree everywhere: chapter 5
counts where they differ.

::: tip What you learned
`distance` is unsigned and squared unless asked otherwise; `sign_algorithm` picks the inside test;
`compute_distance` adds the closest facet, feature, point and normal. Reference:
[Generalized winding number](/guide/advanced/winding-number), [Sign algorithms at sharp features](/guide/advanced/sign-benchmark).
:::

Next: [5. Batched queries](./05-batched).
