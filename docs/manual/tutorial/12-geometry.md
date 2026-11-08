# 12. Geometry processing and analysis

Not every body arrives as a scanned mesh. A machined part is built from primitives, and before it goes to a solver
its surface is remeshed, checked against the design, and measured: where its walls are thinner than a cell, which
faces carry which boundary condition, what it weighs. This chapter does it on a bracket: a base, a round boss on top,
a bore drilled close to one end. The program builds everything itself from `cube-01.stl`.

<<< @/examples/snippets/tut12_geometry-build.f90

- `box` and `cylinder` are helpers of the program (after `contains`): `box` resizes and moves `cube-01.stl`,
  `cylinder` builds a 32-sided prism facet by facet and passes it to `adopt_facets`. Both name their facets with
  `set_patch(facets=, name=)`, so every facet of the part remembers where it came from.
- `boolean(other=, op=, status=)` replaces the surface by the union (`BOOL_UNION`), intersection, difference
  (`BOOL_DIFFERENCE`) or symmetric difference with `other`. Both must be closed. The arrangement is exact: cut points
  are implicit points and every decision is an exact predicate, so the result is watertight with no tolerance to tune.
  The patch tables of the two inputs are merged by name.

## Remeshing

A boolean result is correct but badly shaped: long facets from the box faces, fans from the prism caps. A solver
wants edges of about one length.

<<< @/examples/snippets/tut12_geometry-remesh.f90

- `isotropic_remesh(target_length=, iterations=, preserve_features=, status=)` splits edges longer than 4/3 of the
  target, collapses those shorter than 4/5, flips edges towards valence 6 and relaxes the vertices, projecting them back
  onto the input surface. `preserve_features=.true.` locks the vertices on creases sharper than 30 degrees, so the
  box edges and the rims stay where they are. Patch boundaries are locked too.
- The remesh works on a copy (`mesh = part`): the part stays as the reference for the checks that follow.

## Checking the change

A remeshed surface is a different surface. How far it moved is a number with a guarantee, not an impression.

<<< @/examples/snippets/tut12_geometry-deviation.f90

- `deviation(other, lower=, upper=, ...)` brackets the one-sided Hausdorff distance h(self, other), the largest
  distance from a point of self to other: `lower` from samples of self (vertices, edge midpoints, facet interiors),
  `upper` certified by a branch and bound that splits the facets of self until `upper - lower` is below `tolerance`
  (default 1e-3 of the bounding-box diagonal). `per_vertex` gives the distance at each vertex; `mean` and `rms` are
  also available.
- The two directions answer two questions: h(mesh, part) is how far the mesh strays from the part, h(part, mesh)
  what of the part the mesh misses (a cut corner, a lost rim). `hausdorff(other)` is the larger of the two.
- Where two surfaces share a plane, their facets cross in it and the branch and bound refines along the crossings,
  at a cost of about 1/tolerance. The part and the mesh share every flat face, so the second call takes
  `tolerance=5.e-3_R8P`, which keeps the run short.

## Analysing the part

<<< @/examples/snippets/tut12_geometry-patches.f90

- `get_patches_number()`, `patch_name(i)` and `facet_patch(i)` read the patch table. It also holds the empty
  patches the inputs were made with (the header of `cube-01.stl`, the unnamed prisms), which the loop skips. The
  ASCII writer saves one `solid <name>` block per non-empty patch, and the queries return the patch of the closest
  facet (`query_many(patches=)`), which is how a solver assigns boundary conditions.

<<< @/examples/snippets/tut12_geometry-features.f90

- `feature_edges(edges=, kinds=, corners=, polyline_start=, polyline_vertices=, angle=, status=)` returns the creases
  sharper than `angle` (default 30 degrees; `FEATURE_EDGE_SHARP`), the boundary and non-manifold edges, the corners
  where they meet, and the edges chained into polylines between corners (CSR arrays: polyline p has the vertices
  `polyline_vertices(polyline_start(p):polyline_start(p+1)-1)`).
- Comparing them before and after the remesh is a check the deviation cannot make: a folded facet in a flat face is
  at distance zero from the part, but it adds a crease.

![the remeshed bracket coloured by patch, the base in sepia, the boss in teal, the bore in ember, its 292 feature edges in gold along the 12 edges of the base and the rims of the boss and the bore, and its 8 corners, the corners of the base](../../examples/images/tut12-part.svg){.fossil-figure}

<<< @/examples/snippets/tut12_geometry-thickness.f90

- `thickness(thickness=, gap=, method=, invalid=, status=)` gives per facet how thick the body is behind it, and the
  gap to the next wall in front. `THICKNESS_RAY` is the median of a cone of 31 rays (shape diameter);
  `THICKNESS_MEDIAL` the diameter of the largest inscribed ball. `invalid` counts the facets where too few rays
  counted, which flags an open or misoriented surface.
- `thin_regions(thickness=, gap=, h=, factor=, region=, region_area=)` groups the facets thinner than `factor * h`
  into edge-connected regions: `region(i)` is the region of facet i (0 when it is not thin), `region_area(r)` the
  area of region r. With h the cell size of the grid to come, these are the places where a wall falls between cells.

![the remeshed bracket seen from its x = 0 end, each facet coloured by its ray thickness from teal (thin) to ember (0.5 and more), the facets thinner than h = 0.1 outlined in gold on the end face over the web left by the bore (the second region, the bore wall facing it, is hidden inside the bore)](../../examples/images/tut12-thickness.svg){.fossil-figure}

<<< @/examples/snippets/tut12_geometry-mass.f90

- `mass_properties(mass=, center=, inertia=, principal_moments=, principal_axes=, density=, about=, status=)` gives
  the exact mass properties of a closed body (`STATUS_INVALID_INPUT` otherwise), as in
  [chapter 11](./11-moving-bodies). With the density of aluminium in SI units, the lengths are metres.

## Running it

<<< @/examples/output/tut12.ansi{ansi}

- The booleans give a watertight part of 524 facets and volume 1.00297; the remesh 2218 facets with no edge longer
  than 0.3057, still watertight.
- The vertices of the mesh lie on the part (3.92e-13 at most: the remesh projects them onto it). The facets between
  them stray by at most 3.56e-3, and the part is at most 6.15e-3 from the mesh. Both lower bounds, about 1.2e-3, are
  of the order of the sagitta of the 32-sided prisms, r(1 − cos(π/32)) = 1.7e-3 for the boss: the mesh facets do not
  follow every kink of the polygons. The volume changes by 4.94e-5 of itself.
- The three patches came through both booleans and the remesh: base, boss, bore.
- The part has 260 feature edges, 8 corners and 16 polylines: the 12 edges of the base, and the top and bottom rims
  of the boss and of the bore, which have no corners (their facets meet at 11.25 degrees, below the threshold). The
  mesh has more feature edges, 292, because the remesh split them, but the same 8 corners and 16 polylines.
- Two regions are thinner than h = 0.1, both on the web between the bore and the x = 0 end: the end face (14 facets)
  and the bore wall facing it (62 facets). The thinnest value, 0.0646, is the web of 0.06 seen through the cone:
  on a plate of thickness t the ray method gives 1.072 t (see [Thickness](/guide/advanced/thickness)). A grid with
  cells of 0.1 does not resolve this web.
- In aluminium the bracket weighs 2707.87 kg, with its centre of mass at x = 1.1704, towards the boss and away from
  the bore, and on the plane of symmetry y = 0.5.

The first figure shows the limit of this remesh: fans of thin facets remain around the corners of the base and
along the rims, where the vertices are locked.

::: tip What you learned
`boolean` builds a watertight part from closed primitives, `set_patch` names its faces, `isotropic_remesh` with
`preserve_features` evens it out, `deviation` certifies how far the result moved, and `feature_edges`, `thickness`
with `thin_regions`, and `mass_properties` measure it. Reference: [Boolean operations](/guide/advanced/booleans),
[Isotropic remeshing](/guide/advanced/isotropic-remesh), [Surface deviation](/guide/advanced/surface-deviation),
[Named patches](/guide/advanced/patches), [Feature edges](/guide/advanced/feature-edges),
[Thickness and gaps](/guide/advanced/thickness), [Mass properties](/guide/advanced/mass-properties).
:::

Next: [13. Structured blocks and overset grids](./13-blocks-overset).
