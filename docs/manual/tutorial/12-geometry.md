# 12. Geometry processing and analysis

Not every body arrives as a scanned mesh. A machined part is built from primitives, and before it goes to a solver
its surface is remeshed, checked against the design, and measured: which faces carry which boundary condition, where
its walls are thinner than a cell, what it weighs. This chapter does it on a bracket: a base, a round boss on top, a
bore drilled close to one end. Every program builds the part itself from `cube-01.stl`, in six steps, each a short
run that starts where the previous one ended.

![the bracket built step by step: the base, the boss prism, their union, the bore prism, the difference, the remesh after 2 and 6 iterations, its feature edges and corners, its ray thickness, its principal axes of inertia](../../examples/images/anim-tut12-bracket.svg){.fossil-figure}

## Building the part

The primitives are closed surfaces. The base is `cube-01.stl` resized to 2 × 1 × 0.5; the boss and the bore are
32-sided prisms built facet by facet.

<<< @/examples/snippets/tut12_a_build-primitives.f90

`box` and `cylinder` are helpers of the program, after its `contains`. The prism:

<<< @/examples/snippets/tut12_a_build-cylinder.f90

- `adopt_facets(facets=)` makes a surface from an array of `facet_object`: it builds the vertex pool, the
  connectivity and the tree, as `load_from_file` does. The caps are fans, outward by the vertex order.
- `set_patch(facets=, name=)` puts the facets in a named patch, created if new. Every facet keeps its patch through
  the operations that follow, so the finished part still knows which face came from which primitive.

The boss goes on with a union:

<<< @/examples/snippets/tut12_a_build-union.f90

- `boolean(other=, op=, status=)` replaces the surface by the union (`BOOL_UNION`), intersection (`BOOL_INTERSECT`),
  difference (`BOOL_DIFFERENCE`) or symmetric difference (`BOOL_SYMDIFF`) with `other`. Both must be closed and
  outward. The patch tables of the two are merged by name.
- The arrangement is exact: cut points are implicit points, every decision is an exact predicate, and there is no
  tolerance to tune. The program checks it against the volume the polygonal solids must have, 1 + A₃₂ · 0.4 with
  A₃₂ = 16 r² sin(2π/32) the area of the 32-gon.

The bore comes out with a difference:

<<< @/examples/snippets/tut12_a_build-difference.f90

<<< @/examples/output/tut12-build.ansi{ansi}

- Each primitive is watertight, and so is each result: 248 facets after the union, 524 after the difference.
- The volumes differ from the exact polygonal volumes by −8.88e-16 and −1.55e-15: round-off of the volume sum, no
  error from the cuts.
- The part carries three patches, base, boss and bore. The figure shows them, and the cut curves in gold: the edges
  where a facet of one patch meets a facet of another, the rims where the boss stands on the base and where the
  bore enters and leaves it.

![the bracket after the union and the difference, 524 facets coloured by patch, the base in sepia with long fans of facets from its six faces, the boss in teal, the bore in ember, the 204 cut edges between patches in gold along the rims of the boss and the bore](../../examples/images/tut12-build.svg){.fossil-figure}

## Remeshing

A boolean result is correct but badly shaped: the faces of the base are a few long triangles fanned out to the
rims, and the edges run from 0.0038 to 1.118. A solver wants edges of about one length.

<<< @/examples/snippets/tut12_b_remesh-remesh.f90

- `isotropic_remesh(target_length=, iterations=, preserve_features=, status=)` runs, per iteration, four passes:
  split the edges longer than 4/3 of the target, collapse those shorter than 4/5, flip edges towards valence 6, and
  relax the vertices tangentially, projecting them back onto the input surface.
- `preserve_features=.true.` locks the vertices on creases sharper than 30 degrees, so the box edges and the rims
  stay where they are. Patch boundaries are locked too.
- The remesh works on a copy (`mesh = part`): the part stays as the reference for the checks of the next steps.

Then the same remesh towards a longer target, 0.4, on the result:

<<< @/examples/snippets/tut12_b_remesh-coarsen.f90

<<< @/examples/output/tut12-remesh.ansi{ansi}

- The facets grow from 524 to 2218 and the longest edge falls from 1.118 to 0.3057 over six iterations. That is
  still above 4/3 of the target (0.267): the long edges of the base need more passes.
- The shortest edge stays at 0.0038. Edges shorter than 4/5 of the target are not collapsed: the collapse test
  rejects almost every interior edge ([issue #79](https://github.com/szaghi/FOSSIL/issues/79)). The remesh refines,
  it does not coarsen. The last line shows it: the remesh towards 0.4 leaves the 2218 facets as they are.

![the boolean result on the left, 524 facets in long fans, and its remesh on the right, 2218 facets: finer and more even on the faces of the base, with fans of thin facets left around the locked corners and rims](../../examples/images/tut12-remesh.svg){.fossil-figure}

The figure shows the limit of this remesh: fans of thin facets remain around the corners of the base and along the
rims, where the vertices are locked and the short edges cannot collapse.

## Checking the change

A remeshed surface is a different surface. How far it moved is a number with a guarantee, not an impression.

<<< @/examples/snippets/tut12_c_deviation-deviation.f90

- `deviation(other, lower=, upper=, mean=, rms=, per_vertex=, tolerance=, status=)` brackets the one-sided Hausdorff
  distance h(self, other), the largest distance from a point of self to other. `lower` is the largest distance at
  samples of self (vertices, edge midpoints, facet interiors); `upper` is certified by a branch and bound that splits
  the facets of self until `upper - lower` is below `tolerance` (default 1e-3 of the bounding-box diagonal).
  `per_vertex` is the distance at each vertex, `mean` and `rms` are estimates over the surface.
- The two directions answer two questions: h(mesh, part) is how far the mesh strays from the part, h(part, mesh)
  what of the part the mesh misses (a cut corner, a lost rim). `hausdorff(other)` is the larger of the two.
- Where two surfaces share a plane their facets cross in it, and the branch and bound refines along the crossings at
  a cost of about 1/tolerance. The part and the mesh share every flat face, so the second call takes
  `tolerance=5.e-3_R8P`, which keeps the run short.

<<< @/examples/output/tut12-deviation.ansi{ansi}

- The vertices of the mesh lie on the part, 3.92e-13 at most: the remesh projects them onto it. The facets between
  them stray by at most 3.56e-3 (certified), at least 1.15e-3 (sampled); the mean distance is 7.22e-6.
- The part is at most 6.15e-3 from the mesh. Both lower bounds, about 1.2e-3, are of the order of the sagitta of the
  32-gon of the boss, r(1 − cos(π/32)) = 1.7e-3: a mesh facet that spans a kink of the polygon cuts the corner.
- The volume changes by −4.94e-5 of itself.

The figure colours each facet of the mesh by the largest distance to the part at seven points of it. The distance is
zero on the flat faces and largest on needle facets of the boss wall that bridge two sides of the 32-gon; the
largest value, 1.15e-3, is the lower bound of the run.

![the boss of the remeshed bracket, each facet coloured from cream (on the part) to ember (1.15e-3 away): the flat faces cream, ember streaks on the side wall of the boss where needle facets bridge the kinks of the 32-gon, the largest marked](../../examples/images/tut12-deviation.svg){.fossil-figure}

## Patches and feature edges

A solver assigns boundary conditions by patch, and a grid generator keeps its lines on the creases.

<<< @/examples/snippets/tut12_d_features-patches.f90

- `get_patches_number()`, `patch_name(i)` and `facet_patch(i)` read the patch table. The table also holds the empty
  patches the inputs were made with (the header of `cube-01.stl`, the unnamed prisms), which the loop skips.
- The ASCII writer saves one `solid <name>` block per non-empty patch, and the point queries return the patch of the
  closest facet (`query_many(patches=)`, [chapter 8](./08-ib-queries)).

<<< @/examples/snippets/tut12_d_features-features.f90

- `feature_edges(edges=, kinds=, corners=, polyline_start=, polyline_vertices=, angle=, status=)` returns the creases
  sharper than `angle` (default 30 degrees; `FEATURE_EDGE_SHARP`), the boundary and non-manifold edges, the corners
  where they meet, and the edges chained into polylines between corners. The polylines are CSR arrays: polyline p has
  the vertices `polyline_vertices(polyline_start(p):polyline_start(p+1)-1)`; a polyline whose first vertex is its last
  is a closed loop.
- Comparing them before and after the remesh is a check the deviation cannot make: a folded facet in a flat face is
  at distance zero from the part, but it adds a crease.

<<< @/examples/output/tut12-features.ansi{ansi}

- The three patches came through both booleans and the remesh: base 1324 facets, boss 454, bore 440.
- The part has 260 feature edges, 8 corners and 16 polylines: 12 chains between corners, the edges of the base, and 4
  closed loops without a corner, the top and bottom rims of the boss and of the bore. Their facets meet at 11.25
  degrees, below the threshold, so the rims have no corners.
- The mesh has more feature edges, 292, because the remesh split them, but the same 8 corners and 16 polylines, all
  sharp: no crease was added or lost.

![the remeshed bracket coloured by patch, the base in sepia, the boss in teal, the bore in ember, its 292 feature edges in gold along the 12 edges of the base and the rims of the boss and the bore, and its 8 corners, the corners of the base](../../examples/images/tut12-part.svg){.fossil-figure}

## Thin walls

The bore is 0.06 from the x = 0 end. On a grid of cells of 0.1, that web falls between cells.

<<< @/examples/snippets/tut12_e_thickness-thickness.f90

- `thickness(thickness=, gap=, method=, invalid=, status=)` gives per facet how thick the body is behind it, and the
  gap to the next wall in front. `THICKNESS_RAY` is the median of a cone of 31 rays (shape diameter);
  `THICKNESS_MEDIAL` the diameter of the largest inscribed ball. `invalid` counts the facets where too few rays
  counted, which flags an open or misoriented surface.
- `thin_regions(thickness=, gap=, h=, factor=, region=, region_area=)` groups the facets thinner than `factor * h`
  into edge-connected regions: `region(i)` is the region of facet i (0 when it is not thin), `region_area(r)` the
  area of region r, largest first. With h the cell size of the grid to come, these are the places where a wall falls
  between cells.

<<< @/examples/output/tut12-thickness.ansi{ansi}

- Two regions are thinner than h = 0.1, both on the web: the end face (14 facets) and the bore wall facing it
  (62 facets). No facet is invalid.
- The thinnest value, 0.0646, is the web of 0.06 seen through the cone: on a plate of thickness t the ray method
  gives 1.072 t (see [Thickness](/guide/advanced/thickness)).

![the remeshed bracket seen from its x = 0 end, each facet coloured by its ray thickness from teal (thin) to ember (0.5 and more), the 76 facets thinner than h = 0.1 outlined in gold on the end face over the web left by the bore (the second region, the bore wall facing it, is hidden inside the bore)](../../examples/images/tut12-thickness.svg){.fossil-figure}

## Mass properties

<<< @/examples/snippets/tut12_f_mass-mass.f90

- `mass_properties(mass=, center=, inertia=, principal_moments=, principal_axes=, density=, about=, status=)` gives
  the exact mass properties of a closed body (`STATUS_INVALID_INPUT` otherwise), as in
  [chapter 11](./11-moving-bodies): the mass, the centre of mass, the inertia tensor about it, and its principal
  moments (ascending) with their axes as columns. With the density of aluminium in SI units, the lengths are metres.

<<< @/examples/output/tut12-mass.ansi{ansi}

- The bracket weighs 2707.87 kg. Its centre of mass is at x = 1.1704, towards the boss and away from the bore (the
  base alone has it at x = 1), on the plane of symmetry y = 0.5, and at z = 0.3186, lifted by the boss.
- The principal axes are y and two axes in the xz plane tilted by 0.094 rad (5.4 degrees) from x and z: the part is
  symmetric about y = 0.5 only. The smallest moment, 346.283, is about the long axis.
- The mass of the remesh differs from that of the part by −4.94e-5, the volume change of step 3: a uniform density
  scales the volume.

![the remeshed bracket translucent, its centre of mass at (1.170, 0.500, 0.319) beside the centre of the base, and its three principal axes drawn over their radii of gyration sqrt(I / m), labelled with the moments 346, 891 and 1000](../../examples/images/tut12-mass.svg){.fossil-figure}

::: tip What you learned
`boolean` builds a watertight part from closed primitives, exactly; `set_patch` names its faces; `isotropic_remesh`
with `preserve_features` refines it towards one edge length (it does not coarsen yet); `deviation` certifies how far
the result moved; and `feature_edges`, `thickness` with `thin_regions`, and `mass_properties` measure it.
Reference: [Boolean operations](/guide/advanced/booleans), [Isotropic remeshing](/guide/advanced/isotropic-remesh),
[Surface deviation](/guide/advanced/surface-deviation), [Named patches](/guide/advanced/patches),
[Feature edges](/guide/advanced/feature-edges), [Thickness and gaps](/guide/advanced/thickness),
[Mass properties](/guide/advanced/mass-properties).
:::

Next: [13. Structured blocks and overset grids](./13-blocks-overset).
