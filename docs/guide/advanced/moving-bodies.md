---
title: Moving bodies
---

# Moving bodies

## What it does

A rigidly moving body needs its geometry at every new pose. Rebuilding the
tree every time step is the cost to avoid. FOSSIL offers two mechanisms, and
the kind of query decides which one to use:

| Mechanism | Use it for | Cost | Agreement with a mesh moved and rebuilt |
|---|---|---|---|
| **Body-frame queries**: pass `pose=` | distance, IB queries, inside tests (CPU and device) | one 3 × 3 product per point | distances ≤ 1e-16 (L + \|t\|), signs equal away from the surface |
| **Moved copy**: `move(pose)` with a tree refit | grid-aligned queries: box tagging, line crossings, cut cells | O(N) (a refit, not a rebuild) | bit for bit (feet at edge and vertex ties to 4 ulps) |

`state_changes` flags the points (cell centres) that change state between
two poses: fresh (solid → fluid) and dead (fluid → solid).

![a NACA 0012 section seen from above at the start of a step (dark brown) and at its end (yellow), pitched 3 degrees nose up and moved 0.02 forward and 0.01 up, inside the outline of a grid of cell centres of spacing 0.02: only the 292 centres near the wall that the band filter sends to a second query are drawn, grey; among them 59 fresh centres (teal) along the lower aft surface the section leaves and 54 dead centres (ember) along the upper surface it covers](../../examples/images/ref-moving-bodies-states.svg){.fossil-figure}

<<< @/examples/snippets/fig_ref_moving_bodies_states-states.f90

<<< @/examples/output/figref-moving-bodies-states.ansi{ansi}

Of the 1452 centres, 1160 are decided by one query: their distance at the old pose exceeds how far they move in
the body frame. The other 292 are queried again.

`swept_cells_2d` / `swept_cells_3d` give the swept geometry a conservative
moving cut-cell update needs over a step (see [Swept geometry](#swept-geometry)).

## Example

The posed queries, `state_changes` and the 2-D swept area of a dragon section are in
[tutorial chapter 11](/manual/tutorial/11-moving-bodies), and a posed query of a rotated cube in the
[cookbook](/manual/cookbook#query-a-rotated-body-without-rebuilding-it). This program checks the pose algebra, compares
the three ways to query a moved body, and runs the swept geometry of one step of the unit cube on each of its paths.

<<< @/examples/snippets/ref_moving_bodies.f90

<<< @/examples/output/ref_moving_bodies.ansi{ansi}

- **Pose algebra.** `pose_compose(pose_inverse(p1), p1)` is the identity to 1.11·10⁻¹⁶, and `to_body` undoes
  `to_world` to 5.55·10⁻¹⁷.
- **Screw motion.** At s = 1/2 the rotation, squared, is R₁ to 1.11·10⁻¹⁶, and the translation is exactly t₁/2: the
  angle and the translation are linear in s. The motion is not a power of a pose, so `pose_compose(half, half)` is
  not p1 when the translation has a component off the axis. At s = 1 `pose_interpolate` returns p1 bit for bit.
- **Three ways to query.** The posed query (`pose=`, tree as built) and the moved copy (`move`, tree refitted)
  agree to 1.11·10⁻¹⁶ on 27 points, 4 of them inside. The refitted copy and the same copy re-analyzed from scratch
  (`analyze`, tree rebuilt) agree bit for bit.
- **Translation.** The cube faces cross the grid planes x = 0.03, 1.03 and y = 0.03, 1.03 during the step, and the
  fluid areas of those planes jump at the crossings. The exact path places its breakpoints there: 14 evaluations on
  the busiest plane, residual 2.52·10⁻¹⁷. The swept volumes add up to -5.90·10⁻¹⁷: a rigid motion keeps the volume.
- **Gauss, the negative control.** `exact=.false.` integrates the same step with 4 Gauss nodes across the jumps:
  residual 1.10·10⁻³; `error_estimate`, the largest face change from 4 to 8 nodes, is 1.13·10⁻².
- **Rotation.** The same shift with a 0.2 rad turn about z takes the exact rotation path: 508 evaluations on the
  busiest plane, residual 4.34·10⁻¹⁷. 12 cells are cleared (no fluid at the start, fluid at the end) and 9 covered.
  Summed over the block, which holds the body at both poses, the projected walls cancel to 8.327·10⁻¹⁷.

The 2-D swept geometry over several steps: a square of side 0.5 rolls along a grid of spacing 0.1, turning 0.15 rad
and advancing 0.14 per step while it bobs up and down, one `swept_cells_2d` call per step between consecutive poses.

<<< @/examples/snippets/anim_ref_moving_bodies_swept-step.f90

![a square of side 0.5 rolling to the right over a 24 by 12 grid of spacing 0.1 in 10 steps, each frame one step: the cells shaded by their swept area over the cell area, ember behind the square where fluid is gained, teal ahead of it where fluid is lost, the square at the start of the step in dark brown and at its end in yellow; every step takes the exact rotation path, with a largest residual of 1.04e-16](../../examples/images/anim-ref-moving-bodies-swept.svg){.fossil-figure}

<<< @/examples/output/anim-ref-moving-bodies-swept.ansi{ansi}

Every step takes the exact rotation path (`sw%exact`); the swept areas of a step add up to at most 3.56·10⁻¹⁶ and the
largest uniform-flow residual over the 10 steps is 1.04·10⁻¹⁶.

## Poses

- **`rigid_pose_t`** holds `r(3, 3)` and `t(3)`, with x_world = R x_body + t.
  - Constructors: `pose_identity`, `pose_from_axis_angle(axis, angle, translation)`,
    `pose_from_quaternion(q, translation)`.
  - Algebra: `pose_compose(a, b)` (apply b, then a) and `pose_inverse`.
  - Maps: `to_body`, `to_world`, `rotate`, which are device-callable.
- **`rigid_pose_2d_t`** (`pose_2d(angle, translation)`) is the in-plane pose
  of a curve.
- **`pose_interpolate(p0, p1, s)`** (2-D and 3-D) is the screw motion between
  two poses: the translation and the rotation angle, about the fixed axis of
  R₁R₀ᵀ, are linear in s ∈ [0, 1]. The end poses are returned exactly.

The pose argument is on:
- `distance_many`, `query_many`, `is_point_inside` and their device twins
  `distance_many_device`, `query_many_device`;
- for surfaces and curves alike.

Closest points, second foot points and normals come back in the world frame.

![the unit cube drawn by its edges at six poses of the screw motion from the identity to a turn of 120 degrees about z with translation (1.8, 0.6, 1.2), at s = 0, 1/5, ..., 1, grey up to the yellow end pose; the body-frame origin moves on a straight teal line, one corner on an ember arc; a teal cube drawn at pose_compose(half, half), with the end rotation but displaced from the end cube by 0.949](../../examples/images/ref-moving-bodies-screw.svg){.fossil-figure}

<<< @/examples/output/figref-moving-bodies-screw.ansi{ansi}

The body-frame origin moves on the straight line s t₁ while the body turns about it at a constant rate. Composing the
half-way pose with itself gives the end rotation to 1.11·10⁻¹⁶ but misses the end position by 0.949: applying the
half pose twice turns the first half-step's translation too, which the screw motion does not.

## How it works

- **Body-frame queries.** The points are mapped by x̃ = Rᵀ (x − t); the
  results are mapped back. Distances are invariant. The device twins map on
  the host with the CPU code, so the posed device results equal the posed
  CPU results bit for bit.

  ![the unit cube at rest drawn by its edges and, translucent, at a pose turned 35 degrees about z and lifted by 0.25; a ring of 40 world points around the posed cube, each joined by a yellow segment to its closest point on the posed cube, teal where the point is outside, ember where a corner of the cube reaches past the ring](../../examples/images/ref-moving-bodies-posed.svg){.fossil-figure}

  The tree of the figure is the one built at rest; `query_many(pose=)` returns the closest points on the posed cube,
  in the world frame. 28 points are outside, 12 inside, with signed distances from -0.0684 to 0.0999.
- **Move and refit.**
  - The vertices move, and everything derived from coordinates is
    recomputed exactly as `analyze` does.
  - The SAH BVH keeps its topology: node boxes are recomputed bottom-up as
    unions of facet boxes, and the payloads are rebuilt. Curves refit their
    segment tree the same way.
  - Topology, shells and bodies are kept. An octree is rebuilt.
  - The facets equal those of the mesh moved and analyzed from scratch, so
    the queries agree with a rebuild.
  - Exception: a foot on an edge or vertex shared by tied facets is
    parametrized differently by each; which facet is reported depends on the
    tree, so its bits can differ by a few ulps.
- **Fresh and dead points.**
  - A point at the two body-frame positions b₀ and b₁ can change state only
    if the segment b₀b₁ reaches the surface.
  - So when |φ(b₀)| > |b₀ − b₁| the second query is skipped. That filter is
    exact.
  - States come from the signed distance with the chosen sign algorithm.

## Swept geometry

A conservative cut-cell update over a step with a moving wall needs, per
cell, the old and new fluid volumes. Per grid face it needs the fluid
aperture averaged over the step (Falcovitz et al. 1997). Three decisions
follow from the fixed grid (issue #41):

- **Nothing sweeps through a fixed face.**
  - The volume the wall sweeps through a cell is its fluid volume change,
    by Reynolds transport.
  - So `swept` = (α₁ − α₀)·|cell| comes exactly from the cut cells at the
    two poses, in 2-D and 3-D.
  - The geometric conservation law holds by construction.
  - Natarajan et al. (2022) omit this term and lose order.
- **The only time integral is the face aperture.**
  - `face_u` / `face_v` (2-D) and `face_x` / `face_y` / `face_z` (3-D)
    are the fluid length or area of every grid face, averaged over the
    step.
  - `wall` = −Σ_f Ā_f n_f is the time-averaged projected wall area of
    the cell, which carries the wall pressure and work. It uses the same
    quantities as the face fluxes.
- **Motion.** By default the motion is `pose_interpolate`. A caller that
  knows the true motion, such as a 6-DOF integrator, can instead pass the
  poses at the Gauss nodes (`gauss_legendre_01`).

![a 12 x 12 grid seen from above, the middle layer of a 12 x 12 x 3 block, each cell shaded by the volume the unit cube sweeps through it over one step of a screw motion (0.35 rad about a tilted axis plus a shift), divided by the cell volume: ember where the cell gains fluid on the trailing side, teal where it loses fluid on the leading side, cream elsewhere; the section of the cube at mid-layer at the start in dark brown and at the end in yellow](../../examples/images/ref-moving-bodies-swept.svg){.fossil-figure}

In the figure the step takes the exact rotation path, with the cube's faces x, y = ±0.5 starting on grid planes
(issue #84): the residual is at most 9.69·10⁻¹⁷ in every cell, and the swept volumes add up to -3.84·10⁻¹⁷.

```fortran
call swept_cells_2d(curve, pose0, pose1, x, y, sw, nq, poses, exact, status, check_input)
call swept_cells_3d(surface, pose0, pose1, x, y, z, sw, nq, poses, exact, status, check_input)
```

**Integration.**
- **2-D pure translation: exact.** Both poses share their rotation.
  - Each face is an axis-aligned segment moving against the curve, so
    its fluid length is piecewise linear in time.
  - The breakpoints are the instants when a wall vertex crosses the
    face's line, or when an end node of the face crosses a wall segment.
  - Two-point Gauss on every piece is exact.
- **3-D pure translation: exact.** The same idea, one grid plane at a time.
  - In the frame of the body, a grid plane moves by −s d with its grid
    lines. Every vertex of a face's fluid polygon moves linearly in time
    while the combinatorics hold, so the area is quadratic in time.
  - The breakpoints of a plane are the instants when a mesh vertex
    crosses it, a mesh edge crosses one of its grid lines, or one of its
    grid nodes crosses a triangle.
  - Two-point Gauss on every piece is exact. One evaluation (a section
    of the shifted plane and the 2-D cut cells on the shifted grid)
    gives every face of the plane.
- `exact=.false.` forces the Gauss path for a translation (both
  dimensions); `sw%exact` reports which path ran and `sw%evaluations`
  how many aperture evaluations it took (per face in 2-D, per plane in
  3-D).
- **2-D rotation: exact** (issue #69, see *Exact rotation* below). The screw motion with a rotation takes
  certified breakpoints per face and Gauss–Legendre doubled to round-off on every piece; the residual
  falls to round-off (1.7e-16 to 6.3e-16 on the cases below, where 4-node Gauss left 1.7e-4 to
  1.4e-2). `sw%evaluations` is the largest number of evaluations per face (92 for the square, 1838
  for the 200-vertex NACA 0012).
- **3-D rotation: exact** (issue #69). The same scheme one grid plane at a time: a plane's breakpoints
  are a mesh vertex on it, a mesh edge meeting one of its grid lines, one of its grid nodes on a
  triangle; one evaluation slices a moved copy and gives every face of the plane. Residual 1.4e-16 to
  1.8e-15 on the cases below, where 4-node Gauss left 1.0e-3 to 9.3e-3.
- **Otherwise (explicit node poses, or `exact=.false.`): Gauss–Legendre at `nq` nodes**
  (default 4), on copies moved to the node poses. With the screw motion, `error_estimate` is the
  largest face change between `nq` and `2 nq` nodes.

**Exact rotation (design, issue #69).** The Gauss path converges algebraically under rotation (rate 2–3,
see the validation table), because a face's fluid measure has kinks inside the step. The exact path
removes them, as the translation path does, and is built as follows.

- *Motion.* The screw motion of `pose_interpolate` has the angle linear in s ∈ [0, 1] about a fixed
  axis (the plane normal in 2-D) and the translation linear in s: a body point q is at
  X(s) = R_k(sΔθ) R₀ q + t₀ + s Δt. It is not a rotation about a fixed pivot, so the events below are
  transcendental in general.
- *Events.* A face's fluid measure is smooth between the instants when
  - a wall vertex crosses the face's line (2-D) or plane (3-D);
  - an end node of the face crosses a wall segment (2-D), or a grid node of the plane crosses a
    triangle (3-D);
  - a mesh edge crosses one of the plane's grid lines (3-D).

  Each is a zero of an orientation of moving and fixed points, and with X(s) as above every one has
  the form

  g(s) = (p₀ + p₁ s) cos φ + (q₀ + q₁ s) sin φ + k₀ + k₁ s,   φ = φ₀ + s Δθ,

  with coefficients from the input points (the rotation acts on differences of body points, the
  translation enters linearly), plus a check that the event lies on the face (the vertex inside the
  span, the node's foot inside the segment or triangle).
- *Certified isolation.* |g''(s)| ≤ M = Δθ² max|(p, q)| + 2|Δθ| |(p₁, q₁)|, the maximum over [0, 1].
  On an interval of length h whose end values have the same sign and exceed M h²/8 in magnitude, g has
  no zero (g lies above its chord minus M h²/8). Otherwise the interval is halved; a sign change is
  refined by bisection to one ulp of s. No pair of zeros is missed. An interval narrower than 1e-10
  without a certificate (a tangency, a zero at an end of the step) gives one breakpoint; an interval on
  which |g| stays at round-off throughout (a zero zone: a mesh face starting on a grid plane) gets
  breakpoints at its two ends only; breakpoints within 1e-10 of each other merge. A kink that close to
  a breakpoint changes its piece's Gauss value by ~slope × 1e-20. If the 4096-interval budget of one
  event runs out, `sw%exact` is false.
- *Integration.* On each piece between breakpoints the fluid measure is analytic, so Gauss–Legendre
  converges exponentially. The node count on a piece doubles from 2 until two successive values agree
  to round-off of the face measure (pieces shorter than 1e-6 to 64 ε of the measure times 1e-6: slices
  within 1e-10 of a mesh face carry ~1e-7 relative noise); `sw%evaluations` reports the largest count
  (per face in 2-D, per plane in 3-D). At 64 nodes the target becomes the round-off of one evaluation:
  64 ε of the coordinate magnitude (2-D) or of that times the face side (3-D), no longer times the
  piece length (issues #84, #86). Near an instant when a wall side lies along a face's line, or a body
  face on a grid plane (a side or face starting there, a corner leaving a line tangentially), the
  crossings of a moved copy carry round-off amplified by 1/angle, and successive Gauss values stop
  agreeing at a floor that does not shrink with the piece. A piece that misses this target too marks
  the result `sw%exact = .false.`, never silently. The flag still catches missing breakpoints: with
  every breakpoint dropped, the dragon sections and the cube of issue #84 all come out not exact
  (residuals 1e-6 to 2e-4). Each evaluation moves a
  copy of the body to the node's pose and measures the face exactly as the cut cells do, so apertures
  and fractions share one convention.
- *What changes.* The swept volumes are exact already (cut-cell volumes at the two poses), and `wall`
  closes the uniform flow by construction, so a conservative update is consistent either way
  (Falcovitz's requirement). The exact path makes the time-averaged apertures exact under rotation too,
  and with them the rigid-field `residual`, which falls to round-off. Explicit node poses (`poses=`)
  keep the Gauss path; `exact=.false.` forces it.

**Check.** `residual` is the volume balance of the body's rigid velocity
field: swept volume plus the time integral of the field's outflow through
the fluid part of every face. It vanishes when a uniform flow moving with
the body is reproduced exactly, Falcovitz's design requirement. It is not
computed (0) when the caller passes the node poses. The residual also
checks the input: the exact paths assume a wall that does not cross itself
(see the limitation below).

`cleared` and `covered` flag the cells with no fluid at one pose and
fluid at the other.

### On the device

`swept_cells_3d_device(surface, pose0, pose1, x, y, z, sw, ...)` and `swept_cells_2d_device(curve, pose0, pose1, x, y,
sw, ...)` (issue #100) take the arguments of `swept_cells_3d` / `swept_cells_2d` and return the same `sw`, bit for
bit, on every path: exact translation, exact rotation, Gauss and explicit poses.

**In 3-D:**
- The end poses go through `cut_cells_3d_device(..., alpha_only=.true.)`.
- The faces are batched: the host builds each (plane, Gauss node) as `swept_cells_3d` does, with the breakpoints, the
  events and `isolate_zeros` on the host (they are at most 1 % of a call). It slices every node through a plane sweep,
  and one multi-block launch per chunk cuts all the sections on their own grids.
- Each face then adds over its nodes in node order on one device thread: the same additions in the same order.
- The optional `chunk_cells` caps the device memory of a launch set without changing a bit.

**In 2-D** only the end poses move to the device. The faces, measured at 94–99.9 % of a call, stay on the host, so the
twin gains nothing measurable today.

Measured on an RTX 4070, best of 5 after a warm-up:

| body | cells | path | CPU [s] | device [s] |
|---|---|---|---|---|
| sphere (marching cubes) | 32³ | translation | 10.4 | 6.1 |
| sphere | 32³ | rotation | 15.0 | 8.3 |
| sphere | 64³ | translation | 56.8 | 19.8 |
| sphere | 64³ | Gauss, nq = 4 | 2.37 | 0.59 |
| cone | 32³ | rotation | 3.95 | 1.38 |
| cone | 64³ | translation | 17.5 | 3.8 |
| cone | 64³ | Gauss, nq = 4 | 1.67 | 0.38 |

That is 1.7–4.7×. The exact paths are dominated by the number of (plane, node) sections: 650–15 400 per plane at 64³,
and the exact rotation costs 250–600 s per host call there. The Gauss path is two to three orders of magnitude cheaper
than the exact rotation at 64³.

## Results

| Case | Result |
|---|---|
| body-frame vs moved and rebuilt: cube (also translated by 108), cone, L-prism, bunny, NACA 0012, annulus | distances ≤ 1e-16 (L + \|t\|), closest points ≤ 6.5e-15 (L + \|t\|), signs and inside flags equal |
| posed device vs posed CPU | every output bit for bit (device path run on the host) |
| move (refit) vs moved and rebuilt: cube, cone, L-prism, bunny, octree cube, NACA 0012, annulus | geometry, distances, inside flags, box codes, line crossings, cut cells bit for bit |
| refit time / rebuild time | 0.03 (bunny) to 0.4 (tiny meshes) |
| fresh / dead: cube translated, rotated about an off-centre axis, general poses; square, pitching NACA 0012 | 0 errors against the analytic change of membership; the band skips 9–60 % of the second queries |
| swept, 2-D translation: square (every face vs analytic), non-convex heptagon, heptagon at a rotated pose | faces to 8.3e-17; uniform-flow residual ≤ 1.7e-16 in every cell; walls sum to ≤ 1.4e-16 |
| swept, 2-D negative control: Gauss with 4 nodes, no breakpoints | residual 9.0e-3 |
| swept, 2-D rotation, Gauss (`exact=.false.`): off-centre square, pitching NACA 0012; nq = 2, 4, 8, 16 | residual 5.6e-3 → 1.6e-4 (rate ≈ 2.3) and 1.5e-3 → 1.6e-6 (rate ≈ 3.4) |
| swept, 3-D rotation, exact (issue #69): cube 0.3 rad about (1, 1, 1) + shift (faces starting on grid planes), L-prism 0.5 rad about a tilted axis + shift, 32-face cone pitching 0.4 rad | residual 1.4e-16, 1.8e-15, 5.0e-16 (6³ validation: 4.3e-17, 8.3e-16, 5.0e-16); walls sum ≤ 1.8e-15; Gauss nq = 8 ≥ 4× closer than nq = 2 |
| swept, 2-D rotation, exact (issue #69): square 0.4 rad off-centre + shift, square 1.2 rad + shift, NACA 0012 pitching 0.15 rad | residual 1.7e-16, 1.8e-16, 6.3e-16; walls sum ≤ 1.9e-16; Gauss nq = 16 ≥ 4× closer to the exact faces than nq = 4 |
| swept, rotation with the wall starting on grid lines / planes (issues #84, #86; `fossil_test_swept_rotation_onplane`): square [−0.5, 0.5]² on a 0.5 grid, 0.1 rad + shift; centred unit cube 0.35 rad about (0.3, 0.2, 1) + shift, faces y = ±0.5 on grid planes | `exact` (both were not before); residual 1.3e-16, 4.4e-16; walls sum ≤ 5.6e-17 |
| swept, 2-D rotation of dragon sections, z = −0.35, −0.20, −0.10, −0.05, π/18 on a 0.025 grid (`fossil_val_swept_rotation_2d`) | all `exact` (before: one or two of the four not, which ones depending on the round-off of the placed mesh; turned in 10-degree steps, most steps not); residual ≤ 1.2e-17 |
| swept, 3-D translation: unit cube (every face vs analytic) | faces to 5.6e-17; uniform-flow residual 2.5e-17; walls sum to ≤ 8.3e-17 |
| swept, 3-D translation: cube and L-prism at a pose rotated about (1, 1, 1), 32-face cone | uniform-flow residual 7.6e-17, 3.0e-16, 1.9e-16 |
| swept, 3-D translation without one breakpoint class (vertex / edge-line / node-triangle) | residual up to 5.7e-4 / 1.4e-3 / 1.3e-5: each class is needed |
| swept, 3-D negative control: unit cube translating, Gauss, nq = 2, 4, 8, 16 | faces vs analytic 2.5e-2 → 2.3e-3; residual 1.3e-2 → 7.1e-4 |
| swept, 3-D: cube rotating about (1, 1, 1) while translating; nq = 2, 4, 8 | residual 1.1e-2 → 1.6e-3 |
| swept: explicit node poses = screw poses | face apertures bit for bit |
| swept on the device vs the host: every output of `swept_cells_3d_t` / `swept_cells_2d_t`, exact translation, exact rotation, Gauss, explicit poses; chunk boundaries (`chunk_cells`) | 0 mismatches (`swept cells 3-D` / `2-D device vs serial mismatches`, `fossil_test_swept_3d_device_twin`, `fossil_test_swept_2d_device_twin`) |

## Limitations

- **Swept geometry: a self-intersecting wall breaks the balance.** The swept
  volumes come from the cut cells (Green's theorem: a region counts by its
  winding number) and the face apertures from line crossings (parity). Where
  the wall crosses itself the two count the knot differently, so the cells
  around the knot get a residual that is not round-off, even on the exact
  translation path. Measured on a section of the (non-manifold) dragon that
  crosses itself once at a fold: residual 5.3e-6 in the cell the knot enters;
  with the knot removed, 7.8e-18 over the same 25 shifts. Pass
  `check_input=.true.` to `swept_cells_2d` / `swept_cells_3d`: the input is
  checked first (one `self_intersections`, or `find_self_intersections` in
  3-D, which costs more), `sw%self_intersecting` is set and a warning goes to
  stderr; the result is computed anyway.
- **Swept geometry: rotations converge slowly across jumps.**
  - A face's fluid length (area) jumps when the face's line (plane)
    leaves the body. Gauss across a jump is O(1/nq).
  - Rotation in 2-D gives rates of about 2–3; in 3-D about first order
    (cube: residual 1.1e-2 → 1.6e-3 from nq = 2 to 8).
  - Only translations are exact: under a rotation the polygon vertices
    no longer move linearly, so there is no finite set of pieces on which
    a fixed rule is exact.
- **Swept geometry: the cost.** The 2-D exact path costs O((faces +
  nodes) × segments). The 3-D exact path costs, per grid plane, a scan of
  the triangles in the slab the plane sweeps, plus 2 (breakpoints + 1)
  sections and 2-D cut-cell passes; the breakpoints grow with the mesh
  vertices and edges in that slab, so a fine mesh and a long step make
  it expensive (`sw%evaluations` reports it). The Gauss paths cost one
  moved copy and one face pass per node, plus as many again for the
  error estimate.
- **Rigid bodies only.** Deforming bodies, and refitting under deformation,
  are out of scope.
- **Points near the surface.** `state_changes` decides states with the
  point-query sign, so a point within round-off of the wall can differ from
  the cut-cell parity.
- **GPU not validated yet.** The posed device paths run on the host in
  `tests-gnu`; the GPU benches print `posed device vs serial mismatches` for
  the manual gate.
