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
two poses: fresh (solid → fluid; newborn in Bayyuk et al. 1993, emerging in Schneiders et al. 2013) and dead
(fluid → solid; vanishing, submerging).

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
aperture averaged over the step (Falcovitz et al. 1997; Murman et al. 2003, §4). These apertures are those of the
unsplit motion: Falcovitz et al. split it by direction, and Bennett et al. (2018, §2) find that directional swept
volumes do not add up to the unsplit one, so a solver that splits by direction cannot use them as they are. Three
decisions follow from the fixed grid (issue #41):

- **Nothing sweeps through a fixed face.**
  - The volume the wall sweeps through a cell is its fluid volume change,
    by Reynolds transport.
  - So `swept` = (α₁ − α₀)·|cell| comes exactly from the cut cells at the
    two poses, in 2-D and 3-D.
  - The geometric conservation law holds by construction for the geometry: volume change = swept volume (Murman et
    al. 2003, §4). A scheme's GCL also depends on how its time integrator uses them (Schneiders et al. 2016, §3.6).
  - Natarajan et al. (2022) divide by the old volume fraction with no swept term (§2.1) and lose up to 0.92 % of the
    mass at 32³, 0.31 % at 64³ (§3.6); the paper does not isolate the cause. They attribute their first order on a
    receding piston to the fresh-cell fill (§3.2.2), untested: their flux redistribution is itself first order at
    cut cells (Berger & Giuliani 2021, §1), a disagreement the llm-wiki records as open.
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
(issue #84): the residual is at most 9.69·10⁻¹⁷ in every cell, and the swept volumes add up to -3.14·10⁻¹⁷.

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
  larger of the face change between `nq` and `2 nq` nodes and the largest |residual| divided by the
  largest rigid-field speed on the grid (issue #104): a face's Gauss error enters its cells' balance
  times that speed, so the estimate is in face units and never below what the residual shows.
  `refine=.true.` (screw motion only) integrates the faces of every cell with an event inside the step
  per event interval on the cell's own wall elements, the piece map's machinery: the residual drops to
  round-off, and the faces of the other cells keep their Gauss values (no kink there).

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

In a whole time step (move, tag, cut + links, swept) the swept geometry stays the largest device stage: 37–85 % of
the step at 64³–128³. The step gains 5.4–9.3× there; the table is in
[GPU queries](/guide/advanced/device-distance#a-moving-body-time-step).

## Piece map

A cut cell can hold several fluid pieces, and over a step the wall can move through it. A solver must then know
which fluid at the end of the step came from which fluid at its start (issue #102). `swept_cells_2d(..., map=)`
answers with the **tracks** of every cell, and `swept_cells_3d(..., map=)` does the same in 3-D.

**In the literature.** Bayyuk et al. (1993) merge cells so that no combined cell changes its intersection topology
within a step, which needs the body wider than the cut cells and the step shorter than the smallest uncut cell.
Meinke et al. (2013, §3.2) handle cells cut several times or split into unconnected parts. Murman et al. (2003, §4)
and Badia et al. (2023, §2) treat a cut cell as a space-time cell. None of the sources in the llm-wiki bundle
`immersed-boundary-methods` is recorded as linking a cell's pieces through a step.

**Definition.** Take the cell's fluid over the whole step, as a region of space and time:
{(p, s) : p in the cell, p in the fluid at pose s, s in [0, 1]}. Its connected components are the tracks.
- **Pieces:** every piece of the cut cells at pose0, and every piece at pose1, lies in exactly one track.
- **Fresh and dead:** a track with no old piece is fresh; one with no new piece is dead.
- **Splits and merges:** several pieces on one side make a split or a merge.
- **Transient:** a track with no piece at either end is fluid that crosses the cell within the step.

<<< @/examples/snippets/anim_ref_moving_bodies_tracks-step.f90

![a thin plate 0.9 by 0.03, tilted 60 degrees, moving right by 0.16 per step and turning 0.04 rad over a 24 by 12 grid of spacing 0.1, each frame one step: the cells with several tracks in yellow (18 to 21 per step), along the band the plate swept between its start (grey) and its end (red); there the fluid ahead of the plate was pushed out of the cell, a dead track, while the fluid behind it grew, a continuing one; cells the wall touches with one continuing track in grey](../../examples/images/anim-ref-moving-bodies-tracks.svg){.fossil-figure}

<<< @/examples/output/anim-ref-moving-bodies-tracks.ansi{ansi}

**The end-pose overlap is not enough.** Linking a piece at t with a piece at t + dt when they overlap gets two cases
wrong, and both happen where a body is thin compared with its step:

| In one cell, over one step | Tracks | Overlap of the end poses |
|---|---|---|
| The plate pushes the fluid ahead of it out of the cell, and the fluid behind it fills the cell | a dead track and a continuing one | one merged piece: the fluid on both sides of the wall mixes |
| A thin plate crosses the whole cell | a dead track and a fresh one | "continuing" |
| A narrow gap moves more than its own width | continuing | dead + fresh |

Published schemes keep a step below one cell: Falcovitz et al. (1997, §3) assume no cell goes from fully exposed to
fully covered in one step, Bayyuk et al. (1993) bound the step by the smallest uncut cell, and Mittal et al. (2008,
§2.2) get one layer of fresh cells per step. The plate rows above and the 2-cell bench steps lie outside those
assumptions. On the bench bodies (`fossil_val_device_chain`'s ellipse, NACA 0012, sphere and cone, at displacements
of 0.5 and 2 cells per step), the overlap disagrees with the tracks in 0.02–0.14 % of the cells the wall touches. Most of these
cells are at the NACA trailing edge and on the cone's base rim.

**How it is built** (in 2-D):
- **Events.** The topology of a cell's fluid changes only when the wall's contact with the cell boundary changes:
  - a wall vertex crossing a face line;
  - a cell corner crossing a wall segment.

  These are the instants the exact swept faces already isolate. Two parts of a rigid curve keep their distance, so
  they cannot meet inside a cell.
- **Pieces.** They are evaluated at the two end poses, which must reproduce the end cuts bit for bit, and once inside
  every interval between events.
- **Links.** Pieces of consecutive evaluations that share a boundary feature (a wall segment, a fluid corner) are one
  track. A feature present on both sides of an event keeps a non-degenerate part in the cell, and the fluid next to
  it is continuous.
- **The trivial case.** A cell with no event has one track per piece, unchanged. "One piece at each end" is not
  enough: on the 2-D bench motions, 148 such cells have two tracks.

**What a track carries:**
- its old and new pieces, numbered as `cut_cells` numbers them on a moved copy, plus the inverse maps
  `old_track` / `new_track`;
- V0 and V1;
- the time-averaged fluid length of each face that belongs to it;
- its wall −Σ_f Ā_f n_f;
- `outflow`: the time-integrated outflow of the step's rigid field through each face part, the consistency term that
  keeps a uniform flow moving with the body exact (under rotation the velocity varies along a face, so the apertures
  alone do not give it);
- its residual, V1 − V0 + Σ_f outflow: the rigid-field balance of `residual`, per track.

The per-track residual is round-off on the exact paths: at most 1.3·10⁻¹⁶ in the run above. Summed over a cell's
tracks, the measures give the cell's `swept_cells_2d_t` values: copied in a one-track cell, and to round-off in the
others.

**On the Gauss path** (`exact=.false.`) a cell with several tracks still integrates per interval between its events,
because a track lives only between events and Gauss nodes spread over the whole step would miss a short-lived one. Its
tracks' residuals stay at round-off, and their sums differ from the cell's Gauss values by the quadrature error.
A one-track cell copies the cell's Gauss values, quadrature error included.

**A conservative update through the map.** `fossil_val_moving_conservation` (validation) advects a passive scalar
across 100 steps of an ellipse and a NACA 0012 at 40²–96², by translation and by the screw motion, at 0.5 and 2 cells
per step. The update has four parts:
- **Face fluxes:** one flux F_f per grid face, split to each side's tracks by their own rigid outflow `outflow`, so
  both sides exchange exactly F_f.
- **Track update:** each track takes its old pieces' mass plus its fluxes. A fresh track fills by its own inflow,
  with no interpolation.
- **Dead tracks:** a dead track hands its remainder to a neighbouring new piece.
- **New pieces:** a track's mass is shared among its new pieces by area.
- **In the literature:** Murman et al. (2003, §4) integrate an emerging cell without reconstructing it, as a fresh
  track does here. A dead track's remainder is the vanishing-cell state of Bayyuk et al. (1993), not guaranteed to
  vanish; Schneiders et al. (2016, §3.8) give it to the neighbours by volume weights.

Results:
- **Mass:** constant to 1.5·10⁻¹⁴ (relative) in every case.
- **A uniform scalar carried by the body's own rigid field:** stays uniform to 2.4·10⁻¹² of a cell on the exact
  paths. That is the per-track residual: a uniform flow moving with the body is exact, as Falcovitz et al. (1997)
  require, here shown numerically.

It is a demonstration, not a solver. Its first-order implicit upwind smears the two scalars of the leakage control
(one on each side of the body). Tracks leak less than the end-pose overlap and than a per-cell lumped update, but only
by 0.1–1 %.

**Tolerances.** The construction is exact, with two explicit thresholds:
- every event carries the round-off of its instant: the explicit translation formulas their magnitudes over their
  denominator (times 64 eps), the screw events the width `isolate_zeros` localises them to (issue #104). Two events
  closer than the sum of their widths are one instant in double precision; the state between two farther apart is
  evaluated at its middle, however short (`narrow_cells` counts the cells with one shorter than 1e-10): a plate
  2e-11 thick crossing a cell reads dead + fresh. In the bench every pair closer than 1e-10 is round-off (1e-17 to
  2.5e-13 apart on translations, dyadic widths of uncertified zeros on rotations): no bench cell changed;
- a piece with no known feature is dropped only if its face lengths are below 1e-9 of the cell side.

- **Exact tangential contact under rotation.** A wall feature can touch a cell feature tangentially and leave
  without crossing, a double root of an event function. When that root lies within the round-off of the event
  function, it cannot be certified. Its cell is then marked `TRACK_UNTRUSTED` and counted in `contact_cells`, never
  linked by guess. Off by more than about 1e-12, the two roots are certified and the tracks are right. Translation
  cannot produce it, because its event functions are linear.

A cell whose evaluations cannot be matched is also flagged `TRACK_UNTRUSTED`: one track with all its pieces and the
cell's values, so a conservative update stays exact there. Apart from the constructed exact-contact cases, there were
none in the tests, the benches or the conservation runs.

**In 3-D,** `swept_cells_3d(..., map=)` returns a `piece_map_3d_t` with the same layout: six faces, a wall of three
components, `outflow(6, nt)`.
- **Events:** a mesh vertex crossing one of the cell's face planes, a mesh edge crossing one of its edges, a cell
  corner crossing a facet. They come from the per-plane breakpoints of the exact swept faces.
- **Features:** facet ids (a triangle clipped to a box is convex and bounds one piece) and fluid corners.
- **Results:** on the sphere and the cone, under translation, exact rotation and Gauss rotation, there are 0
  untrusted cells, and the per-track residuals stay at most 2.7·10⁻¹⁵ on the exact paths.

**Cost:** in 3-D, 3–13 % on top of a translating `swept_cells_3d` at 32³–64³ and 49–281 % on top of a Gauss rotation,
where the swept call itself is cheap (0.15–1.7 s). In 2-D, 1–22 % on top of `swept_cells_2d` at 256²–1024², and 0 to 281 cells with several tracks per step (the
NACA 0012 under a Gauss rotation at 1024²). The map is not combined with `poses=` (it needs the screw motion).

### The map on the device

`swept_cells_2d_device(..., map=)` and `swept_cells_3d_device(..., map=)` (issue #105) return the map of
`swept_cells_2d(..., map=)` / `swept_cells_3d(..., map=)`, bit for bit. Every other output is that of the call without
`map`.

- **On the device:** the per-item passes. The one-cell evaluations, one item per (cell, evaluation), each with its own
  workspace sized by a count pass; the tracks of every cell; the measure rounds, one item per (cell, interval, node),
  summed in node order; the assembly by count, scan and fill. In 3-D the posed facets, the plan's facets and the
  evaluations stay resident for the whole build.
- **On the host:** the plan (candidate cells, events, the pose of every evaluation and node, the Gauss nodes), the
  corner flags (2-D) or corner states (3-D), the scans and the convergence decisions of each round. The device `cos`,
  `sin` and `atan2` differ from the host's by 1 ulp in 17 %, 14 % and 1.8 % of random arguments (RTX 4070, nvfortran
  26.9, `-gpu=nofma`), so events and poses computed on the device would not be the host's bits. The device reads the
  host's pose table.
- **End poses:** with `map` they are cut in full, by `cut_cells_3d_device` without `alpha_only` in 3-D and on the host
  in 2-D: the map needs their pieces.
- **`chunk_bytes`** (default 2²⁸) caps the workspace of one launch set of items. No output depends on it.

Measured by `fossil_test_distance_bench` on an RTX 4070 shared with other jobs, one call, best of 2, so the ratios are
indicative. "Map" is the Gauss call
(nq = 4) with `map` less the call without; "refine machinery" is the call with `refine=.true.` less the Gauss call:

| body | cells | map, CPU [s] | map, device [s] | speed-up | refine machinery, CPU [s] | device [s] | speed-up |
|---|---|---|---|---|---|---|---|
| sphere (marching cubes) | 32³ | 0.373 | 0.396 | 0.94× | 6.00 | 0.958 | 6.3× |
| sphere | 64³ | 0.635 | 0.116 | 5.5× | 18.7 | 1.75 | 10.7× |
| cone | 32³ | 0.558 | 0.628 | 0.89× | 2.46 | 0.598 | 4.1× |
| cone | 64³ | 3.93 | 3.58 | 1.10× | 17.0 | 4.66 | 3.7× |

- **With `refine`** the per-item passes are most of the work, and they gain 3.7–10.7×.
- **On the 3-D Gauss path** the map gains 0.89–1.10× except on the sphere at 64³ (1.2–2.1× in the device chain,
  over 3 warm steps): the host plan dominates, above all the events of the rotation (`isolate_zeros`), which stay on
  the host for the trig reason above.
- **In 2-D** the device map gains nothing: few cells are evaluated (about 3,000 of a million at 1024²), so there is
  little device work to gain from, and the copies dominate. In `fossil_test_curve_bench` (ellipse) the Gauss map runs
  at 0.74–0.77× the host's at 256² and 1024², the translation map at 0.48× at 1024²; in the device chain 0.8–1.2×.
  The map is timed as the difference of two calls (with `map` and without), so where it is a small part of the call it
  falls inside the timing noise: the bench then prints a negative map time (translation 256², Gauss + refine 1024²).

The chain of a whole time step with the map is in
[GPU queries](/guide/advanced/device-distance#a-moving-body-time-step).

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
| piece map on the device vs the host: every field of `piece_map_3d_t` / `piece_map_2d_t`; sphere and cone at 32³ (exact translation, exact rotation, Gauss, Gauss + `refine`); ellipse at 256² (translation, Gauss, Gauss + `refine`), NACA 0012 rotating at 64²; launch sets of 4 KB (`chunk_bytes`) | 0 mismatches (`piece map 3-D` / `2-D device vs serial mismatches`, `fossil_test_piece_map_3d_device_twin`, `fossil_test_piece_map_2d_device_twin`; the map stage of `fossil_val_device_chain`) |

## Limitations

- **Swept geometry: the Gauss path at large steps.** With `exact=.false.` the face apertures use nq Gauss nodes over
  the whole step. The face integrands kink, and jump, at the events (a vertex crossing a face line, a corner crossing
  a segment, a face leaving the body), and global nodes resolve them at about first order.
  - Measured (issue #104; ellipse and NACA 0012 at 256², sphere and cone at 32³; 0.5, 1 and 2 cells per step): the
    largest per-cell residual is 7·10⁻⁴ to 0.45 of a cell at nq = 4, and still 2·10⁻⁴ to 0.15 at nq = 8.
  - `error_estimate` includes the residual since issue #104; the face change alone understated it a hundredfold.
  - `refine=.true.` makes the residual round-off: 1.1–2.2 times the Gauss call in 2-D (the exact rotation path takes
    15–45 times), 2.7–22 times in 3-D, about the exact path's cost there. It is off by default, and the default is
    an open decision. The exact path stays the reference.
  - The piece map's multi-track cells do not inherit the error: they integrate per event interval on every path.
- **Spurious force oscillations are not addressed.** Exact swept volumes remove the fresh/dead volume source that
  Seo & Mittal (2011, §2) call their primary cause. Schneiders et al. (2013, §4-§5) show oscillations in a
  conservative scheme whose discrete operator switches, a merge target turning on or off; the llm-wiki records the
  disagreement as open.
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
    (cube: residual 1.1e-2 → 1.6e-3 from nq = 2 to 8). The exact paths
    (the default) and `refine` do not have this error.
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
