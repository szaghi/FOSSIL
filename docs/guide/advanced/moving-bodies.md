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

`swept_cells_2d` / `swept_cells_3d` give the swept geometry a conservative
moving cut-cell update needs over a step (see [Swept geometry](#swept-geometry)).

## Usage

```fortran
program moving
use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, pose_identity
use penf,   only : R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: body, moved
type(rigid_pose_t)       :: p0, p1
type(vector_R8P)         :: centres(2)
real(R8P)                :: d(2)
logical                  :: fresh(2), dead(2)

call body%load_from_file(file_name='body.stl', guess_format=.true.)
p0 = pose_identity()
p1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0.1_R8P, translation=[0.05_R8P, 0._R8P, 0._R8P])
centres = [vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(1.2_R8P, 0.3_R8P, 0.4_R8P)]

! body-frame query: the points are in the world frame, the tree stays as built
call body%distance_many(points=centres, distances=d, is_signed=.true., pose=p1)
! cells that change state between the two poses
call body%state_changes(pose0=p0, pose1=p1, points=centres, fresh=fresh, dead=dead)
! a moved copy for grid-aligned queries (box tagging, line crossings, cut cells)
moved = body
call moved%move(p1)
endprogram moving
```

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

## How it works

- **Body-frame queries.** The points are mapped by x̃ = Rᵀ (x − t); the
  results are mapped back. Distances are invariant. The device twins map on
  the host with the CPU code, so the posed device results equal the posed
  CPU results bit for bit.
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

```fortran
call swept_cells_2d(curve, pose0, pose1, x, y, sw, nq, poses, exact, status)
call swept_cells_3d(surface, pose0, pose1, x, y, z, sw, nq, poses, status)
```

**Integration.**
- **2-D pure translation: exact.** Both poses share their rotation.
  - Each face is an axis-aligned segment moving against the curve, so
    its fluid length is piecewise linear in time.
  - The breakpoints are the instants when a wall vertex crosses the
    face's line, or when an end node of the face crosses a wall segment.
  - Two-point Gauss on every piece is exact.
- **Otherwise: Gauss–Legendre at `nq` nodes** (default 4), on copies moved
  to the node poses. With the screw motion, `error_estimate` is the
  largest face change between `nq` and `2 nq` nodes.

**Check.** `residual` is the volume balance of the body's rigid velocity
field: swept volume plus the time integral of the field's outflow through
the fluid part of every face. It vanishes when a uniform flow moving with
the body is reproduced exactly, Falcovitz's design requirement. It is not
computed (0) when the caller passes the node poses.

`cleared` and `covered` flag the cells with no fluid at one pose and
fluid at the other.

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
| swept, 2-D rotation: off-centre square, pitching NACA 0012; nq = 2, 4, 8, 16 | residual 5.6e-3 → 1.6e-4 (rate ≈ 2.3) and 1.5e-3 → 1.6e-6 (rate ≈ 3.4) |
| swept, 3-D: unit cube translating; nq = 2, 4, 8, 16 | faces vs analytic 2.5e-2 → 2.3e-3; residual 1.3e-2 → 7.1e-4 |
| swept, 3-D: cube rotating about (1, 1, 1) while translating; nq = 2, 4, 8 | residual 1.1e-2 → 1.6e-3 |
| swept: explicit node poses = screw poses | face apertures bit for bit |

## Limitations

- **Swept geometry: Gauss converges slowly across jumps.**
  - A face's fluid length (area) jumps when the face's line (plane)
    leaves the body. Gauss across a jump is O(1/nq).
  - Rotation in 2-D gives rates of about 2–3.
  - 3-D translation is about first order (faces 2.5e-2 → 2.3e-3 from
    nq = 2 to 16).
  - Only the 2-D translation is exact. A per-plane breakpoint path for 3-D
    translations is a possible extension.
- **Swept geometry: the cost.** The 2-D exact path costs O((faces +
  nodes) × segments). The Gauss paths cost one moved copy and one
  face pass per node, plus as many again for the error estimate.
- **Rigid bodies only.** Deforming bodies, and refitting under deformation,
  are out of scope.
- **Points near the surface.** `state_changes` decides states with the
  point-query sign, so a point within round-off of the wall can differ from
  the cut-cell parity.
- **GPU not validated yet.** The posed device paths run on the host in
  `tests-gnu`; the GPU benches print `posed device vs serial mismatches` for
  the manual gate.
