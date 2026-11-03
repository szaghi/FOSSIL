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

## Results

| Case | Result |
|---|---|
| body-frame vs moved and rebuilt: cube (also translated by 108), cone, L-prism, bunny, NACA 0012, annulus | distances ≤ 1e-16 (L + \|t\|), closest points ≤ 6.5e-15 (L + \|t\|), signs and inside flags equal |
| posed device vs posed CPU | every output bit for bit (device path run on the host) |
| move (refit) vs moved and rebuilt: cube, cone, L-prism, bunny, octree cube, NACA 0012, annulus | geometry, distances, inside flags, box codes, line crossings, cut cells bit for bit |
| refit time / rebuild time | 0.03 (bunny) to 0.4 (tiny meshes) |
| fresh / dead: cube translated, rotated about an off-centre axis, general poses; square, pitching NACA 0012 | 0 errors against the analytic change of membership; the band skips 9–60 % of the second queries |

## Limitations

- **No swept geometry yet.** Swept areas and volumes and the time-averaged
  apertures (#41 step 5) are not implemented. Conservative moving cut cells
  need them.
- **Rigid bodies only.** Deforming bodies, and refitting under deformation,
  are out of scope.
- **Points near the surface.** `state_changes` decides states with the
  point-query sign, so a point within round-off of the wall can differ from
  the cut-cell parity.
- **GPU not validated yet.** The posed device paths run on the host in
  `tests-gnu`; the GPU benches print `posed device vs serial mismatches` for
  the manual gate.
