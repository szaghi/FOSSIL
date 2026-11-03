---
title: Immersed-boundary queries
---

# Immersed-boundary queries

## What it does

Image-point, inverse Lax–Wendroff and wall-model closures of an
immersed-boundary (IB) solver need more than the signed distance φ. They
need:

- the **wall intercept**: the closest point;
- the **wall normal**;
- **which body** the wall belongs to;
- **image or probe points**, with the guarantee that the stencil does not
  reach through another wall.

A normal from ∇φ by finite differences degrades at sharp edges and near the
medial axis (Wang et al. 2022). FOSSIL computes these quantities exactly, in
batches, on the CPU or the device, for surfaces (3-D) and
[curves](/guide/advanced/2d-geometry) (2-D).

| Query | 3-D (`surface_stl_object`) | 2-D (`curve_object`) |
|---|---|---|
| Rich batch: distance, closest point, closest facet or segment, feature, normal, body | `query_many`, `query_many_device` | `query_many`, `query_many_device` |
| One point | `query_point`, `compute_distance(..., closest_point=, normal=)` | `compute_distance(..., normal=, body=)` |
| Probe along the normal, with visibility | `probe_many` | `probe_many` |
| Mirror through the closest point, with visibility | `image_points` | `image_points` |
| Several foot points (medial axis, thin gaps) | `query_many(..., ambiguous=, second_foot_points=)` | same |
| Bodies | `get_bodies_number`, `facet_body`, `get_shells_number`, `facet_shell` | `get_bodies_number`, `loop_body` |
| Vector level set | `distance_field_object%compute(..., store_closest=.true.)` | same, 2-D field |

## Usage

```fortran
program ghost_points
use fossil, only : surface_stl_object, distance_field_object, FEATURE_FACE
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)    :: body
type(distance_field_object) :: field
type(vector_R8P)            :: ghost(2)
real(R8P)                   :: d(2), foot(3, 2), nrm(3, 2), img(3, 2), probe(3, 2), x(32)
integer(I4P)                :: feat(2), bod(2), k
logical                     :: vis(2), amb(2)

call body%load_from_file(file_name='body.stl', guess_format=.true.)
ghost(1) = vector_R8P(0.5_R8P, 0.5_R8P, 0.9_R8P)
ghost(2) = vector_R8P(0.1_R8P, 0.2_R8P, 0.3_R8P)

! wall intercept, normal, feature and body of each ghost point
call body%query_many(points=ghost, distances=d, closest_points=foot, features=feat, normals=nrm, bodies=bod, &
                     ambiguous=amb)
! mirror images (2 foot - x) and probes at sqrt(2) dx, with visibility
call body%image_points(points=ghost, images=img, visible=vis)
call body%probe_many(points=ghost, lengths=[(sqrt(2._R8P) * 0.05_R8P, k = 1, 2)], probe_points=probe, visible=vis)

! vector level set: closest point, facet and normal stored on the band cells
do k = 1, 32
   x(k) = -0.5_R8P + (k - 0.5_R8P) / 16._R8P
enddo
call field%initialize(x=x, y=x, z=x)
call field%compute(surface=body, store_closest=.true.)
endprogram ghost_points
```

## Definitions

- **Normal.** The angle-weighted pseudo-normal of the closest feature
  (Bærentzen–Aanæs): the facet normal, or an edge or vertex pseudo-normal.
  It is the same vector the sign test uses, and it is the reconstruction
  direction at edges that Seshadri & De (2020) recommend. In 2-D it is the
  segment normal or the vertex bisector.
- **Feature.** In 3-D: `FEATURE_FACE`, `FEATURE_EDGE` or `FEATURE_VERTEX`,
  with the local edge or vertex index. In 2-D: the segment, plus the vertex
  when the closest point is one.
- **Body.**
  - In 3-D, a **shell** is a set of edge-connected facets. A shell at even
    nesting depth starts a body; a shell at odd depth (a cavity) joins the
    body that contains it; an island inside a cavity is a body of its own.
  - In 2-D: an outer loop with its holes; open loops are bodies of their
    own.
  - Containment is tested only where one bounding box contains another, so
    scattered particles cost nothing.
- **Image point.** The mirror through the closest point, 2·foot − x.
- **Probe.** The point at a given distance along the normal from the foot
  point: positive outward, negative inward.
- **Visibility.** Whether the segment from the foot point to the image or
  probe point crosses the geometry.
  - **3-D:** a ray from the foot point. Hits within 10⁻⁹ of the box
    diagonal are ignored: those are the facets containing the foot point.
  - **2-D:** exact, with `orient2d`, skipping the segments that contain the
    foot point.
- **Ambiguous foot point.** Another facet or segment whose closest point is
  separate from the foot point (by more than `foot_tolerance`, default
  10⁻⁹ of the box diagonal) lies within (1 + `rel_tolerance`) of the minimum
  distance (default 10⁻⁶). This happens on the medial axis and in thin gaps.
  FOSSIL reports the second foot point and the caller decides; the
  literature gives no rule.

## Wall distance and field diagnostic

- **`wall_distance_many(points, distances, bodies, cutoff)`** (surface, with
  a device twin; curve) returns the unsigned distance to the nearest wall
  and the body it belongs to (issue #42). Turbulence models and overset
  assembly need it (Oh et al. 2024).
  - The tree traversal starts at the cutoff, so every node farther away is
    pruned.
  - Points with no wall closer than the cutoff get the cutoff and body 0.
  - CPU and device give the same results bit for bit.
- **`distance_field_object%gradient_diagnostic(q, flagged, threshold)`**
  returns q = |1 − |∇φ|| per node.
  - ∇φ is the three-point non-uniform difference, exact on linear fields on
    stretched grids.
  - Nodes with an unknown neighbour get q = −1.
  - `flagged` counts the nodes above the threshold (default 0.05, Matsushita
    & Aoki).
  - q is large by nature on the medial axis.

## Curvature at the foot point

Curved slip walls need the wall curvature where the image point projects
(∂ₙv_t = v_t / R_c; Funada & Imamura, Ciallella et al.). `query_many`
returns it as an optional output (issue #42, CPU):

- **Surfaces: `curvatures(2, n)` = (H, K).** These are the mean curvature
  (k₁ + k₂)/2 and the Gaussian curvature k₁k₂.
  - Both are the per-vertex `mean_curvature` (cotangent Laplacian) and
    `gaussian_curvature` (angle defect), interpolated with the foot's
    barycentric weights. They are rebuilt once per call, at O(facets) cost.
  - H > 0 where the surface is convex seen from the side the normals point
    to: a sphere gives (1/R, 1/R²), a cylinder (1/(2R), 0).
- **Curves: `curvatures(n)`.** Each vertex gets its turning angle over the
  mean length of its two segments, and the value is interpolated along the
  segment.
  - It is positive on outer loops and negative on holes.
  - It is 0 at the ends of an open loop.

Both work with `pose=`, since curvature is invariant under rigid motion.

## Grid-adapted markers

Wall models evaluated at surface points (van Noordt et al.) and
direct-forcing force integration need markers spaced like the grid. Direct
forcing needs a spacing of 0.8–1.2 Δx (Cheylan et al. 2023).

- **`surface%surface_markers(dx, centroids, areas, normals, spacing_fraction,
  spacings, iterations, preserve_features)`** remeshes a copy of the
  surface with `isotropic_remesh` and returns the facet centroids, areas
  and unit normals. The surface itself is not changed.
  - The target edge is √3 Δx: centroids of two equilateral neighbours are
    then Δx apart.
  - One calibration pass rescales that target by Δx over the measured mean
    spacing, because the remesher's edges come out short.
  - `spacing_fraction` reports the fraction of nearest-neighbour spacings
    inside the window.
  - Pass `preserve_features=.true.` on bodies with sharp edges.
- **`curve%curve_markers(dx, centroids, lengths, normals, spacing_fraction,
  spacings)`** resamples every loop exactly along the polyline into
  nint(L/Δx) pieces of equal arc length. The markers sit at the piece
  midpoints, on the polyline.

## Placement rules from the literature

All of them are lengths for `probe_many`:

| Rule | Length |
|---|---|
| Image point with a complete interpolation stencil (Wang et al. 2022) | √2 Δx |
| Image line of order N + 1 (Funada & Imamura 2023) | ≥ √d Δx (d = 2, 3) |
| Second image layer (Constant et al. 2024) | H + 2 √d Δx |
| Fixed-distance probe that sidesteps the small-cell problem (van Noordt et al. 2022) | fixed |

## Results

| Case | Result |
|---|---|
| `query_many` vs `compute_distance`, cube and bunny (1728 points each, pseudo-normal and FWN signs) | every output bit for bit |
| `query_many_device` vs `query_many`, cube and bunny (3456 points each) | bit for bit (device path run on the host) |
| curves: NACA 0012, annulus with island, open square (3000 points) | bit for bit, CPU and device path |
| vector level set: bunny, cube, 256-gon | band values = `query_many`; \|φ\| = \|x − closest\| |
| cube ghost points in face regions (1100) | closest + \|φ\| n = 2·closest − x |
| two cubes 0.2 apart | a probe into the other cube, or an image inside it, is reported not visible |
| cube centre, diagonal plane, mid-gap, L-shape medial points | flagged; generic points not |
| foot curvature, icospheres R = 2 (2, 3, 4 subdivisions) | RMS error of H 1.45e-2, 8.8e-3, 5.6e-3; of K 8.9e-3, 4.6e-3, 2.8e-3; maximum relative error 0.12–0.15 |
| foot curvature, cylinder and flat faces | exact to round-off |
| foot curvature, regular N-gons R = 1 | (π/N)/sin(π/N) exactly: the error against 1/R falls by 4.0 per halving of h |
| markers, unit icosphere, Δx = 0.06 | 1962 markers, area within 0.6 %, 91 % of spacings in [0.8, 1.2] Δx |
| markers, unit cube, Δx = 0.1, features kept | area 6 to round-off, 72 % of spacings in the window |
| curve markers, square and 64-gon | on the polyline (distance ≤ 6.6e-17), equal lengths summing to the perimeter |

## Limitations

- **GPU not validated yet.** The device paths (`query_many_device`, 3-D and
  2-D) run bit-exact on the host in `tests-gnu`. The GPU benches print
  `rich device vs serial max err` and `rich id mismatches` for the manual
  gate, which has not run yet.
- **Multiple foot points and foot curvature are CPU only.**
- **3-D foot curvature does not converge pointwise.** The per-vertex
  operators use barycentric areas. On icospheres the RMS error decreases
  with refinement, but the maximum stays near 15 %, at the valence-5
  vertices. Cylinders and planes are exact.
- **Marker spacing follows the remesher.** Most spacings land inside
  [0.8, 1.2] Δx (91 % on a sphere, 72 % on a cube), but not all.
  `spacing_fraction` reports how many do.
- **3-D visibility uses ray casting.** A ray grazing an edge exactly may
  count or miss a hit, as in any ray caster. The 2-D test is exact.
