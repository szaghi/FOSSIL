# 14. Contact: distance and collision

Chapter 11 measured the gap between two bodies at one instant with `distance_to`. A moving-body code advances its
bodies by a time step and needs more: the instant inside the step at which two of them first touch. Sampling the step
with `distance_to` cannot give it, because a fast or thin body passes through the other between two samples. This
chapter throws a small cube, the pebble, at the dragon, finds its first contact with a guarantee, does the same for
several bodies at once, and repeats it in 2-D around the NACA 0012 section.

<<< @/examples/snippets/tut14_contact-setup.f90

- The dragon is placed as in chapter 3, without the repair: `distance_to` and the collision queries take any
  triangulated surface, open or closed. The pebble is `cube-01.stl` scaled to side 0.04 and centred on its own
  origin. A body turns about the origin of its own frame, so a pebble centred there spins in place.

## Sampling misses the contact

<<< @/examples/snippets/tut14_contact-sampling.f90

- The step takes the pebble from the pose `p0` at s = 0 to `p1` at s = 1: 1.2 along x, straight through the dragon's
  body, while turning a quarter turn about y.
- `pose_interpolate(p0, p1, s)` is the pose at s, the motion of every step query: the rotation turns about the axis
  through the body origin by s θ, the translation is linear in s (see
  [continuous collision detection](/guide/advanced/collision#the-motion-over-a-step)).
- `distance_to(..., pose_other=)` gives the gap at that pose with the trees as built (chapter 11). The five samples
  print gaps from 0.059 to 0.419 and no contact: at s = 0.5 the pebble is inside the dragon, 0.059 from its surface:
  between s = 0.25 and 0.75 it entered the dragon and left it again, and no sample saw it.

## The time of first contact

<<< @/examples/snippets/tut14_contact-toi.f90

- `time_of_impact(other, pose_self0, pose_self1, pose_other0, pose_other1, toi, status, ...)` gives the first s in
  [0, 1] at which the two bodies come within the tolerance, by conservative advancement: no point of a body moves
  faster than its bound μ = |Δt| + |θ| max ρ (ρ the distance of a vertex to the rotation axis), so from a gap d the
  bodies cannot touch before s + (d − tol/2) / (μ_A + μ_B). Each step is bounded by the gap it starts from, so the
  pebble cannot tunnel through, whatever its speed.
- `status` is `CCD_CONTACT` (with d(toi) ≤ tol), `CCD_SEPARATE` (toi = 1), `CCD_INITIAL_CONTACT` (toi = 0) or
  `CCD_NOT_CONVERGED` (the cap of `max_iterations`, default 1000, was hit: toi is then a certified lower bound). The
  tolerance defaults to 1e-6 of the larger bounding-box diagonal of the two bodies.
- `distance`, `point_self`, `point_other`, `facet_self`, `facet_other` are the gap, the witness points (world frame)
  and the closest facets at toi; `iterations` counts the advancement steps.
- `directional=.true.` adds the directional bound of C2A: per pair of triangles, the closing speed projected on the
  pair's separating direction. The bound holds for each pair of triangles, which is convex, and its minimum over
  the pairs holds for the whole body. The loop takes the larger of the two certified steps.
  `pairs` counts the triangle pairs it evaluated.

![the placed dragon and the pebble's straight path through its body along x: the pebble at the five sampled instants s = 0, 0.25, 0.5, 0.75, 1 (grey, the one at s = 0.5 hidden inside the dragon), none touching it, and the pebble at the first contact s = 0.3552 that time_of_impact returns (ember), where it first touches the dragon](../../examples/images/tut14-throw.svg){.fossil-figure}

## Several bodies

<<< @/examples/snippets/tut14_contact-several.f90

- `first_impact(surfaces, poses0, poses1, toi, pair, status, tolerance=, directional=, tested=)` gives the earliest
  contact among N bodies over the step and the pair `[i, j]` (i < j) that makes it. Body i moves from `poses0(i)` to
  `poses1(i)`; a body at rest has equal poses.
- Broad phase: the box of each body at s = 0, grown by its μ and by tol/2, contains it over the whole step, so two
  bodies with disjoint grown boxes never come within the tolerance. The other pairs are ranked by a lower bound on
  their time of impact and run `time_of_impact` in that order; the search stops when the next bound is later than the
  best contact found. `tested` counts the pairs that ran.
- Ties go to the smallest pair, so the result equals the loop over all pairs, bit for bit, when both use the same
  tolerance. The default tolerance of `first_impact` is 1e-6 of the largest diagonal of all the bodies, and that of
  `time_of_impact` 1e-6 of the larger of the two: here both are passed `tolerance=1.e-6_R8P`.

## The 2-D twins

<<< @/examples/snippets/tut14_contact-curves.f90

- `curve%distance_to(other, distance=, point_self=, point_other=, segment_self=, segment_other=, intersecting=,
  pose_self=, pose_other=)` is the 2-D gap, with `rigid_pose_2d_t` poses from `pose_2d(angle=, translation=)`.
  Curves are small, so every segment pair is evaluated behind a box test; contact is decided by the exact `orient2d`.
- `curve%time_of_impact` and `first_impact_2d(curves, poses0, poses1, ...)` have the arguments of the 3-D calls, with
  `segment_self`, `segment_other` for the facets. A body is a `curve_object`: here the section of
  `naca0012-ascii.stl` at mid-span (`slice_surface`) and three discs built with `adopt_loops`.

![the NACA 0012 section at rest and three particles of radius 0.03 on straight paths while each turns half a turn: their outlines at s = 0 in grey, at the first contact s = 0.3221 that first_impact_2d returns in teal, the particle that touches the upper side of the section in yellow with the witness point on the section](../../examples/images/tut14-particles.svg){.fossil-figure}

<<< @/examples/snippets/tut14_contact-sliding.f90

- The plain bound cannot tell motion along a wall from motion towards it. A particle sliding along a flat wall at a
  gap of 4 tolerances advances by about (gap − tol/2)/μ per step and hits the cap of 1000 steps; the directional bound
  sees that the particle does not approach the wall and crosses the whole step in one.

## Running it

<<< @/examples/output/tut14.ansi{ansi}

- The samples miss the contact; `time_of_impact` finds it at s = 0.355234 in 5 steps, at a gap of 9.784e-7, below
  the default tolerance (1e-6 of the larger bounding-box diagonal).
- On this head-on throw the directional bound saves one step of five and costs 789114 exact pair evaluations: it
  pays off on sliding contact, not on approach. It stays off by default.
- Among the dragon and three pebbles the first contact is between two pebbles (bodies 3 and 4) at s = 0.087759,
  before either reaches the dragon. `first_impact` ran 2 of the 6 pairs; the loop over all of them gives the same
  toi, bit for bit.
- In 2-D the particle falling onto the upper side of the section touches it at s = 0.322097 (14 steps) and is the first
  contact of the four curves; `first_impact_2d` tested 5 of the 6 pairs.
- Sliding at 4e-6 from the wall, the plain bound stops at the cap with s ≥ 0.004375, a certified lower bound; the
  directional bound returns `CCD_SEPARATE` in 1 step.

::: tip What you learned
A gap sampled in time can miss a contact; `time_of_impact` gives the first contact over a step with a guarantee,
`directional=.true.` for bodies sliding along each other, `first_impact` and `first_impact_2d` the earliest contact
among several bodies, equal to the loop over all pairs. Reference: [Continuous collision
detection](/guide/advanced/collision), [Body–body distance](/guide/advanced/body-distance).
:::

Next: [15. Body-fitted grids](./15-body-fitted).
