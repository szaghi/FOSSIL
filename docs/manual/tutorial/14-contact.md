# 14. Contact: distance and collision

Chapter 11 measured the gap between two bodies at one instant with `distance_to`. A moving-body code advances its
bodies by a time step and needs more: the instant inside the step at which two of them first touch. This chapter
throws a small cube, the pebble, at the dragon. Four runs build the answer step by step:

1. sample the gap over the step and see what sampling misses;
2. advance the pebble by conservative steps, each certified by the gap it starts from, to the time of first contact;
3. find the earliest contact among several bodies at once;
4. do the same in 2-D around the NACA 0012 section, then slide a particle along a wall, where the plain bound
   creeps and the directional bound does not.

<<< @/examples/snippets/tut14a_sampling-setup.f90

- The dragon is placed as in chapter 3, without the repair: `distance_to` and the collision queries take any
  triangulated surface, open or closed. The pebble is `cube-01.stl` scaled to side 0.04 and centred on its own
  origin. A body turns about the origin of its own frame, so a pebble centred there spins in place.
- The step takes the pebble from the pose `p0` at s = 0 to `p1` at s = 1: 1.2 along x, straight through the dragon's
  body, while it turns a quarter turn about y. Every run of the chapter uses this setup.

## Sampling misses the contact

<<< @/examples/snippets/tut14a_sampling-sampling.f90

- `pose_interpolate(p0, p1, s)` is the pose at s. All the step queries use this motion: the rotation turns about the
  axis through the body origin by s θ, and the translation is linear in s (see
  [continuous collision detection](/guide/advanced/collision#the-motion-over-a-step)).
- `distance_to(..., intersecting=, pose_other=)` gives the gap at that pose with the trees as built (chapter 11).

<<< @/examples/snippets/tut14a_sampling-dense.f90

<<< @/examples/output/tut14a.ansi{ansi}

- The five samples print gaps from 0.059 to 0.419 and no contact. At s = 0.5 the pebble is inside the dragon, 0.059
  from its surface. It entered the dragon between s = 0.25 and 0.75 and left again, and no sample saw it.
- 200 samples, 1/199 apart, find two intervals of contact: from s = 0.3568 to 0.4171 and from 0.5578 to 0.6080.
  Each tells only that the contact began within the last spacing. A body thinner than its travel between two
  samples, here 0.006, can still pass through unseen. More samples shrink the gap in time, but they never close it.

![the gap between the dragon and the pebble over the step, sampled at 300 instants (yellow curve, s from 0 to 1 left to right), falling to zero over two intervals of contact (ember on the s axis); the five samples of run 1 (teal) read 0.419, 0.125, 0.059, 0.074, 0.359, all positive, and the first contact time_of_impact returns, s = 0.3552, is marked on the axis](../../examples/images/tut14-gap.svg){.fossil-figure}

## Conservative advancement

The guarantee comes from a speed bound. No point of the pebble moves faster than μ (distance per unit of s), so from
a gap d the two bodies cannot touch before s + d/μ. Conservative advancement steps there, measures the gap again and
repeats until the gap falls below a tolerance.

<<< @/examples/snippets/tut14b_toi-bound.f90

- `motion_bound_3d(p0, p1, x, y, z)` is the bound of the motion of `pose_interpolate`: μ = |Δt| + |θ| max ρ, with
  ρ the distance of a vertex (body frame) from the rotation axis. The pebble's μ is 1.2444: 1.2 of translation and
  0.0444 of spin, a quarter turn at ρ = 0.028. A body at rest has μ = 0, and the two bounds add.
- The tolerance defaults to 1e-6 of the larger bounding-box diagonal of the two bodies, here 1.303e-6.

<<< @/examples/snippets/tut14b_toi-steps.f90

- `max_iterations=k` stops `time_of_impact` after k steps. Its `toi` is then the k-th iterate s_k (a certified lower
  bound on the contact) and `distance` the gap there. The loop prints the iterates one by one.

<<< @/examples/output/tut14b.ansi{ansi}

- The first step covers 0.336 of the 0.355 that separate the pebble from the contact: the pebble flies almost
  straight at the dragon, so its closing speed is almost μ. After that each step leaves 5 to 7 % of the gap, and
  step 5 ends at a gap of 9.784e-7, below the tolerance.
- Each step stays within the gap it starts from, so the pebble cannot tunnel through, whatever its speed. The price
  is that the step shrinks with the gap: a body moving along a wall, not towards it, takes many short steps
  ([Sliding contact](#sliding-contact) below).

![conservative advancement of the pebble towards the dragon's foreleg, one frame per step: the pebble at s_k in ember over the earlier iterates in grey, the closest pair in yellow, and above the path, in ember, the stretch its centre may cover before the next query, (d - tol/2)/mu; step 0 jumps 0.336, then the gaps 2.2e-2, 1.4e-3, 8.5e-5, 5.9e-6 and 9.8e-7 at the contact s = 0.355234](../../examples/images/anim-tut14-advance.svg){.fossil-figure}

## The time of first contact

<<< @/examples/snippets/tut14b_toi-toi.f90

- `time_of_impact(other, pose_self0, pose_self1, pose_other0, pose_other1, toi, status, ...)` runs the loop above
  to its end. The dragon is at rest, so both of its poses are `pose_identity()`. The step is
  (d − tol/2)/(μ_A + μ_B), and the loop stops at a gap d(toi) ≤ tol.
- `status` is `CCD_CONTACT` (with d(toi) ≤ tol), `CCD_SEPARATE` (toi = 1), `CCD_INITIAL_CONTACT` (toi = 0) or
  `CCD_NOT_CONVERGED` (the cap of `max_iterations`, default 1000, was hit, and toi is a certified lower bound).
- `distance`, `point_self`, `point_other`, `facet_self`, `facet_other` are the gap, the witness points (world frame)
  and the closest facets at toi; `iterations` counts the steps. The run prints contact at s = 0.355234 in 5 steps,
  on dragon facet 1085 and pebble facet 4, the two witness points equal to 4 decimals.
- `directional=.true.` adds the directional bound of C2A: per pair of triangles, the closing speed projected on the
  pair's separating direction. The bound holds for each pair of triangles, which is convex, and its minimum over the
  pairs holds for the whole body. The loop takes the larger of the two certified steps. On this head-on throw it saves
  one step of five and costs 789114 exact pair evaluations (`pairs`). It stays off by default.

![the placed dragon and the pebble's straight path through its body along x: the pebble at the five sampled instants s = 0, 0.25, 0.5, 0.75, 1 (grey, the one at s = 0.5 hidden inside the dragon), none touching it, and the pebble at the first contact s = 0.3552 that time_of_impact returns (ember), where it first touches the dragon](../../examples/images/tut14-throw.svg){.fossil-figure}

## Several bodies

Now the dragon is at rest and three pebbles fly over the same step: pebble 2 on the throw above, pebble 3 falling
along z through the dragon, pebble 4 flying along x across pebble 3's path, above the dragon.

<<< @/examples/snippets/tut14c_several-several.f90

- `first_impact(surfaces, poses0, poses1, toi, pair, status, tolerance=, directional=, tested=)` gives the earliest
  contact among N bodies over the step and the pair `[i, j]` (i < j) that makes it. Body i moves from `poses0(i)` to
  `poses1(i)`; a body at rest has equal poses.
- Broad phase: the box of each body at s = 0, grown by its μ and by tol/2, contains the body over the whole step, so
  two bodies with disjoint grown boxes never come within the tolerance. The other pairs are ranked by a lower bound
  on their time of impact and run `time_of_impact` in that order. The search stops when the next bound is later
  than the best contact found, and `tested` counts the pairs that ran.

<<< @/examples/snippets/tut14c_several-loop.f90

- Ties go to the smallest pair, so the result equals the loop over all pairs, bit for bit, when both use the same
  tolerance. The default tolerance of `first_impact` is 1e-6 of the largest diagonal of all the bodies, and that of
  `time_of_impact` 1e-6 of the larger of the two: here both get `tolerance=1.e-6_R8P`.

<<< @/examples/snippets/tut14c_several-near.f90

- `proximity_pairs_of(surfaces, cutoff, pairs, distances)` is the instantaneous screen: the pairs of surfaces closer
  than `cutoff` where they stand, with their gaps, found by a box broad phase and `distance_to`. Here it runs on the
  bodies moved to the first contact. `proximity_pairs(cutoff, pairs, distances)` does the same for the bodies
  (`facet_body`) of one surface. See [body–body distance](/guide/advanced/body-distance).

<<< @/examples/output/tut14c.ansi{ansi}

- Pebbles 3 and 4 touch first, at s = 0.087759, before either reaches the dragon. `first_impact` ran 2 of the 6
  pairs. The loop over all of them finds three contacts (1-2 at 0.355234, the throw of run 2; 1-3 at 0.406194;
  3-4) and the same earliest toi, bit for bit.
- At that instant two pairs are within 0.3: pebbles 3 and 4 at 6.6e-7 (the contact, below the tolerance) and the
  dragon and pebble 4 at 0.258.

![the dragon at rest (1) and three pebbles over one step: their straight paths in grey, each pebble at s = 0 translucent and labelled, and each at the first contact s = 0.0878 that first_impact returns: pebbles 3 and 4 in ember where they touch above the dragon, pebble 2 in teal still far from the dragon on its throw](../../examples/images/tut14-several.svg){.fossil-figure}

## The 2-D twins

<<< @/examples/snippets/tut14d_curves-curves.f90

- `curve%distance_to(other, distance=, point_self=, point_other=, segment_self=, segment_other=, intersecting=,
  pose_self=, pose_other=)` is the 2-D gap, with `rigid_pose_2d_t` poses from `pose_2d(angle=, translation=)`.
  Curves are small, so every segment pair is evaluated behind a box test, and contact is decided by the exact
  `orient2d`.
- `curve%time_of_impact` and `first_impact_2d(curves, poses0, poses1, ...)` take the arguments of the 3-D calls, with
  `segment_self`, `segment_other` for the facets (and `motion_bound_2d` for the bound). A body is a `curve_object`:
  here the section of `naca0012-ascii.stl` at mid-span (`slice_surface`) and three discs built with `adopt_loops`.

<<< @/examples/output/tut14d.ansi{ansi}

- At s = 0 particle 2 is 0.179 above the section, closest to segment 47. It falls onto the upper side and touches it
  at s = 0.322097 after 14 steps. That is the first contact of the four curves; `first_impact_2d` tested 5 of the 6
  pairs.

![the NACA 0012 section at rest and three particles of radius 0.03 on straight paths while each turns half a turn: their outlines at s = 0 in grey, at the first contact s = 0.3221 that first_impact_2d returns in teal, the particle that touches the upper side of the section in yellow with the witness point on the section](../../examples/images/tut14-particles.svg){.fossil-figure}

## Sliding contact

<<< @/examples/snippets/tut14d_curves-sliding.f90

- The plain bound cannot tell motion along a wall from motion towards it: a particle moving at speed μ is assumed to
  close the gap at μ. At a gap of 0.02 it advances by about the gap per step and needs 41 steps to cross a step in
  which it never comes closer (last four lines of the run above).
- At a gap of 4e-6, four tolerances, each step is about 3.5e-6 long: the loop hits the cap of 1000 steps at
  s = 0.004375 and returns `CCD_NOT_CONVERGED` with that certified lower bound.
- The directional bound sees that the particle does not approach the wall and crosses the whole step in one, at
  either gap. Pass `directional=.true.` for bodies that slide or roll along each other: resting contact, particles in
  a channel, a valve against its seat.

![a particle of radius 0.03 sliding 0.8 along a flat wall at a gap of 0.02 over one step: top, the plain bound's 41 iterates, each disc outline about one gap ahead of the previous one, the last in yellow; bottom, the directional bound's single step from s = 0 (grey) to s = 1 (yellow)](../../examples/images/tut14-sliding.svg){.fossil-figure}

::: tip What you learned
A gap sampled in time can miss a contact. `motion_bound_3d` bounds the speed of a body over a step, and
`time_of_impact` advances by the steps that bound certifies, to the first contact. Add `directional=.true.` for
bodies sliding along each other. `first_impact` and `first_impact_2d` give the earliest contact among several bodies,
equal to the loop over all pairs, and `proximity_pairs_of` the pairs close at one instant. Reference:
[Continuous collision detection](/guide/advanced/collision), [Body–body distance](/guide/advanced/body-distance).
:::

Next: [15. Body-fitted grids](./15-body-fitted).
