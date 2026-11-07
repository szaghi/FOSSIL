---
title: Continuous collision detection
---

# Continuous collision detection

A moving-body code advances its bodies by a time step and needs the instant at which two of them first touch inside
that step. [Body–body distance](./body-distance) answers at one instant. Sampling the step with it can miss a
contact: a thin or fast body passes through the other between two samples (**tunnelling**). Issue #73 (roadmap #68,
WP4) adds the **time of first contact with a guarantee**, in the module `fossil_ccd`.

::: info Design and reference (issue #73)
The llm-wiki bundles have no source on body contact: every "collision" there is a collision term of a kinetic or
Riemann solver. The method is conservative advancement *(gk: Mirtich 1996; Zhang, Redon, Lee & Kim 2006, C2A; the FCL
and Bullet libraries)*. Statements from general knowledge are marked *(gk)*. CPU only: one query per pair is a
sequence of distance queries, and the device gains nothing there.
:::

## The motion over a step

A step takes each body from a pose P0 = (R0, t0) at s = 0 to a pose P1 = (R1, t1) at s = 1, along the motion of
`pose_interpolate` ([moving bodies](./moving-bodies), issue #41), the same motion as the swept cut-cell geometry:

- the rotation R(s) = Q(s) R0, where Q(s) turns about the fixed unit axis k of R1 R0ᵀ by the angle s θ;
- the translation t(s) = t0 + s Δ, with Δ = t1 − t0.

A body point x_b (body frame) is at x(s) = Q(s) R0 x_b + t(s). The rotation turns about the axis through the
**body-frame origin**, which moves with t(s). This is not Chasles' screw, whose axis is fixed in space. In 2-D,
Q(s) is the in-plane rotation by s θ, with θ ∈ (−π, π].

### The speed bound

The velocity is dx/ds = θ k × y + Δ with y = Q(s) R0 x_b. Q(s) turns about k, so |k × y| = ρ(x_b), the distance of
R0 x_b to the axis through the body origin, does not depend on s. Hence, for every s,

  |dx/ds| ≤ |Δ| + |θ| ρ(x_b).

ρ is a convex function of x_b (a distance to a line), so its maximum over the body (the hull of the vertices) is
reached at a vertex. The bound of a body over the step is

  **μ = |Δ| + |θ| max_v ρ_v**,

and no point of the body moves farther than μ |s2 − s1| between s1 and s2. In 2-D, ρ_v = |x_v| (the axis is the
body origin). A body whose vertices lie far from its own frame origin has a large μ when it turns: it is its real
motion, not slack. Measured in step 2: the bound is attained at the vertex farthest from the axis.

## Conservative advancement

The distance d(s) between bodies A and B is Lipschitz: each point of A moves at most μ_A δ, each point of B at most
μ_B δ, so

  d(s + δ) ≥ d(s) − (μ_A + μ_B) δ.

From s_k with d_k = d(s_k) > tol, the step

  s_{k+1} = min(1, s_k + (d_k − η tol) / (μ_A + μ_B)),  η = 1/2,

keeps d ≥ η tol on [s_k, s_{k+1}]. The loop stops when d_k ≤ tol (contact at s_k), or at s = 1 with d(1) > tol (no
contact over the step).

**Guarantee (proved by the Lipschitz bound, up to the rounding of the distance query):**

- d ≥ η tol on [0, toi], so toi is no later than the first time the bodies come within η tol, and a fortiori no later
  than the first contact (d = 0);
- d(toi) ≤ tol for a contact;
- no tunnelling, whatever the speed or the thickness: every step is bounded by the distance it starts from;
- at most ⌈(μ_A + μ_B) / ((1 − η) tol)⌉ iterations, since every step with d_k > tol advances at least
  (1 − η) tol / (μ_A + μ_B).

The distance comes from `distance_to(..., pose_self=, pose_other=)` (#67) at the interpolated poses: the trees are
built once in the body frames, and the result is bit-identical to moving the bodies. Its rounding is a few ulps of the
body size, far below the default tolerance.

| Status | Meaning | toi |
|---|---|---|
| `CCD_CONTACT` | the bodies come within `tolerance` | the first iterate with d ≤ tol |
| `CCD_SEPARATE` | no contact over [0, 1] | 1 (d(1) > tol) |
| `CCD_INITIAL_CONTACT` | d(0) ≤ tol | 0 |
| `CCD_NOT_CONVERGED` | the iteration cap (`max_iterations`, default 1000) was hit | the last iterate, still a certified lower bound on the contact time |

The tolerance is optional and defaults to 1e-6 of the larger bounding-box diagonal of the two bodies (the convention
of `deviation`). Grazing approaches, where the bodies close slowly, converge slowly: that is what the cap is for.

### The directional bound (C2A)

The plain bound assumes the worst direction of motion. Along the unit witness direction n of a convex pair, the
separating plane gives d ≥ n · (x_b − x_a), whose closing rate is bounded by the motion projected on n
(Mirtich) *(gk)*. **This holds for convex sets only**: a body that is not convex can close around the other from a
side the plane does not see. FOSSIL's bodies are not convex, so the bound is applied per triangle pair (segment pair
in 2-D), each of them convex, inside the descent of the two trees, as C2A does with per-node bounds *(gk)*:

- a triangle pair at distance d_ij with witness direction n has the closing rate
  c_ij = n · (Δ_A − Δ_B) + |θ_A| ρ_A(T_a) + |θ_B| ρ_B(T_b), with ρ(T) the largest ρ of its vertices and n
  pointing from T_a to T_b (the gap along a fixed n bounds the pair's distance from below at every s, and its rate
  is at least −c_ij);
- its safe step is (d_ij − η tol) / c_ij (unbounded if c_ij ≤ 0), and the step is the minimum over all pairs;
- a node pair whose plain-bound step (box distance over the node μ) is already above the best is pruned.

`directional=.true.` selects it. The plain bound stays the default unless the measures of step 5 show the descent
pays for itself.

## Several bodies

`first_impact(surfaces, poses0, poses1, ...)` gives the earliest contact among N bodies and its pair:

- **broad phase:** the box of each body at s = 0 grown by μ contains it over the whole step (conservative); only pairs
  of overlapping swept boxes go on;
- every remaining pair has a lower bound on its contact time, (box distance at s = 0 − η tol) / (μ_A + μ_B); pairs are
  run in the order of these bounds, and the search stops when the next bound exceeds the best contact found;
- ties go to the lexicographically smallest pair, so the result equals the brute force over all pairs, bit for bit.

## API

| Routine | Arguments |
|---|---|
| `motion_bound_3d(p0, p1, x, y, z)` / `motion_bound_2d(p0, p1, u, v)` | μ of a vertex set |
| `ccd_advance(evaluator, mu, tolerance, toi, status, distance, iterations, max_iterations)` | the loop over an abstract distance evaluator (`ccd_evaluator_t`) |
| `surface%time_of_impact(other, pose_self0, pose_self1, pose_other0, pose_other1, toi, status, ...)` | optional: `tolerance`, `distance`, `point_self`, `point_other`, `facet_self`, `facet_other`, `iterations`, `max_iterations`, `directional` |
| `curve%time_of_impact(...)` | the 2-D twin (`segment_self`, `segment_other`) |
| `first_impact(surfaces, poses0, poses1, toi, pair, status, tolerance)` / `first_impact_2d` | several bodies |

## Proved and measured

- **Proved:** the speed bound, the Lipschitz bound, the guarantee and the iteration bound above (up to the rounding
  of the distance query).
- **Measured (steps 4-7):** toi against a brute-force first contact (dense sampling, then bisection); iterations of
  the plain and the directional bound on head-on, oblique and grazing approaches; time per query.

## Left out

- Exact CCD of straight vertex trajectories (vertex-face and edge-edge coplanarity cubics, Brochu, Edwards & Bridson
  2012) *(gk)*: it fits deforming meshes, not rigid rotation.
- Penetration depth and contact response: the solver's.
- Deformable bodies and self-collision over time.
- Device kernels.
