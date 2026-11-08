---
title: Sign algorithms at sharp features
---

# Sign algorithms at sharp features

## What it measures

Every inside/outside decision in FOSSIL goes through a sign algorithm:
- the signed distance;
- the IB queries;
- the box classification;
- the tag bands;
- the field signs.

Near sharp edges and vertices the literature disagrees on what is safe:
- Seshadri & De reject signs from the distance and use ray tracing.
- Mittal et al. use the normal of the closest element and report no
  misclassification.
- Neither paper measures error rates.

This benchmark measures them (issue #39). It covers every 3-D algorithm,
every 2-D algorithm, and that naive closest-element sign.

| Benchmark | Test | What |
|---|---|---|
| 3-D | `fossil_test_sign_bench` | 9 bodies, 8198 points; 5 dirty variants of each (40990 points) |
| 2-D | `fossil_test_sign_bench_2d` | 12 curves, 15258 points; 3 dirty variants of each |
| Device | `fossil_test_sign_bench_device` | device signs = CPU signs on the benchmark points |

## Method

- **Bodies (3-D):**
  - a cube;
  - wedges of 5°, 30°, 90° and 150°;
  - a plate of thickness 10⁻³;
  - a 10° pyramid;
  - a 15° cone with 32 faces;
  - an L-shaped prism (reentrant edge, saddle vertices).

  All are procedural, watertight and outward-oriented (`fossil_sign_bench_bodies`).
- **Curves (2-D):**
  - a square;
  - wedges of 1° to 150°;
  - an L-shape;
  - a cusp: two circular arcs tangent at the origin, discretized to a 0.01° corner;
  - a C-shape with a 10⁻³ slot;
  - a square annulus;
  - two sections of the 3-D bodies.
- **Points.** Each query point is a feature point x plus δ·d:
  - x is on a face, a flat diagonal, an edge (convex or reentrant), or a
    vertex (convex or saddle); in 2-D, a corner;
  - d lies in the inside or the outside sector of x, at least a sixth of the
    opening angle away from every incident face;
  - δ runs from 10⁻¹ to 10⁻¹² L (L the bounding-box diagonal), and never
    exceeds a quarter of the distance to the nearest face that does not
    contain x.

  So the side is known **by construction**: no 3-D predicate is trusted as
  the reference. In 2-D the exact winding number (exact `orient2d`) agrees
  with the construction on every clean point.
- **Dirty variants.**
  - 3-D: a facet removed (a hole), a facet flipped, a facet duplicated,
    needle triangles (aspect about 10⁶) along an edge, a T-junction.
  - 2-D: a gap (open loop), a duplicated vertex, a segment of 10⁻⁹.
  - The ground truth stays the clean geometry, and 3-D errors are split by
    distance to the defect (near: under 0.1 L).

## Results, 3-D

**Clean bodies**, errors out of 8198 points, all decades:

| Algorithm | face | flat edge | convex edge | reentrant edge | convex vertex | saddle vertex |
|---|---|---|---|---|---|---|
| `SIGN_PSEUDO_NORMAL` | 0 | 0 | 0 | 0 | 0 | 0 |
| `SIGN_RAY_INTERSECTIONS` | 0 | 0 | 0 | 0 | 0 | 0 |
| `SIGN_SOLID_ANGLE` | 0 | 0 | 0 | 0 | 0 | 0 |
| `SIGN_FAST_WINDING_NUMBER` | 0 | 0 | 0 | 0 | 0 | 0 |
| naive closest-facet normal | 0 | 0 | 173 (every decade) | 11 | 405 (every decade) | 0 |

**Dirty variants**, errors near the defect / far from it:

| Defect | pseudo-normal | ray | solid angle | fast winding |
|---|---|---|---|---|
| hole | 227 / 0 | 305 / 0 | 0 / 0 | 0 / 0 |
| flipped facet | 945 / 0 | 0 / 0 | 830 / 0 | 819 / 0 |
| duplicated facet | 123 / 1 | 305 / 0 | 0 / 0 | 0 / 0 |
| needles | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| T-junction | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |

Before the needle fix (below) the pseudo-normal had 2 "far" errors in each
variant, its 2 clean-body errors repeated. The fix also moved 14 errors near
the defects of the plate: its defect facet is a needle, and the foot now lands
on the flipped (or duplicated) needle itself, whose normal is wrong by
construction.

**Cost per query** on these small meshes. Only the ratios mean anything:
the ray and solid-angle tests are O(N), the winding number O(log N).

| Algorithm | µs per query |
|---|---|
| pseudo-normal (fused with the distance search) | 0.14 |
| ray | 0.23 |
| solid angle | 0.23 |
| fast winding number | 0.25 |

## Results, 2-D

**Clean curves:** 0 errors out of 15258 points for the pseudo-normal, the
fast winding number and the exact winding number. That holds for every
feature kind (segment, flat vertex, convex, acute ≤ 30° and reentrant
corner), the 0.01° cusp included, and every decade down to 10⁻¹² L.

**Dirty variants** (errors out of 15258 points):

| Defect | pseudo-normal | fast winding | exact winding |
|---|---|---|---|
| gap (open loop) | 527 | 5 | 7893 (it counts closed loops only) |
| duplicated vertex (zero-length segment) | 57 | 0 | 0 |
| micro-segment (10⁻⁹) | 0 | 0 | 0 |

## What the failures are

- **The distance itself was the first failure.** The first run gave the
  pseudo-normal 2220 errors on the clean bodies, starting at 10⁻⁷ L at edges
  and 10⁻⁹ L on faces.
  - **Cause:** the point–triangle distance came from the expanded quadratic
    form of the Eberly solve. Its absolute error is about ε·L², so it was
    noise below √ε·L ≈ 10⁻⁸ L and could come out negative, which flipped the
    sign.
  - **Fix:** it is now the norm of the residual, accurate to about ε·L.
    The 2220 errors fell to 2 (`fossil_test_distance_tiny` guards it). This
    is the reason to measure before trusting a default.
- **Needle triangles (fixed).** The last 2 pseudo-normal errors sat on the
  1 × 10⁻³ side faces of the plate, at δ = 10⁻¹⁰ and 10⁻¹¹ L.
  - **Cause:** the Eberly solve goes through the Gram system, whose error
    grows like ε/|n|² (n = e₁₂ × e₁₃). It put the foot on the edge,
    1.41·10⁻¹⁰ away, instead of the face point 10⁻¹⁰ away. There the edge
    pseudo-normal is orthogonal to p − foot, so the sign flipped.
  - **Fix:** triangles of aspect above about 100 (det < 10⁻⁴ max(a, c)²)
    get the foot from triple products with n, with error about ε/|n|. If
    that foot lies outside the triangle, the nearest clamped edge projection
    is used instead.
  - Other triangles keep the Eberly solve, bit for bit: the robust path costs
    20 % more per query on dragon-fine.
  - Errors went from 2 to 0. `fossil_test_distance_tiny` now includes the
    plate, which gives 2 sign errors with the old solve.
- **The naive closest-facet normal is wrong at edges and vertices at every
  distance.** It fails at about 5 % of convex-edge points and 10 % of
  convex-vertex points, at every distance from the largest δ down. So the method as Mittal et al. describe it does
  misclassify. The angle-weighted pseudo-normal (Bærentzen–Aanæs) is what
  makes the distance-based sign sound.

  ![the section of a 30-degree wedge at its apex edge, 72 points on a circle of radius 0.05 about the edge coloured by the side each sign gives (teal inside, ember outside): the naive closest-facet normal puts an arc of 12 points behind the apex inside, the pseudo-normal none](../../examples/images/ref-sign-benchmark-naive.svg){.fossil-figure}

  Behind an acute edge the closest point is the edge itself, and the normal cone there spans 150°: the normal of
  either facet makes an obtuse angle with part of it. The edge pseudo-normal is the bisector, at most 75° from any
  direction of the cone.
- **Holes and duplicated facets** break the ray test: rays through the
  hole, double hits. They also break the pseudo-normal, because the
  connectivity near the defect is wrong. They do not break the solid angle
  or the winding number: next to a missing planar facet the winding number
  tends to 0.5 from the correct side, so the 0.5 threshold still separates
  the two sides.

  ![the section y = 0.5 of the wedge with the upper-face facet on its apex edge removed (the removed part in ember), a grid of 960 points coloured by the side each algorithm gives, one frame per algorithm: the pseudo-normal 111 wrong, the ray test 33, the solid angle and the fast winding number none](../../examples/images/ref-sign-benchmark-hole.svg){.fossil-figure}

  On the holed wedge the pseudo-normal labels a block of outside points above and left of the apex as inside: their
  closest point is on the apex edge, now a boundary edge with one facet, so its pseudo-normal is the lower face's
  normal alone, the naive sign of the figure above. The ray test fails inside, next to the apex, where the rays
  leave through the hole without a crossing.
- **A flipped facet** breaks every sign test except the ray parity near the
  facet. That is the right failure: the orientation is the input. Running
  `sanitize` (BFS orientation repair) removes it.
- **2-D:** with closed, repaired curves, every algorithm is exact on these
  cases. A gap defeats the exact winding number (it ignores open loops) and
  most of the pseudo-normal, but hardly the fast winding number. A
  zero-length segment gives the pseudo-normal an undefined normal; `repair`
  merges such vertices.

## Example

The method on one body: a wedge of 30° built in code (8 facets), its upper-face facet on the apex edge removed (a
hole) or flipped. The points are $x + \delta d$, as in the benchmark: $x$ on the upper face or on the apex edge,
$d$ in the inside or the outside sector at least a sixth of its opening away from every face, $\delta$ from 10⁻¹ to
10⁻¹². The program is `docs/examples/src/ref_sign_benchmark.f90`; the loop and the naive sign:

<<< @/examples/snippets/ref_sign_benchmark-bench.f90

<<< @/examples/snippets/ref_sign_benchmark-naive.f90

<<< @/examples/output/ref_sign_benchmark.ansi{ansi}

- **Clean:** the four algorithms make no error on the 204 points; the naive closest-facet normal makes 36, 3 at
  every decade from 10⁻¹ to 10⁻¹². Its error does not shrink with δ: it is a wrong direction, not round-off.
- **Hole:** the pseudo-normal (60) and the ray test (71) fail; the solid angle and the fast winding number do not,
  as in the table above.
- **Flipped facet:** every test fails except the ray parity, which does not read the orientation (117, 94 and 83
  errors for the pseudo-normal, the solid angle and the winding number). `sanitize` repairs the orientation first.
- `is_point_inside(point=, sign_algorithm=)` is the call; `compute_distance(..., facet_index=, closest_point=)`
  gives what the naive sign needs. The naive counts on the dirty variants include its 36 clean errors.

The recipe [Get a reliable sign on a mesh with holes](/manual/cookbook#get-a-reliable-sign-on-a-mesh-with-holes) is
the practical side of this page: pseudo-normal against winding number on an open cube.

## Device

The device paths (`distance_many_device`, `query_many_device`, 3-D and 2-D,
pseudo-normal and fast winding number) give the same signed distances as
the CPU on all benchmark points: 0 mismatches, run on the host. The GPU
benches print `sign bench device vs serial mismatches` for the manual gate,
which has not run yet.

## Recommendation

- **Distance and IB queries** (`distance`, `distance_many`, `compute_distance`,
  `query_many`, the band of `distance_field_object`): keep
  `SIGN_PSEUDO_NORMAL` as the default.
  - It comes free with the closest-point search, which these queries do
    anyway.
  - It is exact on watertight, oriented meshes down to 10⁻¹² L, needles
    of aspect 10³ included.
  - On geometry that cannot be sanitized (holes, duplicates), pass
    `SIGN_FAST_WINDING_NUMBER`: it was right on every dirty variant except
    the flipped facet.
- **Inside/outside without a distance** (`box_classify_many`, `tag_block`):
  the sign costs one point query per block or connected component, so its
  price does not matter. These default to `SIGN_FAST_WINDING_NUMBER`, for
  the robustness on dirty input (decided in issue #39). Pass
  `SIGN_PSEUDO_NORMAL` for the previous behaviour.
- **2-D:** the same split. On curves, call `repair` before any query.

## Limitations

- **Small meshes.** The cost figures are relative. The ray and solid-angle
  tests do not scale (O(N) per query).
- **No GPU run yet.** The device = CPU result holds on the host only until
  the manual GPU gate runs.
