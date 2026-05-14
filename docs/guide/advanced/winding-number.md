---
title: Generalized winding number
---

# Generalized / fast winding number

## What it does

Returns a continuous scalar `w(q)` that measures **how many times the
input surface wraps around a query point** `q`. For a closed,
outward-oriented surface:

- `w(q) ≈ 1.0` strictly inside the volume
- `w(q) ≈ 0.0` strictly outside
- intermediate values on the surface boundary

The point of the *generalized* winding number is what happens **when the
surface isn't watertight**: holes, missing facets, near-coincident
duplicate triangles, mixed-orientation patches. The legacy sign
algorithms in FOSSIL (`SIGN_RAY_INTERSECTIONS`, `SIGN_SOLID_ANGLE`,
`SIGN_PSEUDO_NORMAL`) all degrade unpredictably on dirty input — they'll
report some points as inside when they shouldn't be, or oscillate near
holes. The generalized winding number degrades **gracefully**: a query
point on the symmetry plane of a cube-with-one-face-deleted returns
~0.5 (half-wrapped), not garbage.

This is the *highest-leverage* of FOSSIL's advanced primitives. It directly fixes a
known correctness gap in FOSSIL's existing signed-distance API on the
dirty STLs that motivate `sanitize` in the first place.

## Pipeline

The exact formulation is

$$
w(q) = \frac{1}{4\pi} \sum_{f \in F} \Omega_f(q)
$$

where $\Omega_f(q)$ is the *signed solid angle* subtended by triangle
$f$ at point $q$ (Van Oosterom & Strackee 1983). The sum is over all
facets in the surface.

For a query at $q$:

1. **Brute force** (`beta ≤ 0`): evaluate the per-facet solid angle
   directly. O(N) per query, exact within floating-point precision.
2. **Hierarchical Barnes-Hut** (default, `beta = 2.0`): for each AABB
   tree node, precompute its **dipole moment**
   $\sum_f \mathrm{area}_f \cdot \hat{n}_f$ and centroid. If the query
   point is far from the node (distance > `beta` × node diagonal),
   approximate the entire subtree's contribution with the dipole.
   Otherwise recurse. Reduces a 70k-facet bunny query from ~70 µs to
   ~5 µs without measurable accuracy loss.

The same Barnes-Hut traversal works on both the SAH BVH and the octree
(tree kind is irrelevant to the algorithm — child enumeration goes
through `enumerate_children`).

## API

```fortran
real(R8P) function surface%winding_number(point, beta)
   type(vector_R8P), intent(in)           :: point  ! query point
   real(R8P),        intent(in), optional :: beta   ! Barnes-Hut admissibility (default 2.0)
```

- **`point`** is the 3D query position.
- **`beta`** controls the multipole admissibility ratio. Larger β
  means stricter admissibility (more recursion, more accuracy, slower).
  Default `beta = 2.0` is a balanced setting; production use rarely
  needs tuning. `beta ≤ 0` forces the exact O(N) sum.

The returned `w` is the unrounded scalar. To use it as an
inside/outside test, threshold at `0.5`:

```fortran
inside = surface%winding_number(point=q) > 0.5_R8P
```

## Example

Closed cube, far-field decay, and the dirty-mesh graceful-degradation
property all in one program:

```fortran
program ex_winding_number
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none

type(surface_stl_object) :: cube
type(vector_R8P)         :: p_in, p_out, p_far
real(R8P)                :: w_in, w_out, w_far

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)

p_in  = vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)   ! cube centroid
p_out = vector_R8P(2.0_R8P, 2.0_R8P, 2.0_R8P)   ! outside the bbox
p_far = vector_R8P(1.0e3_R8P, 0._R8P, 0._R8P)   ! far field

w_in  = cube%winding_number(point=p_in)         ! ≈ 1.0
w_out = cube%winding_number(point=p_out)        ! ≈ 0.0
w_far = cube%winding_number(point=p_far)        ! < 1e-10

print '(A,F12.9)', 'inside  : w = ', w_in
print '(A,F12.9)', 'outside : w = ', w_out
print '(A,ES12.4)', 'far     : w = ', w_far
endprogram ex_winding_number
```

To force exact (non-hierarchical) evaluation — useful for ground-truth
testing or for inputs where the AABB tree isn't trustworthy — pass
`beta = -1.0` or any non-positive value:

```fortran
program ex_winding_number_exact
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none

type(surface_stl_object) :: cube
type(vector_R8P)         :: q
real(R8P)                :: w_fast, w_exact

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
q = vector_R8P(0.3_R8P, 0.3_R8P, 0.3_R8P)

w_fast  = cube%winding_number(point=q, beta=2._R8P)   ! Barnes-Hut, default
w_exact = cube%winding_number(point=q, beta=-1._R8P)  ! exact O(N) sum

print '(A,F12.9,A,F12.9)', 'fast = ', w_fast, '   exact = ', w_exact
endprogram ex_winding_number_exact
```

The two are identical to ~1e-12 absolute on small inputs (a 12-facet
cube has nothing to truncate); on the dragon mesh the gap is ~1e-3 at
β=2.0 and ~1e-5 at β=4.0.

## Cube reference image

![Cube fixture used in the examples](/pictures/advanced/winding-cube.png)

The cube is the simplest test case: 12 facets, three orthogonal
diagonal-split faces, axis-aligned. A query at the centroid `(0.5, 0.5,
0.5)` returns `w = 1.0` to floating-point precision; a query at a
corner returns `w = 0.125` (one octant of the full $4\pi$ solid angle,
which is what you should expect from this formulation at a vertex
singularity).

## Known limitations

- **β tuning is empirical**. β = 2.0 is fine for typical CAD STL with
  few thousand to ~100k facets. For dense meshes (~1M+ facets) where
  query throughput matters, β = 1.5 trades a few percent accuracy for
  ~2× speedup. For ground-truth tests, force exact (`beta ≤ 0`).
- **At surface vertices and edges**, `w` returns the local solid-angle
  fraction (e.g. 1/8 at a cube corner), not 0 or 1. This is **the
  correct mathematical answer** — the "inside vs. outside" question is
  ill-posed exactly on the surface. Threshold at 0.5 to get a clean
  binary classifier; query a small offset away from the surface if you
  need a definite answer at a known vertex.
- **Non-watertight inputs return values in (0, 1)** — that's the whole
  point, but if your downstream code assumed `w ∈ {0, 1}` it'll need
  updating. Threshold at 0.5 to recover the binary classifier.

## See also

- [`distance` and the `SIGN_*` constants](/guide/api/surface-stl-object#distance) —
  legacy sign algorithms (kept for back-compat; prefer winding number on
  dirty input).
- [Constants](/guide/api/constants) — `SIGN_PSEUDO_NORMAL` is FOSSIL's
  default for `distance`'s sign; winding number is its own scalar query
  with no constant needed.
- [`aabb_tree_object`](/guide/api/aabb-tree-object) — the acceleration
  structure that carries the per-node dipole moments used by the
  hierarchical traversal.

## References

- Jacobson, Kavan & Sorkine-Hornung, *Robust Inside-Outside Segmentation
  Using Generalized Winding Numbers*, SIGGRAPH 2013. The original
  algorithm — the formulation FOSSIL implements.
- Barill, Dickson, Schmidt, Levin & Jacobson, *Fast Winding Numbers for
  Soups and Clouds*, SIGGRAPH 2018. The Barnes-Hut hierarchical
  evaluation that makes large meshes tractable.
- Van Oosterom & Strackee, *The solid angle of a plane triangle*, IEEE
  TBME 30(2), 1983. The closed-form per-facet solid angle.
