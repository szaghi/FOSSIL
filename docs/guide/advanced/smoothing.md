---
title: Mesh smoothing
---

# Mesh smoothing (explicit Laplacian + Taubin)

## What it does

Attenuates high-frequency surface noise on a triangle mesh by moving
each vertex toward the average position of its mesh neighbours. Two
flavours are exposed today behind the same TBP, picked by a `method`
argument:

- **Explicit Laplacian** (`SMOOTH_METHOD_EXPLICIT`) — `iterations`
  single-step applications of `V_new = V + λ ΔV`, where `ΔV` is the
  uniform-weight Laplacian of the vertex field. **Shrinks the volume.**
  Cheap and useful as a building block, but the wrong primitive for
  production denoising of a closed solid.
- **Taubin λ\|μ** (`SMOOTH_METHOD_TAUBIN`, default) — alternates a
  positive `λ` smoothing step with a negative `μ` counter-shrinking step.
  Designed by Taubin (1995) to have a *band-pass* frequency response:
  high-frequency surface variation gets attenuated while the
  low-frequency volume stays nearly invariant. **The right choice for
  CFD-grade STL denoising.**

A third variant — **implicit Laplacian** (`(M − τL) V_new = M V`,
unconditionally stable at any step size, requires sparse SPD
factorisation) — is deferred until a sparse linear solver is bound to
FOSSIL. The public API is forward-compatible with it: when the solver
lands, the new variant ships as `SMOOTH_METHOD_IMPLICIT` and no
caller-side change is needed.

## Why uniform weights, not cotangent

A natural reading of the [cotangent Laplacian](/guide/advanced/cotangent-laplacian)
page and the Desbrun–Meyer–Schröder–Barr 1999 paper is that the *right*
Laplacian for explicit smoothing is the area-normalised cotangent
`M⁻¹ L`. FOSSIL ships the uniform-weight Laplacian instead,
deliberately. The reason:

- `M_ii` is small on dense triangulations (it's the *third of the
  incident triangle area* per vertex). On a 100k-facet bunny,
  `min(M_ii)` is of order `10⁻⁶`.
- For *explicit* time-stepping, the stability bound is roughly
  `λ ≤ 2 / λ_max(M⁻¹ L)`. With small `M_ii`, `λ_max(M⁻¹ L)` is huge, so
  any "meaningful" `λ` blows the iteration up.
- The uniform-weight Laplacian (`(1/|N(i)|) Σ_{j ∈ N(i)} (V_j − V_i)`)
  is unconditionally stable for any `λ ≤ 1`. It loses the discrete-
  geometry consistency of `M⁻¹ L` but is the right operator for
  *iterative noise filtering*, which is what `surface%smooth` is for.

If you need the area-normalised mean-curvature flow (e.g. for
implicit-fairing experiments), build `(L, M)` via the
[cotangent Laplacian](/guide/advanced/cotangent-laplacian) TBP and
time-step it yourself once the solver lands. The shipped TBP is the
production-grade denoiser.

## Pipeline

Per smoothing step (one for explicit, two — `λ` then `μ` — for Taubin):

1. **Per-triangle edge accumulation.** For each triangle, each of its
   three edges contributes to the `(V_j − V_i)` accumulator of both
   endpoints and to both endpoints' neighbour-count accumulator. No
   adjacency list is materialised — every edge is implicitly visited
   twice (once per incident triangle), which double-counts both the
   accumulator and the count *symmetrically*, so the per-vertex
   division cancels the factor and the result is the correct uniform-
   weight Laplacian.

2. **Per-vertex update.**

   $$V_i^{\text{new}} = V_i + \lambda \, \frac{1}{|N(i)|} \sum_{j \in N(i)} (V_j - V_i)$$

3. **Taubin band-pass guard.** For `SMOOTH_METHOD_TAUBIN`, the second
   step uses `μ < −λ`. This is the band-pass condition from Taubin
   1995: without it, the iteration is still a low-pass filter (it just
   shrinks more slowly), which is almost certainly not what the caller
   wants. The TBP **rejects with `SMOOTH_STATUS_BAD_INPUT`** when
   `μ ≥ −λ` rather than silently doing the wrong thing.

4. **Re-adopt after the last iteration.** Smoothing mutates the
   per-facet vertex caches in place. After return the TBP calls
   `adopt_facets` so the vertex pool, connectivity, AABB tree, and
   pseudo-normals are coherent with the smoothed positions. No further
   user action is needed before calling `distance`, `winding_number`,
   `mean_curvature`, etc.

## API

```fortran
call surface%smooth(method, lambda, mu, iterations, status)
   class(surface_stl_object), intent(inout)         :: self
   integer(I4P),              intent(in),  optional :: method        ! SMOOTH_METHOD_*
   real(R8P),                 intent(in),  optional :: lambda
   real(R8P),                 intent(in),  optional :: mu            ! TAUBIN only
   integer(I4P),              intent(in),  optional :: iterations
   integer(I4P),              intent(out), optional :: status        ! SMOOTH_STATUS_*
```

Defaults:

| Argument     | Default                       |
|--------------|-------------------------------|
| `method`     | `SMOOTH_METHOD_TAUBIN`        |
| `lambda`     | `SMOOTH_DEFAULT_LAMBDA = 0.5` |
| `mu`         | `SMOOTH_DEFAULT_MU = −0.53`   |
| `iterations` | `SMOOTH_DEFAULT_ITERATIONS = 5` |

Notes on `iterations`:

- For `SMOOTH_METHOD_EXPLICIT`, `iterations` is the number of `λ`-steps
  applied. `mu` is silently ignored.
- For `SMOOTH_METHOD_TAUBIN`, `iterations` is the number of *complete
  `(λ, μ)` pairs*. The default value of 5 thus produces 10 Laplacian
  applications.

Status codes:

- `SMOOTH_STATUS_OK` — smoothing completed normally.
- `SMOOTH_STATUS_BAD_INPUT` — empty surface, `iterations ≤ 0`, unknown
  `method`, or Taubin called with `μ ≥ −λ`.
- `SMOOTH_STATUS_DEGENERATE` — at least one triangle had area below
  tolerance during the per-step Laplacian assembly; the iteration
  proceeded on the rest. The surface is still re-adopted at the end.

## Examples

### Taubin denoising of a scanned bunny

```fortran
program ex_smooth_taubin
use fossil
use penf, only : I4P, R8P
implicit none

type(surface_stl_object) :: bunny
real(R8P)                :: v_before, v_after, a_before, a_after
integer(I4P)             :: status

call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
v_before = bunny%get_volume()
a_before = bunny%get_area()

! 5 Taubin (lambda, mu) pairs = 10 Laplacian applies. Defaults: lambda=0.5, mu=-0.53.
call bunny%smooth(method=SMOOTH_METHOD_TAUBIN, status=status)
v_after = bunny%get_volume()
a_after = bunny%get_area()

print '(A,ES12.5,A,ES12.5)', 'volume before / after = ', v_before, ' / ', v_after
print '(A,ES12.5,A,ES12.5)', 'area   before / after = ', a_before, ' / ', a_after
print '(A,F6.3,A)',           'volume drift = ', 100._R8P * abs(v_after - v_before) / v_before, ' %'
print '(A,F6.3,A)',           'area drop    = ', 100._R8P * (a_before - a_after) / a_before, ' %'
! Expected: volume drift < 1 %, area drop ~ 0.5-1 % (visible smoothing).
endprogram ex_smooth_taubin
```

### Explicit shrink for a building block

```fortran
program ex_smooth_explicit
use fossil
use penf, only : I4P, R8P
implicit none

type(surface_stl_object) :: noisy
integer(I4P)             :: status

call noisy%load_from_file(file_name='noisy_input.stl', guess_format=.true.)

! One pass at half-strength, only as part of a larger pipeline
! (e.g. before decimation, where local shrinkage doesn't matter).
call noisy%smooth(method=SMOOTH_METHOD_EXPLICIT, lambda=0.25_R8P, &
                  iterations=1_I4P, status=status)
call noisy%save_into_file(file_name='preprocessed.stl')
endprogram ex_smooth_explicit
```

### Custom Taubin parameters

```fortran
! Stronger smoothing, more iterations:
call surface%smooth(method=SMOOTH_METHOD_TAUBIN,    &
                    lambda=0.7_R8P, mu=-0.75_R8P,   &  ! mu < -lambda — band-pass holds
                    iterations=10_I4P, status=status)

! Catches the band-pass violation:
call surface%smooth(method=SMOOTH_METHOD_TAUBIN,    &
                    lambda=0.5_R8P, mu=-0.4_R8P,    &  ! mu > -lambda — REJECTED
                    iterations=5_I4P, status=status)
! status == SMOOTH_STATUS_BAD_INPUT, no mutation
```

## Visual reference

Stanford bunny with deterministic per-vertex Gaussian noise added at
~1 % of the bounding-box diagonal — the "input" panel:

![Noisy bunny — input to Taubin smoothing](/pictures/advanced/smoothing-bunny-noisy.png)

After 5 Taubin `(λ, μ)` pairs at the defaults (`λ = 0.5`, `μ = -0.53`):

![Bunny after 5 Taubin pairs](/pictures/advanced/smoothing-bunny-taubin.png)

The high-frequency surface "fuzz" is gone — the silhouette and the
low-frequency shape (ears, haunch, the overall volume) are preserved.
This is the band-pass response Taubin is designed for: it attenuates
the noise band without the global shrinkage that plain explicit
Laplacian smoothing would cause. The fixture and both renders are
reproducible via `scripts/make_doc_images.sh` (the noise seed is
fixed, so the figures are byte-stable across rebuilds).

## Known limitations

- **Implicit variant not yet shipped.** The third method from the §2.3
  spec — `(M − τL) V_new = M V`, unconditionally stable at any `τ`,
  requires a sparse SPD solver. Deferred until a sparse linear solver
  (SuiteSparse or alternative) is bound to FOSSIL. The public API is
  forward-compatible: when the solver arrives, the new variant ships
  as `SMOOTH_METHOD_IMPLICIT` and no caller-side change is needed.
- **Sharp-feature blunting.** Both explicit and Taubin treat every
  vertex uniformly. Vertices on a sharp edge (e.g. cube corners) get
  pulled inward by neighbour averaging and the feature blurs. If
  feature preservation matters, run [§1.7 isotropic
  remeshing](/guide/advanced/isotropic-remesh) (which has a built-in
  dihedral-angle feature lock) instead, or restrict smoothing to a
  feature-aware vertex mask externally before calling `smooth`.
- **Boundary vertices on open meshes** lose neighbours on the boundary
  side, so the uniform-weight Laplacian pulls them inward toward the
  interior centroid. This is rarely an issue for FOSSIL's target
  workload (closed solids from CAD STL) but matters if you smooth an
  open shell.
- **Degenerate-triangle skip.** Triangles below the area tolerance are
  skipped during Laplacian assembly. Their incident vertices simply
  receive less averaging that step; they aren't *removed*. Run [§1.3
  decimation](/guide/advanced/decimate) or `sanitize` before `smooth`
  if degenerate triangles dominate the input.
- **The cube test of `fossil_test_smoothing`** confirms only that the
  degenerate-input case (corner-only mesh, no flat interior) returns
  `SMOOTH_STATUS_OK` without crashing — the cube geometry collapses
  under uniform-weight smoothing because every corner is pulled toward
  the face centroid. Use Taubin's parameters for production inputs;
  the cube is a "doesn't crash on a pathological input" guard, not a
  recommended use case.

## See also

- [cotangent Laplacian](/guide/advanced/cotangent-laplacian) — the
  prerequisite that the *future* implicit smoothing variant depends
  on. The shipped explicit + Taubin variants use a different Laplacian
  (uniform weights, not cotangent) for stability — see "Why uniform
  weights, not cotangent" above.
- [per-vertex curvature](/guide/advanced/curvature) — running a
  few Taubin passes before `mean_curvature` gives noticeably cleaner
  curvature fields without volume loss.
- [isotropic remeshing](/guide/advanced/isotropic-remesh) — has
  built-in feature-preserving tangential relaxation. Use it instead of
  smoothing when sharp edges matter.
- [Constants](/guide/api/constants) — `SMOOTH_STATUS_*`,
  `SMOOTH_METHOD_*`, `SMOOTH_DEFAULT_LAMBDA`, `SMOOTH_DEFAULT_MU`,
  `SMOOTH_DEFAULT_ITERATIONS`.

## References

- Taubin, *A signal-processing approach to fair surface design*,
  SIGGRAPH 1995. The reference paper for the `λ\|μ` non-shrinking
  alternation that FOSSIL ships as the default.
- Desbrun, Meyer, Schröder & Barr, *Implicit fairing of irregular
  meshes using diffusion and curvature flow*, SIGGRAPH 1999. The
  reference for the *implicit* variant deferred until a sparse solver
  is bound; also explains why area-normalised cotangent smoothing is
  unstable on dense meshes for explicit time-stepping (Section 3.2).
- Botsch, Kobbelt, Pauly, Alliez & Lévy, *Polygon Mesh Processing*,
  A K Peters 2010, Chapter 4. Book-length treatment of mesh smoothing
  including the feature-preservation extensions that motivate falling
  back to [isotropic remeshing](/guide/advanced/isotropic-remesh)
  for sharp-feature inputs.
