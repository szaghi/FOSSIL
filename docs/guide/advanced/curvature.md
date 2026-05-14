---
title: Per-vertex curvature
---

# Per-vertex Gaussian + mean curvature

## What it does

Two TBPs that return scalar curvature fields, one value per vertex of
the surface:

- **`gaussian_curvature(K, status)`** — Gaussian curvature `K_i`
  computed from the **angle defect** at vertex `i`. `K > 0` on
  convex caps (sphere, ellipsoid), `K < 0` on saddles, `K = 0` on
  developable surfaces (cylinder, cone, flat plane).
- **`mean_curvature(H, status)`** — signed mean curvature `H_i` from
  the **mean-curvature normal** identity `H_i n_i = (1/2) M⁻¹ (L V)_i`,
  with sign disambiguated by projection onto the per-vertex pseudo-normal.
  `H > 0` on outward-convex regions (sphere surface, cube corners),
  `H < 0` on saddles / concave-outward regions, `H ≈ 0` on flat regions.

Both build on [cotangent Laplacian + barycentric mass](/guide/advanced/cotangent-laplacian).
The mean variant in particular is a one-liner on top of `(L, M)`:
`L V` is the per-vertex mean-curvature normal scaled by `2 M_ii`;
inverting the diagonal mass and halving gives `H n`.

### A correction worth surfacing

The `feat(curvature,surface): per-vertex Gaussian curvature` commit
message claimed that the **mean** variant required a sparse linear
solver. That was wrong. The mass matrix `M` is *diagonal* (barycentric
per §2.1), so `M⁻¹` is per-row scalar division — no factorisation
needed. The mean variant shipped without any external dependency in the
same release as Gaussian. The mean-curvature TBP docstring carries the
retraction inline; this page records it once and moves on.

## Pipeline

### Gaussian — angle defect

For each interior vertex `i`:

$$K_i = \frac{2\pi - \sum_{T \ni i} \theta_i^T}{A_i}$$

where `θ_i^T` is the angle at vertex `i` inside triangle `T`, the sum
runs over all triangles incident to `i`, and `A_i = (1/3) Σ_{T ∋ i} A_T`
is the same barycentric vertex-area used in §2.1's mass matrix.

For boundary vertices the formula is `(π − Σ θ) / A_i` — half the angle
defect, reflecting the fact that the surface doesn't wrap around `i` on
the open side. FOSSIL detects boundary vertices automatically via the
half-edge incidence pool.

The integrated angle defect `Σ_i K_i A_i = 4π χ` (Gauss–Bonnet) for any
closed orientable manifold of Euler characteristic `χ` — `χ = 2` for a
topological sphere, `χ = 0` for a torus. This is a standard regression
oracle for the implementation.

### Mean — mean-curvature normal

For each vertex `i`, compute the action of the cotangent Laplacian on
the coordinate field `V`:

$$(L V)_i = \sum_{j \in N(i)} \frac{\cot \alpha_{ij} + \cot \beta_{ij}}{2} \, (V_j - V_i)$$

then form

$$H_i \, \vec{n}_i = \frac{1}{2 A_i} (L V)_i$$

The magnitude is `H_i = ||H_i n_i||`; the sign is decided by projecting
onto the angle-weighted pseudo-normal computed by `analyze` (this is
the same per-vertex pseudo-normal used by `SIGN_PSEUDO_NORMAL` for
signed-distance queries). Outward-convex projection gives `H > 0`,
inward-concave gives `H < 0`.

Both computations reuse the §2.1 per-triangle cotangent and area
assembly. They run in a single pass over the facet array with O(N)
storage; no sparse matrix is materialised by either TBP — the assembly
is "on the fly" into the output vector.

## API

### Gaussian

```fortran
call surface%gaussian_curvature(K, status)
   class(surface_stl_object),  intent(in)            :: self
   real(R8P), allocatable,     intent(out)           :: K(:)
   integer(I4P),               intent(out), optional :: status   ! CURV_STATUS_*
```

Output:

- **`K(1:n_vertices)`** — per-vertex Gaussian curvature. Indexed by
  pool vertex id (`pool%get_vertex_count()` gives `n_vertices`). Boundary
  vertices use the half-defect formula; interior vertices the full
  formula.

### Mean

```fortran
call surface%mean_curvature(H, status)
   class(surface_stl_object),  intent(in)            :: self
   real(R8P), allocatable,     intent(out)           :: H(:)
   integer(I4P),               intent(out), optional :: status   ! CURV_STATUS_*
```

Output:

- **`H(1:n_vertices)`** — per-vertex signed mean curvature. Positive on
  outward-convex regions, negative on saddle / inward-concave, ≈ 0 on
  flat regions.

Shared status codes:

- `CURV_STATUS_OK` — operator built successfully.
- `CURV_STATUS_BAD_INPUT` — empty surface or uninitialised pool.
- `CURV_STATUS_DEGENERATE_TRIANGLE` — at least one triangle had area
  below tolerance and was skipped. The curvature field is still
  returned; degenerate vertices receive 0.

**Pre-condition.** `analyze` must have run (pseudo-normals populated).
Every public construction path on `surface_stl_object` ensures this —
`load_from_file`, `adopt_facets`, and `sanitize` all call `analyze`
internally. The `mean_curvature` sign-disambiguation step is the only
direct consumer of pseudo-normals among the curvature TBPs; Gaussian
does not depend on them.

## Examples

### Gaussian curvature on a sphere

```fortran
program ex_gaussian_curvature
use fossil
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P
implicit none

type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
real(R8P),          allocatable :: values(:,:,:), curv_k(:)
real(R8P)                       :: x, y, z, dx
integer(I4P), parameter         :: ng = 32_I4P
integer(I4P)                    :: ix, iy, iz, status

! Build a unit-sphere surface: sample the analytic SDF on a 32^3 grid
! spanning [-2, 2]^3, then marching-cubes the iso=0 level set.
allocate(values(ng, ng, ng))
dx = 4._R8P / real(ng - 1_I4P, R8P)
do iz = 1_I4P, ng
   z = -2._R8P + (iz - 1_I4P) * dx
   do iy = 1_I4P, ng
      y = -2._R8P + (iy - 1_I4P) * dx
      do ix = 1_I4P, ng
         x = -2._R8P + (ix - 1_I4P) * dx
         values(ix, iy, iz) = sqrt(x*x + y*y + z*z) - 1._R8P
      enddo
   enddo
enddo
call extract_isosurface(values=values, bmin=vector_R8P(-2._R8P, -2._R8P, -2._R8P), &
                        bmax=vector_R8P(2._R8P, 2._R8P, 2._R8P), iso=0._R8P,       &
                        surface=facets, status=status)
call sphere%adopt_facets(facets=facets)

call sphere%gaussian_curvature(K=curv_k, status=status)
print '(A,I0)',      'curv status = ', status
print '(A,I0)',      'n_vertices  = ', size(curv_k, kind=I4P)
print '(A,ES12.5)',  'mean K      = ', sum(curv_k) / real(size(curv_k, kind=I4P), R8P)
! A unit sphere has K = 1 everywhere; the discrete mean lands near 1
! (tessellation noise spreads the per-vertex values — see Known limitations).
endprogram ex_gaussian_curvature
```

### Mean curvature on a unit sphere

```fortran
program ex_mean_curvature_sphere
use fossil
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P
implicit none

type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
real(R8P),          allocatable :: values(:,:,:), h_curv(:)
real(R8P)                       :: x, y, z, dx
integer(I4P), parameter         :: ng = 32_I4P
integer(I4P)                    :: ix, iy, iz, n_pos, n_total, status

! Same unit-sphere construction as the Gaussian example above.
allocate(values(ng, ng, ng))
dx = 4._R8P / real(ng - 1_I4P, R8P)
do iz = 1_I4P, ng
   z = -2._R8P + (iz - 1_I4P) * dx
   do iy = 1_I4P, ng
      y = -2._R8P + (iy - 1_I4P) * dx
      do ix = 1_I4P, ng
         x = -2._R8P + (ix - 1_I4P) * dx
         values(ix, iy, iz) = sqrt(x*x + y*y + z*z) - 1._R8P
      enddo
   enddo
enddo
call extract_isosurface(values=values, bmin=vector_R8P(-2._R8P, -2._R8P, -2._R8P), &
                        bmax=vector_R8P(2._R8P, 2._R8P, 2._R8P), iso=0._R8P,       &
                        surface=facets, status=status)
call sphere%adopt_facets(facets=facets)

call sphere%mean_curvature(H=h_curv, status=status)
n_total = size(h_curv, kind=I4P)
n_pos   = count(h_curv > 0._R8P)
print '(A,ES12.5,A)', 'mean H = ', sum(h_curv) / real(n_total, R8P), &
                      '   (analytic value for a unit sphere = 1.0)'
print '(A,I0,A,I0,A)', 'sign: ', n_pos, ' / ', n_total, ' vertices have H > 0'
! Expected: ~100% positive — the sphere surface is convex outward everywhere.
endprogram ex_mean_curvature_sphere
```

### Combining both for principal curvature follow-up

The two scalars together carry the full second-fundamental-form
trace + determinant:

```fortran
! Both shipped today:
call surface%gaussian_curvature(K=K, status=stat_K)
call surface%mean_curvature(H=H, status=stat_H)

! Principal curvatures (not yet a TBP — small follow-up):
!     kappa1 = H + sqrt(max(0, H^2 - K))
!     kappa2 = H - sqrt(max(0, H^2 - K))
! Pointwise on every vertex. Discriminant clamp guards floating-point
! noise on near-umbilical (kappa1 = kappa2) points.
```

A `principal_curvatures` TBP is a documented follow-up — small wrapper,
no new infrastructure. Tracked under [issue
#18](https://github.com/szaghi/FOSSIL/issues/18).

## Visual reference

Gaussian curvature `K` on a unit sphere (built by marching cubes, then
isotropically remeshed), coloured through a diverging map:

![Gaussian curvature K on a remeshed unit sphere](/pictures/advanced/curvature-sphere-gaussian.png)

Signed mean curvature `H` on the same surface:

![Signed mean curvature H on a remeshed unit sphere](/pictures/advanced/curvature-sphere-mean.png)

**Read these figures honestly.** The bulk of each sphere reads as a
near-uniform colour — that uniform field *is* the correct answer:
`K ≈ 1` and `H ≈ 1` everywhere on a unit sphere. The speckled
hot-spot clusters are **discrete-curvature artefacts at
tessellation-degenerate vertices**, exactly the coarse-mesh noise
described in [Known limitations](#known-limitations) below. They are
not a rendering bug and they are not a bug in the curvature TBPs —
they are what the discrete angle-defect and mean-curvature-normal
formulae genuinely produce on the sliver triangles that survive even
after a remesh of a marching-cubes sphere. The colour bar is clamped
to the 2nd/98th percentile so the bulk signal stays visible despite
the outliers. On a clean analytic triangulation the field is visibly
uniform; on real CAD STL, run [§1.7 isotropic
remeshing](/guide/advanced/isotropic-remesh) first if the curvature
field looks this noisy. The fixtures and renders are reproducible via
`scripts/make_doc_images.sh`.

## Known limitations

- **Discrete curvature is noisy on coarse triangulations.** Both
  formulae carry an `O(h)` discretisation error where `h` is the local
  edge length; on very-coarse meshes (cube with 12 facets, low-res MC)
  the values are dominated by tessellation artefact rather than smooth
  curvature. This is why the §2.4 tests target sphere + bunny rather
  than the cube — sphere has well-behaved curvature, bunny has enough
  vertices to swamp tessellation noise. Use [§1.7 isotropic
  remeshing](/guide/advanced/isotropic-remesh) or [§1.3
  decimation](/guide/advanced/decimate) first if the input is
  pathological.
- **Mean-curvature sign relies on the pseudo-normal.** If the surface
  has inconsistent normal orientation (one patch flipped), the
  per-vertex pseudo-normal may point inward at some vertices and the
  sign of `H` there will flip. Run `sanitize_normals` (or full
  `sanitize`) before `mean_curvature` if the input is dirty. The
  unsigned magnitude `|H|` is robust regardless.
- **Principal curvatures and directions** are *not yet* a TBP — only
  the trace (`H`) and determinant-shadow (`K`) are returned. Computing
  `(κ₁, κ₂, d₁, d₂)` requires either (a) the algebraic combination
  shown in the example above for the eigenvalues, or (b) Rusinkiewicz
  2004 per-vertex 3×3 covariance + `dsyev` LAPACK eigen-decomposition
  for the directions. The follow-up TBP is small but not yet shipped.
- **Discrete curvature is per-*vertex*, not per-*facet*.** If you want
  a per-facet curvature scalar (e.g. for face-based segmentation), the
  standard recipe is to average the three incident vertex values —
  the curvature TBPs do not provide a per-facet helper.

## See also

- [cotangent Laplacian + barycentric mass](/guide/advanced/cotangent-laplacian)
  — the prerequisite operator. Mean curvature is a one-liner on top of
  `(L, M)`.
- [mesh smoothing](/guide/advanced/smoothing) — the shipped
  smoothing variants reduce high-frequency surface noise; running
  `mean_curvature` after a few Taubin passes gives noticeably cleaner
  `H` fields without volume loss.
- [`sanitize_normals`](/guide/api/surface-stl-object#sanitize_normals)
  — ensures the pseudo-normals used by mean-curvature sign
  disambiguation are consistently outward-oriented.
- [Constants](/guide/api/constants) — `CURV_STATUS_*`.

## References

- Meyer, Desbrun, Schröder & Barr, *Discrete differential-geometry
  operators for triangulated 2-manifolds*, VisMath 2003. The reference
  formulation for both `K_i` (angle defect) and `H_i n_i`
  (mean-curvature normal), including the mixed-area mass discussion
  that informs the §2.1 follow-up.
- Rusinkiewicz, *Estimating curvatures and their derivatives on
  triangle meshes*, 3DPVT 2004. The per-vertex 3×3 covariance method
  for principal curvatures and directions — the documented follow-up.
- Crane, *Discrete Differential Geometry: An Applied Introduction*,
  SIGGRAPH course notes 2013–2023. The pedagogical reference covering
  both Gauss–Bonnet on triangle meshes and the `H n = (1/2) M⁻¹ L V`
  identity used by the FOSSIL mean-curvature TBP.
- do Carmo, *Differential Geometry of Curves and Surfaces*,
  Prentice-Hall 1976. The smooth-surface theory the discrete formulae
  approximate — useful when you need to remember why
  `H = (κ₁ + κ₂)/2` and `K = κ₁ κ₂`.
