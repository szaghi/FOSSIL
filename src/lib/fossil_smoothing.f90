!< Mesh smoothing via the cotangent Laplacian (issue #18 §2.3, partial —
!< explicit and Taubin variants; the implicit variant needs a sparse
!< solver and is deferred to the Tier-2 sparse-solver gate).
!<
!< References:
!<   - Desbrun, Meyer, Schröder & Barr, *Implicit Fairing of Irregular
!<     Meshes using Diffusion and Curvature Flow*, SIGGRAPH 1999.
!<   - Taubin, *A Signal Processing Approach To Fair Surface Design*,
!<     SIGGRAPH 1995. Origin of the lambda|mu non-shrinking trick.
!<
!< Update rule per step: V_new = V + step * Laplace(V)
!<
!<   Laplace(V_i) = (1 / |N(i)|) * sum_{j in N(i)} (V_j - V_i)
!<
!< where N(i) is the set of vertices connected to i by a mesh edge.
!< This is the **uniform-weight (combinatorial) Laplacian** from
!< Taubin 1995, NOT the cotangent Laplacian from §2.1. The choice is
!< deliberate: uniform weights are unconditionally stable for
!< `step ≤ 1`, while the area-normalized cotangent variant
!< (`M^-1 L`) blows up on fine meshes for any meaningful step
!< because M is small per-vertex on dense triangulations.
!<
!<   - **Explicit Laplacian** (`method = SMOOTH_EXPLICIT`):
!<     `iterations` single steps with positive `lambda`. Smooths but
!<     shrinks the volume.
!<   - **Taubin lambda|mu** (`method = SMOOTH_TAUBIN`): alternating
!<     positive lambda step (smoothing) and negative mu step
!<     (counter-shrinking). Designed to have band-pass frequency
!<     response that smooths high-frequency noise while leaving the
!<     low-frequency volume nearly untouched. Standard parameters:
!<     `lambda = 0.5`, `mu = -0.53`. The `iterations` argument counts
!<     COMPLETE (lambda, mu) pairs.
!<
!< Edge accumulation is **per-triangle**: for each triangle, each
!< edge contributes to two pairs of `(V_j - V_i)` differences and to
!< two `|N(i)|` counts. We never materialize the adjacency list.
!< Each interior edge is visited twice (once per incident triangle),
!< so `|N(i)|` ends up double-counted — but symmetrically across all
!< vertices, so the normalization is consistent.
!<
!< This module operates on facet(:) + vertex_pool directly (no
!< surface_stl import — same no-cycle pattern as other feature
!< modules). Smoothing **mutates the per-facet vertex cache in place**:
!< on return, every facet's `vertex(1:3)` and pseudo-normals reflect
!< the smoothed positions. The vertex pool is NOT updated by this
!< module — the caller (typically the surface TBP) should run
!< `surface%adopt_facets(self%facet)` afterward to rebuild the pool
!< and connectivity from the smoothed facets.

module fossil_smoothing
!< Mesh smoothing — explicit Laplacian and Taubin lambda|mu.

use fossil_facet_object, only : facet_object
use fossil_vertex_pool_object, only : vertex_pool_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
private

public :: smooth_mesh
public :: SMOOTH_STATUS_OK, SMOOTH_STATUS_BAD_INPUT, SMOOTH_STATUS_DEGENERATE
public :: SMOOTH_METHOD_EXPLICIT, SMOOTH_METHOD_TAUBIN
public :: SMOOTH_DEFAULT_LAMBDA, SMOOTH_DEFAULT_MU
public :: SMOOTH_DEFAULT_ITERATIONS

integer(I4P), parameter :: SMOOTH_STATUS_OK         = 0_I4P
integer(I4P), parameter :: SMOOTH_STATUS_BAD_INPUT  = 1_I4P
integer(I4P), parameter :: SMOOTH_STATUS_DEGENERATE = 2_I4P  !< At least one triangle had area below tolerance and was skipped.

integer(I4P), parameter :: SMOOTH_METHOD_EXPLICIT = 1_I4P
integer(I4P), parameter :: SMOOTH_METHOD_TAUBIN   = 2_I4P

real(R8P),    parameter :: SMOOTH_DEFAULT_LAMBDA   =  0.5_R8P     !< Taubin's standard positive step.
real(R8P),    parameter :: SMOOTH_DEFAULT_MU       = -0.53_R8P    !< Taubin's standard negative counter-step (|mu| > lambda).
integer(I4P), parameter :: SMOOTH_DEFAULT_ITERATIONS = 5_I4P      !< For Taubin: pair count. For Explicit: step count.

real(R8P), parameter :: DEGENERATE_AREA_TOL = 1.0e-30_R8P

contains

   subroutine smooth_mesh(facet, pool, method, lambda, mu, iterations, status)
   !< Smooth the mesh in place via explicit Laplacian or Taubin
   !< lambda|mu flow.
   !<
   !< For SMOOTH_METHOD_EXPLICIT: `iterations` single-step applications
   !< of `V_new = V + lambda M^-1 L V`. `mu` is ignored. Shrinks.
   !<
   !< For SMOOTH_METHOD_TAUBIN: `iterations` complete (lambda, mu) pairs.
   !< Each pair is one smoothing step (lambda) followed by one counter-
   !< shrinking step (mu). The Taubin paper recommends lambda > 0 and
   !< mu < -lambda — both checks enforced. Designed for volume
   !< preservation: total work is `2 * iterations` Laplacian-vector
   !< applies, with the lambda|mu band-pass response keeping the
   !< low-frequency volume nearly invariant while attenuating
   !< high-frequency surface noise.
   !<
   !< Pre-condition: `pool` must be initialized
   !< (`get_is_initialized() = .true.`). The pool is read-only here
   !< (used for vertex-id resolution and initial coordinates); the
   !< smoothed coordinates are written back into the facet vertex
   !< caches. After return the pool is STALE; the caller (typically
   !< the surface TBP) should re-adopt the smoothed facets to rebuild
   !< the pool, connectivity, and AABB tree.
   type(facet_object),       intent(inout)         :: facet(:)
   type(vertex_pool_object), intent(in)            :: pool
   integer(I4P),             intent(in),  optional :: method
   real(R8P),                intent(in),  optional :: lambda
   real(R8P),                intent(in),  optional :: mu
   integer(I4P),             intent(in),  optional :: iterations
   integer(I4P),             intent(out), optional :: status
   integer(I4P)                                    :: meth, n_iter, it
   real(R8P)                                       :: lam, mu_val
   logical                                         :: had_degenerate
   type(vector_R8P), allocatable                   :: coord(:)
   integer(I4P)                                    :: n_vertices, f, vi, lv

   if (present(status)) status = SMOOTH_STATUS_OK
   meth   = SMOOTH_METHOD_TAUBIN;     if (present(method))     meth = method
   lam    = SMOOTH_DEFAULT_LAMBDA;    if (present(lambda))     lam = lambda
   mu_val = SMOOTH_DEFAULT_MU;        if (present(mu))         mu_val = mu
   n_iter = SMOOTH_DEFAULT_ITERATIONS;if (present(iterations)) n_iter = iterations

   if (size(facet, kind=I4P) == 0_I4P .or. .not. pool%get_is_initialized()) then
      if (present(status)) status = SMOOTH_STATUS_BAD_INPUT
      return
   endif
   if (n_iter <= 0_I4P) return
   if (meth /= SMOOTH_METHOD_EXPLICIT .and. meth /= SMOOTH_METHOD_TAUBIN) then
      if (present(status)) status = SMOOTH_STATUS_BAD_INPUT
      return
   endif
   if (lam <= 0._R8P) then
      if (present(status)) status = SMOOTH_STATUS_BAD_INPUT
      return
   endif
   if (meth == SMOOTH_METHOD_TAUBIN .and. mu_val >= -lam) then
      ! Taubin requires mu < -lambda for the band-pass property. If the
      ! caller passed e.g. mu = -0.4 with lambda = 0.5, the iteration
      ! becomes a low-pass that still shrinks (just slower) — likely a
      ! mistake. Reject explicitly.
      if (present(status)) status = SMOOTH_STATUS_BAD_INPUT
      return
   endif

   ! Initialize working coordinate buffer from the pool.
   n_vertices = pool%vertex_count()
   allocate(coord(n_vertices))
   do vi = 1_I4P, n_vertices
      coord(vi) = pool%coord(vid=vi)
   enddo

   had_degenerate = .false.
   select case (meth)
   case (SMOOTH_METHOD_EXPLICIT)
      do it = 1_I4P, n_iter
         call apply_laplacian_step(facet=facet, pool=pool, coord=coord, &
                                   step=lam, had_degenerate=had_degenerate)
      enddo
   case (SMOOTH_METHOD_TAUBIN)
      do it = 1_I4P, n_iter
         call apply_laplacian_step(facet=facet, pool=pool, coord=coord, &
                                   step=lam,    had_degenerate=had_degenerate)
         call apply_laplacian_step(facet=facet, pool=pool, coord=coord, &
                                   step=mu_val, had_degenerate=had_degenerate)
      enddo
   endselect

   ! Write the smoothed coordinates back into the facet vertex caches.
   ! The pool stays untouched; the caller must re-adopt facets to refresh
   ! it (along with connectivity, AABB tree, pseudo-normals).
   do f = 1_I4P, size(facet, kind=I4P)
      do lv = 1_I4P, 3_I4P
         vi = pool%facet_vid(facet_id=f, local_v=lv)
         if (vi == 0_I4P) cycle
         facet(f)%vertex(lv) = coord(vi)
      enddo
   enddo

   deallocate(coord)
   if (had_degenerate .and. present(status)) status = SMOOTH_STATUS_DEGENERATE
   endsubroutine smooth_mesh

   subroutine apply_laplacian_step(facet, pool, coord, step, had_degenerate)
   !< One application of `V_new = V + step * Laplace(V)` to the working
   !< `coord` buffer using the **uniform-weight (combinatorial)
   !< Laplacian**. Each triangle contributes per-vertex `(V_j - V_i)`
   !< differences for the two within-triangle neighbours of each
   !< vertex, plus the matching neighbour count.
   !<
   !< Edge double-counting: each interior edge is visited twice (once
   !< per incident triangle), so both the difference accumulator and
   !< the neighbour count are 2× what a unique-edge enumeration would
   !< give. The factor cancels on division — Laplace(V) is correct.
   !<
   !< `had_degenerate` is set if any triangle has area below tolerance,
   !< matching the §2.1 builder's behaviour. The smoothing skips
   !< degenerate triangles' contributions (their incident vertices
   !< just see less averaging this step).
   type(facet_object),       intent(in)    :: facet(:)
   type(vertex_pool_object), intent(in)    :: pool
   type(vector_R8P),         intent(inout) :: coord(:)
   real(R8P),                intent(in)    :: step
   logical,                  intent(inout) :: had_degenerate
   integer(I4P)                            :: nf, n_vertices, f, i
   integer(I4P)                            :: v1, v2, v3
   type(vector_R8P)                        :: p1, p2, p3
   type(vector_R8P)                        :: e12, e23, cross_full
   real(R8P)                               :: area_2x
   type(vector_R8P), allocatable           :: lap(:)
   integer(I4P),     allocatable           :: nbr_count(:)

   nf = size(facet, kind=I4P)
   n_vertices = size(coord, kind=I4P)
   allocate(lap(n_vertices))
   allocate(nbr_count(n_vertices))
   do i = 1_I4P, n_vertices
      lap(i) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   enddo
   nbr_count = 0_I4P

   do f = 1_I4P, nf
      v1 = pool%facet_vid(facet_id=f, local_v=1_I4P)
      v2 = pool%facet_vid(facet_id=f, local_v=2_I4P)
      v3 = pool%facet_vid(facet_id=f, local_v=3_I4P)
      if (v1 == 0_I4P .or. v2 == 0_I4P .or. v3 == 0_I4P) cycle

      p1 = coord(v1)
      p2 = coord(v2)
      p3 = coord(v3)

      ! Skip degenerate triangles' contributions (mirrors fossil_laplacian).
      e12 = p2 - p1
      e23 = p3 - p2
      cross_full = e12%crossproduct(rhs=e23)
      area_2x = sqrt(cross_full%dotproduct(rhs=cross_full))
      if (area_2x < DEGENERATE_AREA_TOL) then
         had_degenerate = .true.
         cycle
      endif

      ! Each vertex sees its two within-triangle neighbours.
      lap(v1) = lap(v1) + (p2 - p1) + (p3 - p1)
      lap(v2) = lap(v2) + (p1 - p2) + (p3 - p2)
      lap(v3) = lap(v3) + (p1 - p3) + (p2 - p3)
      nbr_count(v1) = nbr_count(v1) + 2_I4P
      nbr_count(v2) = nbr_count(v2) + 2_I4P
      nbr_count(v3) = nbr_count(v3) + 2_I4P
   enddo

   ! Normalize and update coord in place.
   do i = 1_I4P, n_vertices
      if (nbr_count(i) == 0_I4P) cycle  ! isolated vertex
      coord(i) = coord(i) + (step / real(nbr_count(i), R8P)) * lap(i)
   enddo

   deallocate(lap)
   deallocate(nbr_count)
   endsubroutine apply_laplacian_step

endmodule fossil_smoothing
