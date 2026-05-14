!< Per-vertex Gaussian curvature via the angle-defect formula
!< (issue #18 §2.4, partial — Gaussian variant only; the mean-curvature
!< variant via L · V / (2 M) needs a sparse-matrix solve and is deferred
!< to the Tier-2 sparse-solver gate).
!<
!< Definition:
!<
!<   K_i = (2*pi - sum_{T containing i} theta_T(i))  / A_i      interior
!<       = (pi   - sum_{T containing i} theta_T(i))  / A_i      boundary
!<
!< where theta_T(i) is the angle at vertex i within triangle T, and
!< A_i = (1/3) * sum_{T containing i} area(T) is the barycentric area.
!<
!< Boundary detection uses `facet%fcon_edge`: a vertex is on the
!< boundary iff at least one edge incident on it has no opposite-side
!< facet. For closed solids (cube, sphere) no vertex is on the
!< boundary.
!<
!< Sign convention:
!<   - K > 0 at convex vertices (cube corners, sphere surface)
!<   - K < 0 at saddle / hyperbolic vertices
!<   - K ~= 0 at locally-flat regions
!<
!< Gauss-Bonnet check: for a closed orientable surface of genus g,
!<   sum_i K_i * A_i = 2 * pi * (2 - 2g)
!< In particular, sum = 4*pi for a topological sphere (cube, sphere,
!< bunny) and 0 for a torus. This is the load-bearing global test of
!< the discrete formulation.
!<
!< References:
!<   - Meyer, Desbrun, Schröder & Barr, *Discrete Differential-Geometry
!<     Operators for Triangulated 2-Manifolds*, VisMath 2002. §4.3 for
!<     the angle-defect Gaussian curvature.
!<   - Crane, *Discrete Differential Geometry: An Applied Introduction*,
!<     SIGGRAPH course 2013. Chapter 5.
!<
!< Module operates on facet(:) + vertex_pool directly (no surface_stl
!< import — same no-cycle pattern as other feature modules).

module fossil_curvature
!< Per-vertex Gaussian curvature via angle defect.

use fossil_facet_object, only : facet_object
use fossil_vertex_pool_object, only : vertex_pool_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
private

public :: build_gaussian_curvature
public :: build_mean_curvature
public :: CURV_STATUS_OK, CURV_STATUS_BAD_INPUT, CURV_STATUS_DEGENERATE_TRIANGLE

integer(I4P), parameter :: CURV_STATUS_OK                   = 0_I4P
integer(I4P), parameter :: CURV_STATUS_BAD_INPUT            = 1_I4P  !< Empty facet array or pool not initialized.
integer(I4P), parameter :: CURV_STATUS_DEGENERATE_TRIANGLE  = 2_I4P  !< At least one triangle has area below tolerance; affected vertices skip that triangle.

real(R8P), parameter :: DEGENERATE_AREA_TOL = 1.0e-30_R8P
real(R8P), parameter :: PI = 4._R8P * atan(1._R8P)

contains

   subroutine build_gaussian_curvature(facet, pool, K, status)
   !< Compute per-vertex Gaussian curvature `K(1:n_vertices)` via the
   !< angle-defect formula.
   !<
   !< Vertices whose accumulated barycentric area is zero (no incident
   !< triangle survived the degenerate filter) get `K = 0` — defensible
   !< default since the angle defect is undefined for an isolated vertex.
   !<
   !< Boundary vertices use the `(pi - sum angles)` numerator; interior
   !< vertices use `(2*pi - sum angles)`. Boundary detection via
   !< `facet%fcon_edge == 0` on either incident edge.
   type(facet_object),       intent(in)            :: facet(:)
   type(vertex_pool_object), intent(in)            :: pool
   real(R8P), allocatable,   intent(out)           :: K(:)
   integer(I4P),             intent(out), optional :: status
   integer(I4P)                                    :: n_vertices, nf, f
   integer(I4P)                                    :: v1, v2, v3
   type(vector_R8P)                                :: p1, p2, p3
   type(vector_R8P)                                :: e12, e23, e31
   type(vector_R8P)                                :: cross_full
   real(R8P)                                       :: area_2x, area_tri
   real(R8P)                                       :: ang1, ang2, ang3
   real(R8P), allocatable                          :: angle_sum(:)  !< Per-vertex sum of incident triangle angles.
   real(R8P), allocatable                          :: area_third(:) !< Per-vertex 1/3 * incident-triangle-area sum.
   logical, allocatable                            :: on_boundary(:)
   integer(I4P)                                    :: i, e1
   logical                                         :: had_degenerate

   if (present(status)) status = CURV_STATUS_OK
   nf = size(facet, kind=I4P)
   if (nf == 0_I4P .or. .not. pool%get_is_initialized()) then
      allocate(K(0))
      if (present(status)) status = CURV_STATUS_BAD_INPUT
      return
   endif

   n_vertices = pool%vertex_count()
   allocate(K(n_vertices))
   allocate(angle_sum(n_vertices))
   allocate(area_third(n_vertices))
   allocate(on_boundary(n_vertices))
   K           = 0._R8P
   angle_sum   = 0._R8P
   area_third  = 0._R8P
   on_boundary = .false.
   had_degenerate = .false.

   ! ---- Per-triangle pass: accumulate angles, areas, and boundary flags.
   do f = 1_I4P, nf
      v1 = pool%facet_vid(facet_id=f, local_v=1_I4P)
      v2 = pool%facet_vid(facet_id=f, local_v=2_I4P)
      v3 = pool%facet_vid(facet_id=f, local_v=3_I4P)
      if (v1 == 0_I4P .or. v2 == 0_I4P .or. v3 == 0_I4P) cycle

      p1 = pool%coord(vid=v1)
      p2 = pool%coord(vid=v2)
      p3 = pool%coord(vid=v3)

      e12 = p2 - p1
      e23 = p3 - p2
      e31 = p1 - p3
      cross_full = e12%crossproduct(rhs=e23)
      area_2x = sqrt(cross_full%dotproduct(rhs=cross_full))
      if (area_2x < DEGENERATE_AREA_TOL) then
         had_degenerate = .true.
         cycle
      endif
      area_tri = 0.5_R8P * area_2x

      ! Per-vertex angles via atan2(|cross|, dot) — well-conditioned
      ! across the full [0, pi] range. The two edges meeting at each
      ! vertex are oriented outward from that vertex.
      ang1 = angle_at(u=e12,       v=-1._R8P * e31)  ! at v1: edges to v2 (e12) and to v3 (-e31)
      ang2 = angle_at(u=-1._R8P*e12, v=e23)            ! at v2: edges to v1 (-e12) and to v3 (e23)
      ang3 = angle_at(u=-1._R8P*e23, v=e31)            ! at v3: edges to v2 (-e23) and to v1 (e31)

      angle_sum(v1)  = angle_sum(v1)  + ang1
      angle_sum(v2)  = angle_sum(v2)  + ang2
      angle_sum(v3)  = angle_sum(v3)  + ang3
      area_third(v1) = area_third(v1) + area_tri / 3._R8P
      area_third(v2) = area_third(v2) + area_tri / 3._R8P
      area_third(v3) = area_third(v3) + area_tri / 3._R8P

      ! Boundary detection: any incident edge with fcon_edge == 0 flags
      ! the two vertices it touches.
      do e1 = 1_I4P, 3_I4P
         if (facet(f)%fcon_edge(e1) /= 0_I4P) cycle
         ! Edge e1 is disconnected. It touches vertices at local positions
         ! (e1) and (mod(e1, 3) + 1).
         select case (e1)
         case (1_I4P); on_boundary(v1) = .true.; on_boundary(v2) = .true.
         case (2_I4P); on_boundary(v2) = .true.; on_boundary(v3) = .true.
         case (3_I4P); on_boundary(v3) = .true.; on_boundary(v1) = .true.
         endselect
      enddo
   enddo

   ! ---- Per-vertex pass: compute K via the angle-defect formula.
   do i = 1_I4P, n_vertices
      if (area_third(i) <= 0._R8P) then
         K(i) = 0._R8P
         cycle
      endif
      if (on_boundary(i)) then
         K(i) = (PI - angle_sum(i)) / area_third(i)
      else
         K(i) = (2._R8P * PI - angle_sum(i)) / area_third(i)
      endif
   enddo

   deallocate(angle_sum, area_third, on_boundary)
   if (had_degenerate .and. present(status)) status = CURV_STATUS_DEGENERATE_TRIANGLE
   endsubroutine build_gaussian_curvature

   subroutine build_mean_curvature(facet, pool, H, status)
   !< Compute per-vertex **signed** mean curvature `H(1:n_vertices)` via
   !<
   !<   H_i * n_i  =  (1/2) M^-1 (L V)_i
   !<
   !< where `L` is the cotangent Laplacian, `M` is the diagonal
   !< barycentric mass, and `V` is the per-vertex coordinate field.
   !<
   !< Sign convention: H > 0 at convex-outward vertices (sphere
   !< surface, cube corners), H < 0 at saddle / concave-outward
   !< vertices. Sign is set by projecting the mean-curvature normal
   !< vector onto the per-vertex angle-weighted pseudo-normal
   !< (`facet%vertex_pnormal`, populated by `surface%analyze`); the
   !< magnitude is `(1/2) ||M^-1 L V||`.
   !<
   !< **Important note retracting an earlier claim**: the §2.4 commit
   !< message that shipped Gaussian curvature stated mean curvature
   !< needed a sparse linear solver. That was wrong. The mass matrix
   !< `M` is diagonal (barycentric per the §2.1 builder), so `M^-1`
   !< is just per-row scalar division. Mean H is shippable today
   !< without any external solver.
   !<
   !< Pre-conditions:
   !<   - `pool` initialized.
   !<   - `analyze` has run (pseudo-normals populated). All public
   !<     surface construction paths ensure both.
   !<
   !< Vertices with no surviving incident triangle (e.g. all triangles
   !< degenerate) get `H = 0` — defensible default.
   type(facet_object),       intent(in)            :: facet(:)
   type(vertex_pool_object), intent(in)            :: pool
   real(R8P), allocatable,   intent(out)           :: H(:)
   integer(I4P),             intent(out), optional :: status
   integer(I4P)                                    :: n_vertices, nf, f, i
   integer(I4P)                                    :: v1, v2, v3
   type(vector_R8P)                                :: p1, p2, p3
   type(vector_R8P)                                :: e12, e23, e31, cross_full
   real(R8P)                                       :: area_2x, area_tri
   real(R8P)                                       :: cot1, cot2, cot3
   type(vector_R8P), allocatable                   :: LV(:)
   real(R8P),        allocatable                   :: area_third(:)
   type(vector_R8P), allocatable                   :: pnormal_at(:)  !< Per-pool-vertex pseudo-normal proxy.
   logical, allocatable                            :: pnormal_set(:)
   type(vector_R8P)                                :: Hn, n_proxy
   real(R8P)                                       :: mag, sign_dot
   logical                                         :: had_degenerate

   if (present(status)) status = CURV_STATUS_OK
   nf = size(facet, kind=I4P)
   if (nf == 0_I4P .or. .not. pool%get_is_initialized()) then
      allocate(H(0))
      if (present(status)) status = CURV_STATUS_BAD_INPUT
      return
   endif

   n_vertices = pool%vertex_count()
   allocate(H(n_vertices))
   allocate(LV(n_vertices))
   allocate(area_third(n_vertices))
   allocate(pnormal_at(n_vertices))
   allocate(pnormal_set(n_vertices))
   H = 0._R8P
   do i = 1_I4P, n_vertices
      LV(i) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
      pnormal_at(i) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   enddo
   area_third  = 0._R8P
   pnormal_set = .false.
   had_degenerate = .false.

   ! Single per-triangle pass: assemble L*V, accumulate areas, and
   ! capture one pseudo-normal per pool vertex (any incident facet's
   ! vertex_pnormal — all are angle-weighted toward the same outward
   ! direction, so picking the first is fine for sign assignment).
   do f = 1_I4P, nf
      v1 = pool%facet_vid(facet_id=f, local_v=1_I4P)
      v2 = pool%facet_vid(facet_id=f, local_v=2_I4P)
      v3 = pool%facet_vid(facet_id=f, local_v=3_I4P)
      if (v1 == 0_I4P .or. v2 == 0_I4P .or. v3 == 0_I4P) cycle

      p1 = pool%coord(vid=v1)
      p2 = pool%coord(vid=v2)
      p3 = pool%coord(vid=v3)

      e12 = p2 - p1
      e23 = p3 - p2
      e31 = p1 - p3
      cross_full = e12%crossproduct(rhs=e23)
      area_2x = sqrt(cross_full%dotproduct(rhs=cross_full))
      if (area_2x < DEGENERATE_AREA_TOL) then
         had_degenerate = .true.
         cycle
      endif
      area_tri = 0.5_R8P * area_2x

      cot1 = (-(e31%dotproduct(rhs=e12))) / area_2x
      cot2 = (-(e12%dotproduct(rhs=e23))) / area_2x
      cot3 = (-(e23%dotproduct(rhs=e31))) / area_2x

      LV(v1) = LV(v1) + (0.5_R8P * cot3) * (p2 - p1) + (0.5_R8P * cot2) * (p3 - p1)
      LV(v2) = LV(v2) + (0.5_R8P * cot3) * (p1 - p2) + (0.5_R8P * cot1) * (p3 - p2)
      LV(v3) = LV(v3) + (0.5_R8P * cot2) * (p1 - p3) + (0.5_R8P * cot1) * (p2 - p3)

      area_third(v1) = area_third(v1) + area_tri / 3._R8P
      area_third(v2) = area_third(v2) + area_tri / 3._R8P
      area_third(v3) = area_third(v3) + area_tri / 3._R8P

      ! First-touch capture of pseudo-normals for sign assignment.
      if (.not. pnormal_set(v1)) then
         pnormal_at(v1) = facet(f)%vertex_pnormal(1_I4P)
         pnormal_set(v1) = .true.
      endif
      if (.not. pnormal_set(v2)) then
         pnormal_at(v2) = facet(f)%vertex_pnormal(2_I4P)
         pnormal_set(v2) = .true.
      endif
      if (.not. pnormal_set(v3)) then
         pnormal_at(v3) = facet(f)%vertex_pnormal(3_I4P)
         pnormal_set(v3) = .true.
      endif
   enddo

   ! Per-vertex post-process: H_i * n_i = (1/2) M^-1 LV_i, then signed
   ! magnitude via projection onto pseudo-normal.
   do i = 1_I4P, n_vertices
      if (area_third(i) <= 0._R8P) cycle
      Hn = (0.5_R8P / area_third(i)) * LV(i)
      mag = sqrt(Hn%dotproduct(rhs=Hn))
      if (mag <= 0._R8P) then
         H(i) = 0._R8P
         cycle
      endif
      n_proxy = pnormal_at(i)
      sign_dot = Hn%dotproduct(rhs=n_proxy)
      ! Sign convention: H > 0 outward-convex. With FOSSIL's
      ! positive-semidefinite cotangent Laplacian (L_ii < 0 on the
      ! diagonal), H*n = (1/2) M^-1 L V points INWARD on convex
      ! surfaces — opposite to the outward pseudo-normal. So
      ! `sign_dot < 0` means convex outward → assign H positive.
      if (sign_dot <= 0._R8P) then
         H(i) =  mag
      else
         H(i) = -mag
      endif
   enddo

   deallocate(LV, area_third, pnormal_at, pnormal_set)
   if (had_degenerate .and. present(status)) status = CURV_STATUS_DEGENERATE_TRIANGLE
   endsubroutine build_mean_curvature

   pure function angle_at(u, v) result(theta)
   !< Angle between vectors u and v via atan2(|u x v|, u . v). Well-
   !< conditioned across the full [0, pi] range — beats acos(dot/norms)
   !< which loses precision near 0 and pi.
   type(vector_R8P), intent(in) :: u, v
   real(R8P)                    :: theta
   type(vector_R8P)             :: c
   real(R8P)                    :: cross_mag, dot

   c = u%crossproduct(rhs=v)
   cross_mag = sqrt(c%dotproduct(rhs=c))
   dot = u%dotproduct(rhs=v)
   theta = atan2(cross_mag, dot)
   endfunction angle_at

endmodule fossil_curvature
