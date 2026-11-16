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
use fossil_shape_operator, only : eig_sym2_oac
use penf, only : I4P, R8P
use vecfor, only : vector_R8P
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

implicit none
private

public :: build_gaussian_curvature
public :: build_mean_curvature
public :: build_jet_curvature
public :: menger_curvature, vertex_curvature_2d
public :: CURV_STATUS_OK, CURV_STATUS_BAD_INPUT, CURV_STATUS_DEGENERATE_TRIANGLE, CURV_STATUS_FIT_FAILED

integer(I4P), parameter :: CURV_STATUS_OK                   = 0_I4P
integer(I4P), parameter :: CURV_STATUS_BAD_INPUT            = 1_I4P  !< Empty facet array or pool not initialized.
integer(I4P), parameter :: CURV_STATUS_DEGENERATE_TRIANGLE  = 2_I4P  !< At least one triangle has area below tolerance; affected vertices skip that triangle.
integer(I4P), parameter :: CURV_STATUS_FIT_FAILED           = 3_I4P  !< Jet fit rank deficient at a vertex: H, K NaN there (#91).

real(R8P), parameter :: DEGENERATE_AREA_TOL = 1.0e-30_R8P
real(R8P), parameter :: JET_MERGE_TOL = 1.0e-8_R8P !< Ring vertices closer than this times the ring radius are one point.
integer(I4P), parameter :: JET_MAX_RINGS = 4_I4P   !< Largest ring grown when merged points leave too few for the fit.
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

   subroutine build_jet_curvature(facet, pool, H, K, status, principal, directions)
   !< Per-vertex mean and Gaussian curvature by a local polynomial fit (osculating jet, Cazals &
   !< Pouget 2005; issue #42).
   !<
   !< At each vertex p, in the frame (t1, t2, n) of its angle-weighted pseudo-normal n, the 2-ring
   !< vertices are fitted in the least-squares sense by the height function
   !<   z = d x + e y + a x^2 + b x y + c y^2 [+ cubic terms]
   !< (cubic terms when the 2-ring has >= 12 vertices, else quadratic; Householder QR on coordinates
   !< scaled by the RMS neighbour distance). The curvatures are those of the graph at the origin (the
   !< linear terms absorb a tilted n): with W = sqrt(1 + d^2 + e^2), the second fundamental form
   !< (2a, b, 2c) / W over the first (1 + d^2, d e, 1 + e^2).
   !<
   !< Accuracy: the cubic jet gives curvatures with error O(h^2) on any sampling, the quadratic O(h);
   !< unlike the cotangent / angle-defect operators of `build_mean_curvature` /
   !< `build_gaussian_curvature`, which do not converge pointwise on irregular meshes (on icospheres
   !< their maximum error stays near 15 % at the valence-5 vertices). A vertex with fewer than 6
   !< 2-ring vertices takes the cotangent / angle-defect value. Flat neighbourhoods give 0 exactly.
   !<
   !< Conditioning (issue #91): ring vertices within 1e-8 of the ring radius of the centre or of an earlier
   !< ring vertex are one point (marching cubes leaves vertices 1e-12 apart around a grid node 1e-12 off the
   !< level); if merging removed points and fewer than 12 distinct ones remain, the ring grows by one
   !< neighbourhood at a time, up to the 4-ring. The columns of the least-squares matrix are scaled by powers
   !< of two (exact) before the QR, whose rank check (a diagonal of R below 1e-12 of the largest) is then
   !< relative to equilibrated columns. If the cubic and the quadratic fits are both rank deficient, or
   !< merging left fewer than 6 distinct points, H and K are NaN at that vertex (not a number, rather than
   !< a wrong one) and `status` = CURV_STATUS_FIT_FAILED (a ring on two lines, for instance; a ring with all
   !< heights exactly 0 still gives 0). Without merged points the ring and the fit are those of before.
   !<
   !< Sign convention as `build_mean_curvature`: H > 0 where the surface is convex seen from the side
   !< the normals point to (sphere: H = 1/R, K = 1/R^2).
   !<
   !< Shape operator (roadmap #109 WP1 step 03, optional `principal(2, nv)`, `directions(3, 2, nv)`; H and K unchanged bit
   !< for bit): the fitted graph r(x, y) = x t1 + y t2 + z(x, y) n has, at the origin, the tangents r_x = (1, 0, d),
   !< r_y = (0, 1, e) (frame coordinates), the first fundamental form I = J^T J (J = [r_x r_y]) and the second
   !< II = [[2a, b], [b, 2c]] / W about the graph normal (-d, -e, 1) / W. In the orthonormal tangent basis u1 = r_x / |r_x|,
   !< u2 = N x u1, whose parameter coordinates are c1 = (1, 0) / l1 and c2 = (-d e, 1 + d^2) / (W l1), l1 = sqrt(1 + d^2)
   !< (the first two frame components of a tangent vector are its parameter coordinates, since J's top block is the
   !< identity), the shape operator is the symmetric S_ij = -c_i^T II c_j / hs (sign: convex positive; hs undoes the
   !< coordinate scaling). Its eigenvalues are those of I^-1 II, so k1 + k2 = 2 H and k1 k2 = K to round-off; its
   !< eigenvectors (`fossil_shape_operator%eig_sym2_oac`, closed form, an umbilic takes u1) are the principal directions
   !< of the fitted surface. They lie in the graph's tangent plane, tilted from the pseudo-normal's by the linear terms
   !< (d, e): each is stored projected on the plane normal to the pseudo-normal n and normalised, d1 first, d2 = n x d1
   !< (orthonormal, (d1, d2, n) right-handed). Vertices without a fit (fewer than 6 ring vertices: the cotangent /
   !< angle-defect H, K; a flat rank-deficient ring: 0; a failed fit: NaN) take k = H +- sqrt(max(H^2 - K, 0)) and the
   !< jet frame (t1, t2) as directions; vertices of no facet take zero directions.
   type(facet_object),       intent(in)            :: facet(:) !< Facets.
   type(vertex_pool_object), intent(in)            :: pool     !< Vertex pool.
   real(R8P), allocatable,   intent(out)           :: H(:)     !< Mean curvature per pool vertex.
   real(R8P), allocatable,   intent(out)           :: K(:)     !< Gaussian curvature per pool vertex.
   integer(I4P),             intent(out), optional :: status   !< CURV_STATUS_*.
   real(R8P), allocatable,   intent(out), optional :: principal(:,:)    !< Principal curvatures k1 >= k2 [2, nv].
   real(R8P), allocatable,   intent(out), optional :: directions(:,:,:) !< Principal directions d1, d2 [3, 2, nv].
   integer(I4P), parameter       :: MAXN = 256_I4P              !< Neighbour cap per vertex (2-ring).
   real(R8P),    allocatable     :: hc(:), kc(:)                !< Fallback curvatures.
   real(R8P),    allocatable     :: xyz(:,:), pn(:,:)           !< Coordinates, pseudo-normals.
   logical,      allocatable     :: pn_set(:)                   !< Pseudo-normal captured.
   integer(I4P), allocatable     :: deg(:), off(:), adj(:), stamp(:) !< One-ring CSR, dedup stamps.
   integer(I4P), allocatable     :: fill(:)                     !< CSR fill pointers.
   integer(I4P)                  :: ring(MAXN)                  !< 2-ring of the current vertex.
   real(R8P)                     :: a(MAXN, 9), rhs(MAXN), coef(9) !< Least-squares system.
   real(R8P)                     :: t1(3), t2(3), nn(3), q(3), x, y, z, hs, w2, eg, hh, kk !< Frame, work.
   integer(I4P)                  :: nv, nf, f, i, j, m, v(3), l, ii, jj, np, nr !< Counters.
   logical                       :: ok                          !< Fit succeeded.
   integer(I4P)                  :: st                          !< Builder status.
   type(vector_R8P)              :: pv                          !< Pool vertex.
   integer(I4P)                  :: keep(MAXN), nd, nd_at, lo, hi, layer !< Distinct ring vertices (for nd_at), layers.
   logical                       :: failed                      !< Some vertex fit failed.
   logical                       :: want_sd                     !< Shape operator requested.
   logical,      allocatable     :: fitted(:)                   !< Vertex shape operator from its fit.

   if (present(status)) status = CURV_STATUS_OK
   nf = size(facet, kind=I4P)
   want_sd = present(principal) .and. present(directions)
   if (present(principal)) allocate(principal(2, 0))
   if (present(directions)) allocate(directions(3, 2, 0))
   if (nf == 0_I4P .or. .not. pool%get_is_initialized()) then
      allocate(H(0), K(0))
      if (present(status)) status = CURV_STATUS_BAD_INPUT
      return
   endif
   call build_mean_curvature(facet=facet, pool=pool, H=hc, status=st)
   call build_gaussian_curvature(facet=facet, pool=pool, K=kc)
   if (present(status)) status = st
   nv = pool%vertex_count()
   allocate(H(nv), K(nv), xyz(3, nv), pn(3, nv), pn_set(nv), deg(nv), off(nv + 1), stamp(nv), fill(nv))
   H = hc ; K = kc
   if (want_sd) then
      deallocate(principal, directions) ; allocate(principal(2, nv), directions(3, 2, nv), fitted(nv))
      principal = 0._R8P ; directions = 0._R8P ; fitted = .false.
   endif
   do i = 1_I4P, nv
      pv = pool%coord(i) ; xyz(:, i) = [pv%x, pv%y, pv%z]
   enddo
   pn = 0._R8P ; pn_set = .false. ; deg = 0_I4P
   do f = 1_I4P, nf ! one-ring with duplicates (each edge from both facets), pseudo-normals
      do l = 1_I4P, 3_I4P
         v(l) = pool%facet_vid(facet_id=f, local_v=l)
      enddo
      if (any(v <= 0_I4P)) cycle
      do l = 1_I4P, 3_I4P
         deg(v(l)) = deg(v(l)) + 2_I4P
         if (.not. pn_set(v(l))) then
            pn(:, v(l)) = [facet(f)%vertex_pnormal(l)%x, facet(f)%vertex_pnormal(l)%y, facet(f)%vertex_pnormal(l)%z]
            pn_set(v(l)) = .true.
         endif
      enddo
   enddo
   off(1) = 1_I4P
   do i = 1_I4P, nv
      off(i + 1) = off(i) + deg(i)
   enddo
   allocate(adj(off(nv + 1) - 1_I4P))
   fill = off(1:nv)
   do f = 1_I4P, nf
      do l = 1_I4P, 3_I4P
         v(l) = pool%facet_vid(facet_id=f, local_v=l)
      enddo
      if (any(v <= 0_I4P)) cycle
      do l = 1_I4P, 3_I4P
         adj(fill(v(l))) = v(mod(l, 3_I4P) + 1_I4P) ; adj(fill(v(l)) + 1_I4P) = v(mod(l + 1_I4P, 3_I4P) + 1_I4P)
         fill(v(l)) = fill(v(l)) + 2_I4P
      enddo
   enddo
   stamp = 0_I4P
   failed = .false.
   do i = 1_I4P, nv
      if (.not. pn_set(i)) cycle
      nn = pn(:, i) ; hs = norm2(nn)
      if (hs <= 0._R8P) cycle
      nn = nn / hs
      ! 2-ring, deduplicated with stamps (the vertex itself excluded)
      nr = 0_I4P ; stamp(i) = i
      do ii = off(i), off(i + 1) - 1_I4P
         j = adj(ii)
         if (stamp(j) /= i .and. nr < MAXN) then
            stamp(j) = i ; nr = nr + 1_I4P ; ring(nr) = j
         endif
      enddo
      lo = 1_I4P ; hi = nr ; nd_at = -1_I4P
      do layer = 2_I4P, JET_MAX_RINGS ! layer 2 always; more only while merged points leave < 12 distinct ones
         if (layer > 2_I4P) then
            call distinct_ring(xyz=xyz, centre=i, ring=ring(1:nr), keep=keep, nd=nd) ; nd_at = nr
            if (nd == nr .or. nd >= 12_I4P .or. hi < lo) exit
         endif
         m = nr
         do ii = lo, hi
            do jj = off(ring(ii)), off(ring(ii) + 1) - 1_I4P
               j = adj(jj)
               if (stamp(j) /= i .and. nr < MAXN) then
                  stamp(j) = i ; nr = nr + 1_I4P ; ring(nr) = j
               endif
            enddo
         enddo
         lo = m + 1_I4P ; hi = nr
      enddo
      if (nr < 6_I4P) cycle ! keep the fallback
      if (nd_at /= nr) call distinct_ring(xyz=xyz, centre=i, ring=ring(1:nr), keep=keep, nd=nd)
      if (nd < nr) then ! merged: the fit runs on the distinct points
         if (nd < 6_I4P) then
            H(i) = ieee_value(1._R8P, ieee_quiet_nan) ; K(i) = H(i) ; failed = .true.
            cycle
         endif
         ring(1:nd) = ring(keep(1:nd)) ; nr = nd
      endif
      ! tangent frame
      if (abs(nn(1)) < 0.9_R8P) then
         t1 = [0._R8P, -nn(3), nn(2)] ! e_x x n
      else
         t1 = [nn(3), 0._R8P, -nn(1)] ! e_y x n
      endif
      t1 = t1 / norm2(t1)
      t2 = [nn(2) * t1(3) - nn(3) * t1(2), nn(3) * t1(1) - nn(1) * t1(3), nn(1) * t1(2) - nn(2) * t1(1)]
      hs = 0._R8P
      do ii = 1_I4P, nr
         hs = hs + sum((xyz(:, ring(ii)) - xyz(:, i))**2)
      enddo
      hs = sqrt(hs / nr)
      if (hs <= 0._R8P) cycle
      np = 5_I4P ; if (nr >= 12_I4P) np = 9_I4P
      do ! cubic jet if the 2-ring allows it, else (or if rank deficient) quadratic
         do ii = 1_I4P, nr
            q = (xyz(:, ring(ii)) - xyz(:, i)) / hs
            x = dot_product(q, t1) ; y = dot_product(q, t2) ; z = dot_product(q, nn)
            a(ii, 1:5) = [x, y, x * x, x * y, y * y]
            if (np == 9_I4P) a(ii, 6:9) = [x**3, x * x * y, x * y * y, y**3]
            rhs(ii) = z
         enddo
         call householder_lsq(a, rhs, nr, np, coef, ok)
         if (ok .or. np == 5_I4P) exit
         np = 5_I4P
      enddo
      if (.not. ok) then
         if (any([(abs(dot_product(xyz(:, ring(ii)) - xyz(:, i), nn)), ii = 1_I4P, nr)] > 0._R8P)) then
            H(i) = ieee_value(1._R8P, ieee_quiet_nan) ; K(i) = H(i) ; failed = .true.
         else ! every height exactly 0: flat
            H(i) = 0._R8P ; K(i) = 0._R8P
         endif
         cycle
      endif
      ! graph z = f(x, y) at the origin: f_x = d, f_y = e, f_xx = 2a, f_xy = b, f_yy = 2c
      w2 = 1._R8P + coef(1)**2 + coef(2)**2
      eg = w2 ! (1 + d^2)(1 + e^2) - (d e)^2
      kk = (4._R8P * coef(3) * coef(5) - coef(4)**2) / (w2 * eg)
      hh = ((1._R8P + coef(2)**2) * 2._R8P * coef(3) - 2._R8P * coef(1) * coef(2) * coef(4) + &
            (1._R8P + coef(1)**2) * 2._R8P * coef(5)) / (2._R8P * sqrt(w2) * eg)
      ! z points along the outward normal: a convex surface bends away from it (hh < 0), so H = -hh
      H(i) = -hh / hs
      K(i) = kk / (hs * hs)
      if (want_sd) then
         call jet_shape_operator(coef(1:5), hs, t1, t2, nn, principal(:, i), directions(:, :, i))
         fitted(i) = .true.
      endif
   enddo
   if (want_sd) then ! the vertices without a fit: k from H and K, the jet frame as directions
      do i = 1_I4P, nv
         if (fitted(i)) cycle
         hh = sqrt(max(H(i) * H(i) - K(i), 0._R8P))
         principal(:, i) = [H(i) + hh, H(i) - hh]
         if (.not. pn_set(i)) cycle
         nn = pn(:, i) ; hs = norm2(nn)
         if (hs <= 0._R8P) cycle
         nn = nn / hs
         if (abs(nn(1)) < 0.9_R8P) then
            t1 = [0._R8P, -nn(3), nn(2)]
         else
            t1 = [nn(3), 0._R8P, -nn(1)]
         endif
         t1 = t1 / norm2(t1)
         t2 = [nn(2) * t1(3) - nn(3) * t1(2), nn(3) * t1(1) - nn(1) * t1(3), nn(1) * t1(2) - nn(2) * t1(1)]
         directions(:, 1, i) = t1 ; directions(:, 2, i) = t2
      enddo
   endif
   if (failed .and. present(status)) status = CURV_STATUS_FIT_FAILED
   endsubroutine build_jet_curvature

   pure subroutine jet_shape_operator(coef, hs, t1, t2, nn, k, d)
   !< Principal curvatures and directions of a fitted jet z = d x + e y + a x^2 + b x y + c y^2 + ... at the origin (frame
   !< (t1, t2, nn), coordinates scaled by hs), see `build_jet_curvature`: S_ij = -c_i^T II c_j / hs in the orthonormal
   !< tangent basis (u1, u2) of the graph, `eig_sym2_oac`, the eigenvectors projected on the plane normal to nn.
   real(R8P), intent(in)  :: coef(5)      !< (d, e, a, b, c).
   real(R8P), intent(in)  :: hs           !< Coordinate scale.
   real(R8P), intent(in)  :: t1(3), t2(3), nn(3) !< Frame.
   real(R8P), intent(out) :: k(2)         !< k1 >= k2.
   real(R8P), intent(out) :: d(3, 2)      !< d1, d2 (unit, normal to nn).
   real(R8P)              :: dd, ee, w, l1, p11, p12, p22, s11, s12, s22, c, s, u1(3), u2(3), g(3), gl !< Work.

   dd = coef(1) ; ee = coef(2)
   w = sqrt(1._R8P + dd * dd + ee * ee) ; l1 = sqrt(1._R8P + dd * dd)
   p11 = -2._R8P * coef(3) / (w * hs) ; p12 = -coef(4) / (w * hs) ; p22 = -2._R8P * coef(5) / (w * hs) ! -II / hs
   s11 = p11 / (l1 * l1)
   s12 = (p11 * (-dd * ee) + p12 * (1._R8P + dd * dd)) / (w * l1 * l1)
   s22 = (p11 * (dd * ee)**2 - 2._R8P * p12 * dd * ee * (1._R8P + dd * dd) + p22 * (1._R8P + dd * dd)**2) / (w * w * l1 * l1)
   call eig_sym2_oac(s11, s12, s22, k(1), k(2), c, s)
   u1 = [1._R8P, 0._R8P, dd] / l1                                ! frame coordinates
   u2 = [-dd * ee, 1._R8P + dd * dd, ee] / (w * l1)
   g = c * u1 + s * u2                                            ! d1 in the graph's tangent plane
   g = g(1) * t1 + g(2) * t2                                      ! projected on the plane normal to nn, world coordinates
   gl = norm2(g)
   if (gl > 0._R8P) then
      d(:, 1) = g / gl
   else
      d(:, 1) = t1
   endif
   d(:, 2) = [nn(2) * d(3, 1) - nn(3) * d(2, 1), nn(3) * d(1, 1) - nn(1) * d(3, 1), nn(1) * d(2, 1) - nn(2) * d(1, 1)]
   endsubroutine jet_shape_operator

   pure function menger_curvature(ua, va, ub, vb, uc, vc) result(kappa)
   !< Signed curvature of the circle through three points a, b, c of a planar polyline (roadmap #109 WP1 step 03):
   !< 2 ((b - a) x (c - b)) / (|b - a| |c - b| |c - a|), positive for a left turn (a counter-clockwise outer loop, the
   !< solid on the left: convex seen from outside). 0 if two points coincide. At b with uniform spacing the error is
   !< O(h^2); with spacings h1, h2 it has an O(h2 - h1) term.
   real(R8P), intent(in) :: ua, va, ub, vb, uc, vc !< Points.
   real(R8P)             :: kappa                  !< Curvature.
   real(R8P)             :: lab, lbc, lac          !< Lengths.

   kappa = 0._R8P
   lab = sqrt((ub - ua)**2 + (vb - va)**2) ; lbc = sqrt((uc - ub)**2 + (vc - vb)**2) ; lac = sqrt((uc - ua)**2 + (vc - va)**2)
   if (lab <= 0._R8P .or. lbc <= 0._R8P .or. lac <= 0._R8P) return
   kappa = 2._R8P * ((ub - ua) * (vc - vb) - (vb - va) * (uc - ub)) / (lab * lbc * lac)
   endfunction menger_curvature

   pure function vertex_curvature_2d(u, v) result(kappa)
   !< Signed curvature at the middle of five consecutive polyline vertices p1..p5 (roadmap #109 WP1 step 03): the circle
   !< through p2, p3, p4 (`menger_curvature`) minus its first-order spacing term, kappa_M - kappa_s (h2 - h1) / 3 (h1 =
   !< |p3 - p2|, h2 = |p4 - p3|; the circumcircle of three points at arc lengths -h1, 0, h2 has curvature
   !< kappa + kappa_s (h2 - h1) / 3 + O(h^2)), with the arc-length derivative kappa_s of a cubic jet: in the frame of the
   !< chord t = (p4 - p2) / |p4 - p2| and its left normal (the solid side) the four neighbours are fitted in the
   !< least-squares sense by y = d x + a x^2 + c x^3 (Householder QR on coordinates scaled by their RMS distance), then
   !< kappa_s = (y3 (1 + y1^2) - 3 y1 y2^2) / (1 + y1^2)^3 at x = 0 (y1, y2, y3 the derivatives). Error O(h^2) on any
   !< spacing (the plain circle has an O(h2 - h1) term; the jet's own curvature is O(h^2) too, with a constant about 2-3
   !< times larger on ellipses); equal to the circle on uniform spacing (exact on a uniformly sampled circle). Positive for
   !< a left turn. Rank deficient or degenerate fit: the plain circle.
   real(R8P), intent(in) :: u(5), v(5) !< Vertices p1..p5 (p3 the vertex).
   real(R8P)             :: kappa      !< Curvature.
   real(R8P)             :: tu, tv, l, hs, a(4, 3), b(4), coef(3), y1, y2, y3, h1, h2 !< Frame, scale, system, jet.
   integer(I4P)          :: i, r       !< Counters.
   logical               :: ok         !< Fit succeeded.

   kappa = menger_curvature(u(2), v(2), u(3), v(3), u(4), v(4))
   h1 = sqrt((u(3) - u(2))**2 + (v(3) - v(2))**2) ; h2 = sqrt((u(4) - u(3))**2 + (v(4) - v(3))**2)
   if (h1 == h2) return ! no first-order term
   tu = u(4) - u(2) ; tv = v(4) - v(2) ; l = sqrt(tu * tu + tv * tv)
   if (l <= 0._R8P) return
   tu = tu / l ; tv = tv / l
   hs = 0._R8P
   do i = 1_I4P, 5_I4P
      hs = hs + (u(i) - u(3))**2 + (v(i) - v(3))**2
   enddo
   hs = sqrt(hs / 4._R8P)
   if (hs <= 0._R8P) return
   r = 0_I4P
   do i = 1_I4P, 5_I4P
      if (i == 3_I4P) cycle
      r = r + 1_I4P
      a(r, 1) = ((u(i) - u(3)) * tu + (v(i) - v(3)) * tv) / hs
      b(r) = (-(u(i) - u(3)) * tv + (v(i) - v(3)) * tu) / hs
      a(r, 2) = a(r, 1) * a(r, 1) ; a(r, 3) = a(r, 2) * a(r, 1)
   enddo
   call householder_lsq(a, b, 4_I4P, 3_I4P, coef, ok)
   if (.not. ok) return
   y1 = coef(1) ; y2 = 2._R8P * coef(2) / hs ; y3 = 6._R8P * coef(3) / (hs * hs)
   kappa = kappa - (y3 * (1._R8P + y1 * y1) - 3._R8P * y1 * y2 * y2) / (1._R8P + y1 * y1)**3 * (h2 - h1) / 3._R8P
   endfunction vertex_curvature_2d

   pure subroutine distinct_ring(xyz, centre, ring, keep, nd)
   !< Positions in `ring` of its distinct vertices (issue #91): a vertex within JET_MERGE_TOL of the ring
   !< radius (largest distance from the centre) of the centre or of an earlier kept vertex is dropped.
   !< Ring order is kept, so a ring without such vertices gives keep = 1..size(ring).
   real(R8P),    intent(in)  :: xyz(:,:) !< Coordinates [3, n].
   integer(I4P), intent(in)  :: centre   !< Centre vertex.
   integer(I4P), intent(in)  :: ring(:)  !< Ring vertices.
   integer(I4P), intent(out) :: keep(:)  !< Positions in ring of the distinct vertices [nd].
   integer(I4P), intent(out) :: nd       !< Distinct vertices.
   real(R8P)                 :: tol2     !< Squared merge distance.
   integer(I4P)              :: a, b     !< Counters.
   logical                   :: dup      !< Vertex a merged.

   tol2 = 0._R8P
   do a = 1_I4P, size(ring, kind=I4P)
      tol2 = max(tol2, sum((xyz(:, ring(a)) - xyz(:, centre))**2))
   enddo
   tol2 = JET_MERGE_TOL**2 * tol2
   nd = 0_I4P
   do a = 1_I4P, size(ring, kind=I4P)
      dup = sum((xyz(:, ring(a)) - xyz(:, centre))**2) <= tol2
      do b = 1_I4P, nd
         if (dup) exit
         dup = sum((xyz(:, ring(a)) - xyz(:, ring(keep(b))))**2) <= tol2
      enddo
      if (dup) cycle
      nd = nd + 1_I4P ; keep(nd) = a
   enddo
   endsubroutine distinct_ring

   pure subroutine householder_lsq(a, b, m, n, x, ok)
   !< Least-squares solution of a(1:m, 1:n) x = b(1:m) by Householder QR (m >= n, n <= 9). `ok` is false
   !< if m < n or the system is rank deficient (a diagonal of R below 1e-12 of the largest). The columns
   !< are first scaled by powers of two to unit size (exact, so the solution's bits do not change; issue
   !< #91), which makes the rank check relative to equilibrated columns.
   real(R8P),    intent(inout) :: a(:,:) !< Matrix (overwritten).
   real(R8P),    intent(inout) :: b(:)   !< Right-hand side (overwritten).
   integer(I4P), intent(in)    :: m, n   !< Rows, columns.
   real(R8P),    intent(out)   :: x(:)   !< Solution [n].
   logical,      intent(out)   :: ok     !< Success.
   real(R8P)                   :: alpha, vnorm2, s, rmax !< Work.
   real(R8P)                   :: v(size(a, 1))           !< Householder vector.
   integer(I4P)                :: j, k                    !< Counters.
   integer(I4P)                :: ce(n)                   !< Column scale exponents.

   x = 0._R8P ; ok = .false.
   if (m < n) return
   do j = 1_I4P, n
      ce(j) = exponent(maxval(abs(a(1:m, j))))
      a(1:m, j) = scale(a(1:m, j), -ce(j))
   enddo
   rmax = 0._R8P
   do k = 1_I4P, n
      alpha = -sign(norm2(a(k:m, k)), a(k, k))
      v(k:m) = a(k:m, k) ; v(k) = v(k) - alpha
      vnorm2 = dot_product(v(k:m), v(k:m))
      if (vnorm2 > 0._R8P) then
         do j = k, n
            s = 2._R8P * dot_product(v(k:m), a(k:m, j)) / vnorm2
            a(k:m, j) = a(k:m, j) - s * v(k:m)
         enddo
         s = 2._R8P * dot_product(v(k:m), b(k:m)) / vnorm2
         b(k:m) = b(k:m) - s * v(k:m)
      endif
      rmax = max(rmax, abs(a(k, k)))
   enddo
   do k = 1_I4P, n
      if (abs(a(k, k)) <= 1.e-12_R8P * rmax) return
   enddo
   do k = n, 1_I4P, -1_I4P
      x(k) = (b(k) - dot_product(a(k, k + 1:n), x(k + 1:n))) / a(k, k)
   enddo
   do k = 1_I4P, n
      x(k) = scale(x(k), -ce(k))
   enddo
   ok = .true.
   endsubroutine householder_lsq

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
