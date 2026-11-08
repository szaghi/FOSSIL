!< Ray-mesh intersection queries (issue #18 §2.5): the hit record, the exact ray-facet kernel (issue #82) and the flat
!< O(n) reference driver.

module fossil_ray_query
!< Ray-mesh intersection queries — driver that consumes a `facet_object` array
!< and returns hit lists for first / all / any queries.
!<
!< Exact hits (issue #82). Every decision is exact, so a ray crosses a closed surface once per crossing: a ray through
!< a shared edge or vertex hits one of the facets there, a ray touching the surface at an edge or vertex without
!< crossing it, or lying in a face plane, hits none or two, and the parity of the hit count is the state of the
!< origin. Before #82 Moller-Trumbore with closed barycentric tests hit every facet at a shared edge (12 hits for 2
!< crossings through the face centres of cube-01) and could miss both at round-off.
!<
!< The ray o + t d is the oriented line through o and the far point q = o + s d (rounded once, `ray_prepare`), with s
!< such that q is beyond the box of the surface: every decision is then about two points in doubles.
!<
!< - Axis-aligned rays (two zero components of d) are decided by `tri_line_hit_oac`, the kernel of `line_crossings`
!<   and of the cut cells: the line translated by the global perturbation (eps1, eps2, eps3) of issue #69.
!< - Other rays: the line crosses triangle ABC iff the three Plucker side tests orient3d(o, q, A, B), (o, q, B, C),
!<   (o, q, C, A) have one sign (Woop et al. 2013, made exact). A zero side test (the line meets the edge line) takes
!<   one symbolic perturbation of the ray, o + delta and q + delta' with delta = (e1, e2, e4), delta' = (e8, e16, e32)
!<   (e^k: powers of one infinitesimal, so every product has its own order): orient3d(o + delta, q + delta', X, Y) =
!<   orient3d(o, q, X, Y) + g . delta + g' . delta' + (X - Y) . (delta x delta'), g = (q - Y) x (X - Y), g' = (X - Y) x
!<   (o - Y), whose sign is the sign of the first nonzero term in the order g_x, g_y, g_z, g'_x, -(X - Y)_z, (X - Y)_y,
!<   g'_y, -(X - Y)_x (orders e^1, 2, 4, 8, 10, 12, 16, 20). Each term is an exact orient2d or a comparison, the
!<   expansion is antisymmetric in (X, Y) as orient3d is, and X /= Y ends it: the perturbed line meets no edge, and
!<   the two facets of an edge see it on opposite sides.
!<
!< Both rays: the crossing is ahead iff o and q lie on opposite sides of the facet's plane (orient3d, exact; q is
!< beyond every facet). An origin exactly on the facet's plane is a hit at t = 0 (the crossing is at o): hits at t = 0
!< are kept, so a ray cast from the surface hits the facet it starts on, unless that facet is passed as `exclude`.
!< A facet with collinear vertices is never crossed. The topology is exact, t is rounded: the crossing point is the
!< barycentric combination by the side tests (by the line kernel's weights on axis rays), t its projection on d.
!<
!< The tree traversals prune boxes with the slab test on boxes padded by `ray%pad` (1e-9 of the ray's extent over
!< the surface box, plus the octree's facet reach): a crossing on a box face is never lost to round-off.

use fossil_facet_object, only : facet_object
use fossil_line_kernel,  only : tri_line_hit_oac, first_normal_sign
use fossil_predicates,   only : orient2d_xy, orient3d_xyz
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
private

public :: ray_hit_t
public :: ray_exact_t, ray_prepare, ray_facet_crossing, sort_hits_by_t
public :: ray_intersect_all_flat
public :: RAY_STATUS_OK, RAY_STATUS_BAD_INPUT

integer(I4P), parameter :: RAY_STATUS_OK        = 0_I4P
integer(I4P), parameter :: RAY_STATUS_BAD_INPUT = 1_I4P  !< Empty facet array, or |dir| < tolerance.

real(R8P), parameter :: DIR_TOL = 1.0e-30_R8P  !< |dir|^2 below this is rejected as degenerate.

type :: ray_hit_t
   !< One ray-facet intersection record.
   integer(I4P)     :: facet_id = 0_I4P     !< 1-based index into the facet array passed in.
   real(R8P)        :: t        = 0._R8P    !< Parametric distance along the ray.
   type(vector_R8P) :: point                !< Hit point in world coordinates.
endtype ray_hit_t

type :: ray_exact_t
   !< A ray prepared for the exact kernel (issue #82): origin, far point beyond the surface box, direction.
   real(R8P)    :: o(3) = 0._R8P   !< Origin.
   real(R8P)    :: q(3) = 0._R8P   !< Far point o + s d, rounded, beyond the surface box.
   real(R8P)    :: d(3) = 0._R8P   !< Direction (as given: t is parametric).
   integer(I4P) :: axis = 0_I4P    !< 1, 2, 3 for an axis-aligned ray (two zero components of d), else 0.
   real(R8P)    :: pad  = 0._R8P   !< Pad of the boxes in the slab test.
endtype ray_exact_t

contains

   pure function ray_prepare(ray_origin, ray_direction, bmin, bmax, reach) result(ray)
   !< Prepare a ray for the exact kernel (issue #82): the far point q = o + s d beyond the box [bmin, bmax] (the surface
   !< box), s = 2 R / |d| with R >= the distance from o to any point of the box (s doubled until q /= o), and the pad
   !< of the slab test: 1e-9 R plus `reach` (the octree's facet reach; 0 for the BVH).
   type(vector_R8P), intent(in)           :: ray_origin    !< Ray origin.
   type(vector_R8P), intent(in)           :: ray_direction !< Ray direction (nonzero).
   type(vector_R8P), intent(in)           :: bmin, bmax    !< Box of the surface.
   real(R8P),        intent(in), optional :: reach         !< Facet reach beyond its tree box (octree).
   type(ray_exact_t)                      :: ray           !< Prepared ray.
   real(R8P)                              :: lo(3), hi(3), r, s, dn !< Work.
   integer(I4P)                           :: k             !< Counter.

   ray%o = [ray_origin%x, ray_origin%y, ray_origin%z]
   ray%d = [ray_direction%x, ray_direction%y, ray_direction%z]
   lo = [bmin%x, bmin%y, bmin%z] ; hi = [bmax%x, bmax%y, bmax%z]
   r = sqrt(sum(max(abs(lo - ray%o), abs(hi - ray%o))**2))
   r = r + 1.e-6_R8P * maxval(abs(ray%o)) + tiny(1._R8P)
   dn = sqrt(sum(ray%d**2))
   s = 2._R8P * r / dn
   do k = 1, 64
      ray%q = ray%o + s * ray%d
      if (any(ray%q /= ray%o)) exit
      s = 2._R8P * s
   enddo
   ray%axis = 0_I4P
   if (count(ray%d /= 0._R8P) == 1) ray%axis = int(maxloc(abs(ray%d), dim=1), I4P)
   ray%pad = 1.e-9_R8P * r
   if (present(reach)) ray%pad = ray%pad + reach
   endfunction ray_prepare

   subroutine ray_facet_crossing(ray, facet, hit, t)
   !< Exact crossing of the prepared ray with a facet (issue #82, see the module header): hit and the parameter t >= 0.
   type(ray_exact_t),  intent(in)  :: ray    !< Prepared ray.
   type(facet_object), intent(in)  :: facet  !< Facet.
   logical,            intent(out) :: hit    !< The ray crosses the facet (t >= 0).
   real(R8P),          intent(out) :: t      !< Parameter of the crossing (0 if the origin is on the facet's plane).
   real(R8P)                       :: v(3, 3), p(3), w1, w2, w3, s1, s2, s3, tc, oo, qq, a, b !< Work.
   logical                         :: enter  !< Unused (line kernel's entering flag).
   integer(I4P)                    :: k      !< Axis.

   hit = .false. ; t = 0._R8P
   v(1, 1) = facet%vertex(1)%x ; v(2, 1) = facet%vertex(1)%y ; v(3, 1) = facet%vertex(1)%z
   v(1, 2) = facet%vertex(2)%x ; v(2, 2) = facet%vertex(2)%y ; v(3, 2) = facet%vertex(2)%z
   v(1, 3) = facet%vertex(3)%x ; v(2, 3) = facet%vertex(3)%y ; v(3, 3) = facet%vertex(3)%z
   k = ray%axis
   if (k > 0_I4P) then
      select case (k)
      case (1_I4P) ; a = ray%o(2) ; b = ray%o(3)
      case (2_I4P) ; a = ray%o(3) ; b = ray%o(1)
      case default ; a = ray%o(1) ; b = ray%o(2)
      endselect
      call tri_line_hit_oac(reshape(v, [9]), k, a, b, hit, tc, enter)
      if (.not. hit) return
      hit = .false.
   else
      call side(ray%o, ray%q, v(:, 2), v(:, 3), w1, s1)
      call side(ray%o, ray%q, v(:, 3), v(:, 1), w2, s2)
      if (s1 /= s2) return
      call side(ray%o, ray%q, v(:, 1), v(:, 2), w3, s3)
      if (s3 /= s1) return
      if (first_normal_sign(v) == 0._R8P) return ! collinear vertices: never crossed
   endif
   oo = o3d(v(:, 1), v(:, 2), v(:, 3), ray%o)
   if (oo == 0._R8P) then ! the origin on the facet's plane: the crossing is at o
      hit = .true. ; t = 0._R8P ; return
   endif
   qq = o3d(v(:, 1), v(:, 2), v(:, 3), ray%q)
   if (qq == 0._R8P .or. (oo > 0._R8P .eqv. qq > 0._R8P)) return ! behind the origin (q on the plane: never with a hit)
   hit = .true.
   if (k > 0_I4P) then
      t = (tc - ray%o(k)) / ray%d(k)
   else
      if (w1 + w2 + w3 /= 0._R8P) then
         p = (w1 * v(:, 1) + w2 * v(:, 2) + w3 * v(:, 3)) / (w1 + w2 + w3)
      else
         p = (v(:, 1) + v(:, 2) + v(:, 3)) / 3._R8P
      endif
      t = dot_product(p - ray%o, ray%d) / dot_product(ray%d, ray%d)
   endif
   t = max(t, 0._R8P)
   endsubroutine ray_facet_crossing

   subroutine side(o, q, x, y, w, sg)
   !< Plucker side test of the line o -> q against the edge x -> y: w = orient3d(o, q, x, y) (approximate magnitude,
   !< exact sign) and its sign under the ray's symbolic perturbation (module header): -1 or +1, never 0 for x /= y.
   real(R8P), intent(in)  :: o(3), q(3), x(3), y(3) !< Line, edge.
   real(R8P), intent(out) :: w                      !< orient3d(o, q, x, y).
   real(R8P), intent(out) :: sg                     !< Perturbed sign.
   real(R8P)              :: g(3), h(3), e(3)       !< g = n(y, q, x), g' = n(y, x, o), e = x - y (signs).
   integer(I4P)           :: k                      !< Counter.
   real(R8P), parameter   :: EPS_HALF = epsilon(1._R8P) / 2._R8P                !< Unit roundoff.
   real(R8P), parameter   :: O3D_BOUND = (7._R8P + 56._R8P * EPS_HALF) * EPS_HALF !< Filter bound (orient3d_xyz).
   real(R8P)              :: ad1, ad2, ad3, bd1, bd2, bd3, cd1, cd2, cd3, p1, p2, p3, p4, p5, p6 !< Filter.

   ! the filter of orient3d_xyz inlined (the same double when conclusive; the hot path of the ray queries)
   ad1 = o(1) - y(1) ; ad2 = o(2) - y(2) ; ad3 = o(3) - y(3)
   bd1 = q(1) - y(1) ; bd2 = q(2) - y(2) ; bd3 = q(3) - y(3)
   cd1 = x(1) - y(1) ; cd2 = x(2) - y(2) ; cd3 = x(3) - y(3)
   p1 = bd1 * cd2 ; p2 = cd1 * bd2 ; p3 = cd1 * ad2 ; p4 = ad1 * cd2 ; p5 = ad1 * bd2 ; p6 = bd1 * ad2
   w = ad3 * (p1 - p2) + bd3 * (p3 - p4) + cd3 * (p5 - p6)
   if (abs(w) > O3D_BOUND * ((abs(p1) + abs(p2)) * abs(ad3) + (abs(p3) + abs(p4)) * abs(bd3) + &
                            (abs(p5) + abs(p6)) * abs(cd3))) then
      sg = sign(1._R8P, w) ; return
   endif
   w = orient3d_xyz(o, q, x, y)
   if (w /= 0._R8P) then
      sg = sign(1._R8P, w) ; return
   endif
   do k = 1, 3
      g(k) = normal_sign(y, q, x, k)
      if (g(k) /= 0._R8P) then
         sg = g(k) ; return
      endif
   enddo
   do k = 1, 3
      e(k) = 0._R8P
      if (x(k) > y(k)) e(k) = 1._R8P
      if (x(k) < y(k)) e(k) = -1._R8P
   enddo
   h(1) = normal_sign(y, x, o, 1_I4P)
   sg = h(1)                 ; if (sg /= 0._R8P) return ! e^8
   sg = -e(3)                ; if (sg /= 0._R8P) return ! e^10
   sg = e(2)                 ; if (sg /= 0._R8P) return ! e^12
   sg = normal_sign(y, x, o, 2_I4P) ; if (sg /= 0._R8P) return ! e^16
   sg = -e(1)                                            ! e^20: x /= y, so e(1) /= 0 here
   endsubroutine side

   pure function o3d(a, b, c, d) result(det)
   !< orient3d_xyz(a, b, c, d) with its filter inlined: the same double when the filter is conclusive (the common case),
   !< else the exact `orient3d_xyz`. Exact sign either way; it spares the exact stage's large work array, which a debug
   !< build (-finit-real) initialises on every call.
   real(R8P), intent(in) :: a(3), b(3), c(3), d(3) !< Points.
   real(R8P)             :: det                    !< Its sign is exact.
   real(R8P), parameter  :: EPS_HALF = epsilon(1._R8P) / 2._R8P                !< Unit roundoff.
   real(R8P), parameter  :: O3D_BOUND = (7._R8P + 56._R8P * EPS_HALF) * EPS_HALF !< Filter bound (orient3d_xyz).
   real(R8P)             :: ad(3), bd(3), cd(3), bxcy, cxby, cxay, axcy, axby, bxay, perm !< Work.

   ad = a - d ; bd = b - d ; cd = c - d
   bxcy = bd(1) * cd(2) ; cxby = cd(1) * bd(2)
   cxay = cd(1) * ad(2) ; axcy = ad(1) * cd(2)
   axby = ad(1) * bd(2) ; bxay = bd(1) * ad(2)
   det = ad(3) * (bxcy - cxby) + bd(3) * (cxay - axcy) + cd(3) * (axby - bxay)
   perm = (abs(bxcy) + abs(cxby)) * abs(ad(3)) + (abs(cxay) + abs(axcy)) * abs(bd(3)) + &
          (abs(axby) + abs(bxay)) * abs(cd(3))
   if (abs(det) > O3D_BOUND * perm) return
   det = orient3d_xyz(a, b, c, d)
   endfunction o3d

   function normal_sign(t1, t2, t3, k) result(sg)
   !< Exact sign (-1, 0, +1) of component k of n = (t2 - t1) x (t3 - t1): the orientation of the projection on the
   !< plane of the other two coordinates (cyclic).
   real(R8P),    intent(in) :: t1(3), t2(3), t3(3) !< Triangle.
   integer(I4P), intent(in) :: k                   !< Component.
   real(R8P)                :: sg                  !< Sign.
   integer(I4P)             :: i1, i2              !< Projected coordinates.

   i1 = mod(k, 3_I4P) + 1_I4P ; i2 = mod(k + 1_I4P, 3_I4P) + 1_I4P
   sg = orient2d_xy(t1(i1), t1(i2), t2(i1), t2(i2), t3(i1), t3(i2))
   if (sg /= 0._R8P) sg = sign(1._R8P, sg)
   endfunction normal_sign

   subroutine ray_intersect_all_flat(facet, ray_origin, ray_direction, hits, status, exclude)
   !< Return all ray-facet intersections via flat O(n) loop, sorted ascending by `t`: one hit per crossing, exactly
   !< (issue #82, the module header), the facet `exclude` skipped.
   !<
   !< Hits with `t < 0` (behind the ray origin) are excluded; hits at t = 0 are kept.
   type(facet_object),           intent(in)               :: facet(:)        !< Facet array.
   type(vector_R8P),             intent(in)               :: ray_origin      !< Ray origin.
   type(vector_R8P),             intent(in)               :: ray_direction   !< Ray direction.
   type(ray_hit_t), allocatable, intent(out)              :: hits(:)         !< Sorted hit records.
   integer(I4P),                 intent(out),    optional :: status          !< Status code.
   integer(I4P),                 intent(in),     optional :: exclude         !< Facet to skip (issue #82).
   real(R8P)                                              :: t               !< Crossing parameter.
   logical                                                :: hit             !< Per-facet hit flag.
   real(R8P)                                              :: dir_norm2       !< |ray_direction|^2.
   integer(I4P)                                           :: nf, n_hits, f   !< Counters.
   type(ray_hit_t), allocatable                           :: tmp(:)          !< Pre-sort buffer.
   type(ray_exact_t)                                      :: ray             !< Prepared ray.
   type(vector_R8P)                                       :: bmin, bmax      !< Box of the facets.

   if (present(status)) status = RAY_STATUS_OK
   nf = size(facet, kind=I4P)
   if (nf == 0_I4P) then
      allocate(hits(0))
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif
   dir_norm2 = ray_direction%dotproduct(rhs=ray_direction)
   if (dir_norm2 < DIR_TOL) then
      allocate(hits(0))
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif

   call facets_box(facet, bmin, bmax)
   ray = ray_prepare(ray_origin, ray_direction, bmin, bmax)
   allocate(tmp(nf))   ! upper bound: every facet hit
   n_hits = 0_I4P
   do f = 1_I4P, nf
      if (present(exclude)) then
         if (f == exclude) cycle
      endif
      call ray_facet_crossing(ray, facet(f), hit, t)
      if (.not. hit) cycle
      n_hits = n_hits + 1_I4P
      tmp(n_hits)%facet_id = f
      tmp(n_hits)%t        = t
      tmp(n_hits)%point    = ray_origin + t * ray_direction
   enddo

   allocate(hits(n_hits))
   hits(1:n_hits) = tmp(1:n_hits)
   call sort_hits_by_t(hits)
   endsubroutine ray_intersect_all_flat

   pure subroutine facets_box(facet, bmin, bmax)
   !< Box of the vertices of the facets.
   type(facet_object), intent(in)  :: facet(:)   !< Facets.
   type(vector_R8P),   intent(out) :: bmin, bmax !< Box.
   integer(I4P)                    :: f, k       !< Counters.

   bmin = vector_R8P(huge(1._R8P), huge(1._R8P), huge(1._R8P)) ; bmax = -bmin
   do f = 1, size(facet, kind=I4P)
      do k = 1, 3
         bmin%x = min(bmin%x, facet(f)%vertex(k)%x) ; bmax%x = max(bmax%x, facet(f)%vertex(k)%x)
         bmin%y = min(bmin%y, facet(f)%vertex(k)%y) ; bmax%y = max(bmax%y, facet(f)%vertex(k)%y)
         bmin%z = min(bmin%z, facet(f)%vertex(k)%z) ; bmax%z = max(bmax%z, facet(f)%vertex(k)%z)
      enddo
   enddo
   endsubroutine facets_box

   pure subroutine sort_hits_by_t(hits)
   !< Insertion sort on `t` ascending, then on the facet id (issue #82). n is small (typically 2..10), so the
   !< O(n^2) constant beats the O(n log n) call overhead.
   type(ray_hit_t), intent(inout) :: hits(:)
   type(ray_hit_t)                :: key
   integer(I4P)                   :: i, j

   do i = 2_I4P, size(hits, kind=I4P)
      key = hits(i)
      j = i - 1_I4P
      do while (j >= 1_I4P)
         if (hits(j)%t < key%t .or. (hits(j)%t == key%t .and. hits(j)%facet_id <= key%facet_id)) exit ! t, then id
         hits(j + 1_I4P) = hits(j)
         j = j - 1_I4P
      enddo
      hits(j + 1_I4P) = key
   enddo
   endsubroutine sort_hits_by_t

endmodule fossil_ray_query
