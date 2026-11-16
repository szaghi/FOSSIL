!< FOSSIL, the exact ray-facet crossing (issue #82) on plain arrays: the kernel of the host ray queries and of the segment
!< visibility of the probes, CPU and device (WP1 step 01).
module fossil_ray_kernel
!< FOSSIL, the exact ray-facet crossing (issue #82) on plain arrays: the kernel of the host ray queries and of the segment
!< visibility of the probes, CPU and device (WP1 step 01).
!<
!< The rule is described in the header of `fossil_ray_query`: the ray o + t d is the oriented line from o to the far point
!< q (beyond the surface box, rounded once by `ray_prepare_oac`); axis-aligned rays go through `tri_line_hit_oac` (the
!< line kernel of `line_crossings` and the cut cells), the others through three Plucker side tests made exact (orient3d
!< with one symbolic perturbation of the ray); the crossing is ahead iff o and q lie on opposite sides of the facet's
!< plane, an origin on the plane is a hit at t = 0.
!<
!< Every routine is `!$acc routine seq` on plain reals (no vector_R8P, no class, no allocation, no internal procedure):
!< `fossil_ray_query` (host `intersect_ray_*`) and `fossil_segment_visibility` (host and device probes) run this very
!< code, so the hits agree bit for bit. Small reductions are written as scalar expressions in the order gfortran
!< evaluates the array forms they replaced (`sum`, `dot_product`: from 0, left to right), so the host ray queries are
!< unchanged bit for bit (nvfortran 26.9 got an array-section `maxval` wrong in the host code of a routine-seq kernel).
use fossil_line_kernel, only : tri_line_hit_oac, first_normal_sign
use fossil_predicates,  only : orient2d_xy, orient3d_xyz
use penf,               only : I4P, R8P

implicit none
private
public :: ray_prepare_oac, ray_crossing_oac

contains
   pure subroutine ray_prepare_oac(o, d, lo, hi, reach, q, axis, pad)
   !< The far point q = o + s d beyond the box [lo, hi] (the surface box), s = 2 R / |d| with R >= the distance from o to
   !< any point of the box (s doubled until q /= o); the axis of an axis-aligned ray (two zero components of d; else 0)
   !< and the pad of the slab test, 1e-9 R + reach (issue #82).
   !$acc routine seq
   real(R8P),    intent(in)  :: o(3)          !< Origin.
   real(R8P),    intent(in)  :: d(3)          !< Direction (nonzero).
   real(R8P),    intent(in)  :: lo(3), hi(3)  !< Box of the surface.
   real(R8P),    intent(in)  :: reach         !< Facet reach beyond its tree box (octree; 0 for the BVH).
   real(R8P),    intent(out) :: q(3)          !< Far point.
   integer(I4P), intent(out) :: axis          !< 1, 2, 3 for an axis-aligned ray, else 0.
   real(R8P),    intent(out) :: pad           !< Pad of the boxes in the slab test.
   real(R8P)                 :: m1, m2, m3, r, s, dn, mo !< Work.
   integer(I4P)              :: k, nz         !< Counters.

   m1 = max(abs(lo(1) - o(1)), abs(hi(1) - o(1)))
   m2 = max(abs(lo(2) - o(2)), abs(hi(2) - o(2)))
   m3 = max(abs(lo(3) - o(3)), abs(hi(3) - o(3)))
   r = sqrt(((0._R8P + m1 * m1) + m2 * m2) + m3 * m3)
   mo = max(abs(o(1)), abs(o(2)), abs(o(3)))
   r = r + 1.e-6_R8P * mo + tiny(1._R8P)
   dn = sqrt(((0._R8P + d(1) * d(1)) + d(2) * d(2)) + d(3) * d(3))
   s = 2._R8P * r / dn
   do k = 1, 64
      q(1) = o(1) + s * d(1) ; q(2) = o(2) + s * d(2) ; q(3) = o(3) + s * d(3)
      if (q(1) /= o(1) .or. q(2) /= o(2) .or. q(3) /= o(3)) exit
      s = 2._R8P * s
   enddo
   nz = 0_I4P ; axis = 0_I4P
   do k = 1, 3
      if (d(k) /= 0._R8P) then
         nz = nz + 1_I4P ; axis = int(k, I4P)
      endif
   enddo
   if (nz /= 1_I4P) axis = 0_I4P
   pad = 1.e-9_R8P * r + reach
   endsubroutine ray_prepare_oac

   pure subroutine ray_crossing_oac(o, q, d, axis, tv, hit, t)
   !< Exact crossing of the prepared ray (o, q, d, axis) with the triangle tv = [A, B, C] (9 coordinates): hit and the
   !< parameter t >= 0 (0 if the origin is on the triangle's plane). A triangle with collinear vertices is never crossed.
   !$acc routine seq
   real(R8P),    intent(in)  :: o(3), q(3), d(3) !< Origin, far point, direction.
   integer(I4P), intent(in)  :: axis   !< Axis of an axis-aligned ray, else 0.
   real(R8P),    intent(in)  :: tv(9)  !< Triangle vertices (x, y, z of A, then B, then C).
   logical,      intent(out) :: hit    !< The ray crosses the triangle (t >= 0).
   real(R8P),    intent(out) :: t      !< Parameter of the crossing.
   real(R8P)                 :: v(3, 3), w1, w2, w3, s1, s2, s3, tc, oo, qq, a, b, p1, p2, p3, ws !< Work.
   logical                   :: enter  !< Unused (line kernel's entering flag).
   integer(I4P)              :: k      !< Counter.

   hit = .false. ; t = 0._R8P
   do k = 1, 3
      v(k, 1) = tv(k) ; v(k, 2) = tv(3 + k) ; v(k, 3) = tv(6 + k)
   enddo
   w1 = 0._R8P ; w2 = 0._R8P ; w3 = 0._R8P ; tc = 0._R8P
   if (axis > 0_I4P) then
      select case (axis)
      case (1_I4P) ; a = o(2) ; b = o(3)
      case (2_I4P) ; a = o(3) ; b = o(1)
      case default ; a = o(1) ; b = o(2)
      endselect
      call tri_line_hit_oac(tv, axis, a, b, hit, tc, enter)
      if (.not. hit) return
      hit = .false.
   else
      call side_oac(o, q, v(:, 2), v(:, 3), w1, s1)
      call side_oac(o, q, v(:, 3), v(:, 1), w2, s2)
      if (s1 /= s2) return
      call side_oac(o, q, v(:, 1), v(:, 2), w3, s3)
      if (s3 /= s1) return
      if (first_normal_sign(v) == 0._R8P) return ! collinear vertices: never crossed
   endif
   oo = o3d_oac(v(:, 1), v(:, 2), v(:, 3), o)
   if (oo == 0._R8P) then ! the origin on the triangle's plane: the crossing is at o
      hit = .true. ; t = 0._R8P ; return
   endif
   qq = o3d_oac(v(:, 1), v(:, 2), v(:, 3), q)
   if (qq == 0._R8P .or. (oo > 0._R8P .eqv. qq > 0._R8P)) return ! behind the origin (q on the plane: never with a hit)
   hit = .true.
   if (axis > 0_I4P) then
      t = (tc - o(axis)) / d(axis)
   else
      ws = w1 + w2 + w3
      if (ws /= 0._R8P) then
         p1 = (w1 * v(1, 1) + w2 * v(1, 2) + w3 * v(1, 3)) / ws
         p2 = (w1 * v(2, 1) + w2 * v(2, 2) + w3 * v(2, 3)) / ws
         p3 = (w1 * v(3, 1) + w2 * v(3, 2) + w3 * v(3, 3)) / ws
      else
         p1 = (v(1, 1) + v(1, 2) + v(1, 3)) / 3._R8P
         p2 = (v(2, 1) + v(2, 2) + v(2, 3)) / 3._R8P
         p3 = (v(3, 1) + v(3, 2) + v(3, 3)) / 3._R8P
      endif
      t = (((0._R8P + (p1 - o(1)) * d(1)) + (p2 - o(2)) * d(2)) + (p3 - o(3)) * d(3)) / &
          (((0._R8P + d(1) * d(1)) + d(2) * d(2)) + d(3) * d(3))
   endif
   t = max(t, 0._R8P)
   endsubroutine ray_crossing_oac

   pure subroutine side_oac(o, q, x, y, w, sg)
   !< Plucker side test of the line o -> q against the edge x -> y: w = orient3d(o, q, x, y) (approximate magnitude,
   !< exact sign) and its sign under the ray's symbolic perturbation (`fossil_ray_query` header): -1 or +1, never 0 for
   !< x /= y.
   !$acc routine seq
   real(R8P), intent(in)  :: o(3), q(3), x(3), y(3) !< Line, edge.
   real(R8P), intent(out) :: w                      !< orient3d(o, q, x, y).
   real(R8P), intent(out) :: sg                     !< Perturbed sign.
   real(R8P)              :: g, e1, e2, e3          !< A component of n(y, q, x); signs of e = x - y.
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
      g = normal_sign_oac(y, q, x, k)
      if (g /= 0._R8P) then
         sg = g ; return
      endif
   enddo
   e1 = 0._R8P ; e2 = 0._R8P ; e3 = 0._R8P
   if (x(1) > y(1)) e1 = 1._R8P
   if (x(1) < y(1)) e1 = -1._R8P
   if (x(2) > y(2)) e2 = 1._R8P
   if (x(2) < y(2)) e2 = -1._R8P
   if (x(3) > y(3)) e3 = 1._R8P
   if (x(3) < y(3)) e3 = -1._R8P
   sg = normal_sign_oac(y, x, o, 1_I4P) ; if (sg /= 0._R8P) return ! e^8
   sg = -e3                                    ; if (sg /= 0._R8P) return ! e^10
   sg = e2                                     ; if (sg /= 0._R8P) return ! e^12
   sg = normal_sign_oac(y, x, o, 2_I4P) ; if (sg /= 0._R8P) return ! e^16
   sg = -e1                                                               ! e^20: x /= y, so e1 /= 0 here
   endsubroutine side_oac

   pure function o3d_oac(a, b, c, d) result(det)
   !< orient3d_xyz(a, b, c, d) with its filter inlined: the same double when the filter is conclusive (the common case),
   !< else the exact `orient3d_xyz`. Exact sign either way; it spares the exact stage's large work array, which a debug
   !< build (-finit-real) initialises on every call.
   !$acc routine seq
   real(R8P), intent(in) :: a(3), b(3), c(3), d(3) !< Points.
   real(R8P)             :: det                    !< Its sign is exact.
   real(R8P), parameter  :: EPS_HALF = epsilon(1._R8P) / 2._R8P                !< Unit roundoff.
   real(R8P), parameter  :: O3D_BOUND = (7._R8P + 56._R8P * EPS_HALF) * EPS_HALF !< Filter bound (orient3d_xyz).
   real(R8P)             :: ad1, ad2, ad3, bd1, bd2, bd3, cd1, cd2, cd3 !< Differences.
   real(R8P)             :: bxcy, cxby, cxay, axcy, axby, bxay, perm   !< Work.

   ad1 = a(1) - d(1) ; ad2 = a(2) - d(2) ; ad3 = a(3) - d(3)
   bd1 = b(1) - d(1) ; bd2 = b(2) - d(2) ; bd3 = b(3) - d(3)
   cd1 = c(1) - d(1) ; cd2 = c(2) - d(2) ; cd3 = c(3) - d(3)
   bxcy = bd1 * cd2 ; cxby = cd1 * bd2
   cxay = cd1 * ad2 ; axcy = ad1 * cd2
   axby = ad1 * bd2 ; bxay = bd1 * ad2
   det = ad3 * (bxcy - cxby) + bd3 * (cxay - axcy) + cd3 * (axby - bxay)
   perm = (abs(bxcy) + abs(cxby)) * abs(ad3) + (abs(cxay) + abs(axcy)) * abs(bd3) + &
          (abs(axby) + abs(bxay)) * abs(cd3)
   if (abs(det) > O3D_BOUND * perm) return
   det = orient3d_xyz(a, b, c, d)
   endfunction o3d_oac

   pure function normal_sign_oac(t1, t2, t3, k) result(sg)
   !< Exact sign (-1, 0, +1) of component k of n = (t2 - t1) x (t3 - t1): the orientation of the projection on the
   !< plane of the other two coordinates (cyclic).
   !$acc routine seq
   real(R8P),    intent(in) :: t1(3), t2(3), t3(3) !< Triangle.
   integer(I4P), intent(in) :: k                   !< Component.
   real(R8P)                :: sg                  !< Sign.
   integer(I4P)             :: i1, i2              !< Projected coordinates.

   i1 = mod(k, 3_I4P) + 1_I4P ; i2 = mod(k + 1_I4P, 3_I4P) + 1_I4P
   sg = orient2d_xy(t1(i1), t1(i2), t2(i1), t2(i2), t3(i1), t3(i2))
   if (sg /= 0._R8P) sg = sign(1._R8P, sg)
   endfunction normal_sign_oac
endmodule fossil_ray_kernel
