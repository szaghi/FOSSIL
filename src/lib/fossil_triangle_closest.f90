!< FOSSIL, point-to-triangle parametric solve (shared core).

module fossil_triangle_closest
!< Shared core of the point-to-triangle kernel: solve for the parametric
!< coordinates `(sq, tq)` of the closest point on a triangle to a query point.
!<
!< Extracted from `fossil_facet_object` (issue #20 Step 4) to work around
!< nvfortran 26.1 NVVM codegen bug: when `triangle_closest_st` shares a TU
!< with `triangle_point_distance` and `triangle_point_distance_sq` (both
!< `!$acc routine seq` callers), the device-IR emission emits two
!< structurally-identical-but-nominally-distinct LLVM type aliases
!< (`%struct.DT1_*` and `%struct.DT7_*`) for `vector_R8P`, then declares the
!< VecFor common-block global once with each alias — which NVVM rejects:
!< `parse forward reference and definition of global have different types`.
!<
!< Isolating `triangle_closest_st` in its own module gives nvfortran a clean
!< compilation context where only ONE LLVM type alias is generated for
!< `vector_R8P`, and the duplicate-global emission is avoided.

use penf, only : I4P, R8P
use vecfor, only : vector_R8P, dotproduct_R8P_oac, vector_sub_vector_R8P_oac

implicit none
private
public :: triangle_closest_st

contains
   pure subroutine triangle_closest_st(v1, e12, e13, a, b, c, det, point, sq, tq, distance)
   !$acc routine seq
   !< Solve for the parametric coordinates `(sq, tq)` of the closest point on a
   !< triangle to `point`, and the squared distance there.
   !<
   !< The shared core of the point-to-triangle kernel — `triangle_point_distance`
   !< (rich: + closest point + region) and `triangle_point_distance_sq` (lean:
   !< distance only, issue #19 §B5) are both thin wrappers over this. Keeping the
   !< ~120-line branchy `(sq, tq)` solve in one place is the single source of truth:
   !< the two public faces differ only in what they compute *after* the solve.
   !<
   !< The closest point is `v1 + sq*e12 + tq*e13`; `(sq, tq)` are barycentric-style
   !< coordinates with V1 ↔ (0,0), V2 ↔ (1,0), V3 ↔ (0,1).
   !<
   !< @note Algorithm by David Eberly, Geometric Tools LLC, http://www.geometrictools.com.
   type(vector_R8P), intent(in)  :: v1                               !< Triangle first vertex.
   type(vector_R8P), intent(in)  :: e12                              !< Edge 1-2, `v2 - v1`.
   type(vector_R8P), intent(in)  :: e13                              !< Edge 1-3, `v3 - v1`.
   real(R8P),        intent(in)  :: a                                !< `e12 . e12`.
   real(R8P),        intent(in)  :: b                                !< `e12 . e13`.
   real(R8P),        intent(in)  :: c                                !< `e13 . e13`.
   real(R8P),        intent(in)  :: det                              !< Gram determinant `a*c - b*b`.
   type(vector_R8P), intent(in)  :: point                            !< Point.
   real(R8P),        intent(out) :: sq                               !< Parametric coordinate along `e12`.
   real(R8P),        intent(out) :: tq                               !< Parametric coordinate along `e13`.
   real(R8P),        intent(out) :: distance                         !< Closest squared distance from point to the triangle.
   type(vector_R8P)              :: V1P                              !< `v1 - point`.
   real(R8P)                     :: d, e, f, s, t                    !< Plane equation coefficients.
   real(R8P)                     :: tmp0, tmp1, numer, denom, invdet !< Temporary.

   ! Device-callable VecFor `_oac` API (issue #20 §Step 4): `vector - vector` and
   ! `.dot.` lower to TBP-generic operator dispatch, which nvfortran 26.1 rejects
   ! inside `!$acc routine seq`. The `_oac` free subroutines take `type(...)` (not
   ! `class`) and are device-safe. Bit-exact vs the operator form on CPU.
   call vector_sub_vector_R8P_oac(v1, point, V1P)
   d = dotproduct_R8P_oac(e12, V1P)
   e = dotproduct_R8P_oac(e13, V1P)
   f = dotproduct_R8P_oac(V1P, V1P)
   s = b * e - c * d
   t = b * d - a * e
   if (s+t <= det) then
      if (s < 0._R8P) then
         if (t < 0._R8P) then ! region 4
            if (e < 0._R8P) then
               sq = 0._R8P
               if (c >= -e) then
                  tq = -e / c
               else
                  tq = 1._R8P
               endif
            else
               if (d > 0._R8P) then
                  sq = 0._R8P
               else
                  if (a >= -d) then
                     sq = -d / a
                  else
                     sq = 1._R8P
                  endif
               endif
               tq = 0._R8P
            endif
         else ! region 3
            sq = 0._R8P
            if (e >= 0._R8P) then
               tq = 0._R8P
            else
               if (-e >= c) then
                  tq = 1._R8P
               else
                  tq = -e / c
               endif
            endif
         endif
      elseif (t < 0._R8P) then ! region 5
         if (d >= 0._R8P) then
            sq = 0._R8P
         else
            if (-d >= a) then
               sq = 1._R8P
            else
               sq = -d / a
            endif
         endif
         tq = 0._R8P
      else ! region 0
        invdet = 1._R8P / det
        sq = s * invdet
        tq = t * invdet
      endif
   else
      if (s < 0._R8P) then ! region 2
         tmp0 = b + d
         tmp1 = c + e
         if (tmp1 > tmp0) then
            numer = tmp1 - tmp0
            denom = a - 2._R8P * b + c
            if (numer >= denom) then
               sq = 1._R8P
            else
               sq = numer / denom
            endif
            tq = 1._R8P - sq
         else
            sq = 0._R8P
            if (tmp1 <= 0._R8P) then
               tq = 1._R8P
            else
               if (e >= 0._R8P) then
                  tq = 0._R8P
               else
                  tq = -e / c
               endif
            endif
         endif
      elseif (t < 0._R8P) then ! region 6
         tmp0 = a + d
         tmp1 = b + e
         if (tmp0 > tmp1) then
            numer = b + d - c - e
            denom = a - 2._R8P * b + c
            if (numer >= 0._R8P) then
               sq = 0._R8P
            else
               if (denom > -numer) then
                  sq = -numer / denom
               else
                  sq = 1._R8P
               endif
            endif
            tq = 1._R8P - sq
         else
            if (tmp0 <= 0._R8P) then
               sq = 1._R8P
            else
               if (d >= 0._R8P) then
                  sq = 0._R8P
               else
                  sq = -d / a
               endif
            endif
            tq = 0._R8P
         endif
      else ! region 1
         numer = c + e - b - d
         if (numer <= 0._R8P) then
            sq = 0._R8P
         else
            denom = a - 2._R8P * b + c
            if (numer >= denom) then
               sq = 1._R8P
            else
               sq = numer / denom
            endif
         endif
         tq = 1._R8P - sq
      endif
   endif
   ! Quadratic form Q(sq, tq) = ||point - (v1 + sq*e12 + tq*e13)||^2. Non-negative
   ! by construction (it is a squared norm); no defensive abs() needed — issue #19
   ! §B5 dropped the abs that the original Eberly port carried.
   distance = a * sq * sq + 2._R8P * b * sq * tq + c * tq * tq + 2._R8P * d * sq + 2._R8P * e * tq + f
   endsubroutine triangle_closest_st
endmodule fossil_triangle_closest
