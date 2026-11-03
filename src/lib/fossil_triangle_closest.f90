!< FOSSIL, point-to-triangle parametric solve (shared core).

module fossil_triangle_closest
!< Shared core of the point-to-triangle kernel: solve for the parametric
!< coordinates `(sq, tq)` of the closest point on a triangle to a query point.
!<
!< Extracted from `fossil_facet_object` (issue #20 Step 4) to work around an
!< nvfortran NVVM codegen bug (seen on 26.1, still present on 26.9): when
!< `triangle_closest_st` shares a TU
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
!<
!< The bug is not specific to this routine: adding ANY further
!< `!$acc routine seq` procedure that uses `vector_R8P` to `fossil_facet_object`
!< triggers it again (verified on 26.9). New device routines of that kind must
!< live in their own module, like this one.

use penf, only : I4P, R8P
use vecfor, only : vector_R8P, dotproduct_R8P_dev, vector_sub_vector_R8P_dev

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
   !< distance only, issue #19 §B5) are both thin wrappers over this.
   !<
   !< The closest point is `v1 + sq*e12 + tq*e13`; `(sq, tq)` are barycentric-style
   !< coordinates with V1 ↔ (0,0), V2 ↔ (1,0), V3 ↔ (0,1).
   !<
   !< Conditioning (issue #39). Eberly's solve goes through the Gram system: its barycentric
   !< error is ~eps max(a, c) / det, with det = |e12 x e13|^2. On needle triangles (aspect 1e3)
   !< that put the foot of a point 1e-10 L off a face on the edge, where the edge pseudo-normal
   !< is orthogonal to p - foot, and the sign flipped. Triangles with det < SLIVER max(a, c)^2
   !< (aspect beyond ~100) take `robust_st`, whose error is ~sqrt of that; the others keep the
   !< Eberly solve, which is ~20 % cheaper per query, bit for bit as before.
   type(vector_R8P), intent(in)  :: v1       !< Triangle first vertex.
   type(vector_R8P), intent(in)  :: e12      !< Edge 1-2, `v2 - v1`.
   type(vector_R8P), intent(in)  :: e13      !< Edge 1-3, `v3 - v1`.
   real(R8P),        intent(in)  :: a        !< `e12 . e12`.
   real(R8P),        intent(in)  :: b        !< `e12 . e13`.
   real(R8P),        intent(in)  :: c        !< `e13 . e13`.
   real(R8P),        intent(in)  :: det      !< Gram determinant `a*c - b*b`.
   type(vector_R8P), intent(in)  :: point    !< Point.
   real(R8P),        intent(out) :: sq       !< Parametric coordinate along `e12`.
   real(R8P),        intent(out) :: tq       !< Parametric coordinate along `e13`.
   real(R8P),        intent(out) :: distance !< Closest squared distance from point to the triangle.
   real(R8P), parameter          :: SLIVER = 1.e-4_R8P !< det / max(a, c)^2 below which the robust path runs.

   if (det >= SLIVER * max(a, c)**2) then
      call eberly_st(v1, e12, e13, a, b, c, det, point, sq, tq, distance)
   else
      call robust_st(v1, e12, e13, a, b, c, det, point, sq, tq, distance)
   endif
   endsubroutine triangle_closest_st

   pure subroutine eberly_st(v1, e12, e13, a, b, c, det, point, sq, tq, distance)
   !$acc routine seq
   !< Eberly's parametric solve (David Eberly, Geometric Tools): regions from the Gram system, fast and
   !< accurate on well-shaped triangles; its barycentric error grows like eps / |n|^2 (see triangle_closest_st).
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
   real(R8P)                     :: d, e, s, t                       !< Plane equation coefficients.
   real(R8P)                     :: tmp0, tmp1, numer, denom, invdet !< Temporary.
   real(R8P)                     :: rx, ry, rz                       !< Residual point - closest.

   ! Device-callable VecFor `_oac` API (issue #20 §Step 4): `vector - vector` and
   ! `.dot.` lower to TBP-generic operator dispatch, and nvfortran does not
   ! support type-bound procedure calls inside `!$acc routine seq`. The `_dev` free subroutines take `type(...)` (not
   ! `class`) and are device-safe. Bit-exact vs the operator form on CPU.
   call vector_sub_vector_R8P_dev(v1, point, V1P)
   d = dotproduct_R8P_dev(e12, V1P)
   e = dotproduct_R8P_dev(e13, V1P)
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
   ! Squared norm of the residual (v1 - point) + sq*e12 + tq*e13 (issue #39). The expanded
   ! quadratic form a sq^2 + 2b sq tq + c tq^2 + 2d sq + 2e tq + f, used before, cancels terms
   ! of size |v1 - point|^2: its absolute error is ~eps L^2, so it was noise (and could be
   ! negative, flipping the pseudo-normal sign of is_point_inside) below d ~ sqrt(eps) L ~ 1e-8 L.
   ! The residual form has an absolute error ~eps L in d: accurate down to d ~ 1e-14 L.
   rx = V1P%x + sq * e12%x + tq * e13%x
   ry = V1P%y + sq * e12%y + tq * e13%y
   rz = V1P%z + sq * e12%z + tq * e13%z
   distance = rx * rx + ry * ry + rz * rz
   endsubroutine eberly_st

   pure subroutine robust_st(v1, e12, e13, a, b, c, det, point, sq, tq, distance)
   !$acc routine seq
   !< Closest point by triple products with n = e12 x e13 (error ~eps / |n|), else the nearest clamped edge
   !< projection; for needle triangles (see triangle_closest_st). `b` and `det` are unused.
   type(vector_R8P), intent(in)  :: v1                               !< Triangle first vertex.
   type(vector_R8P), intent(in)  :: e12                              !< Edge 1-2, `v2 - v1`.
   type(vector_R8P), intent(in)  :: e13                              !< Edge 1-3, `v3 - v1`.
   real(R8P),        intent(in)  :: a                                !< `e12 . e12`.
   real(R8P),        intent(in)  :: b                                !< `e12 . e13` (unused).
   real(R8P),        intent(in)  :: c                                !< `e13 . e13`.
   real(R8P),        intent(in)  :: det                              !< Gram determinant `a*c - b*b` (unused).
   type(vector_R8P), intent(in)  :: point                            !< Point.
   real(R8P),        intent(out) :: sq                               !< Parametric coordinate along `e12`.
   real(R8P),        intent(out) :: tq                               !< Parametric coordinate along `e13`.
   real(R8P),        intent(out) :: distance                         !< Closest squared distance from point to the triangle.
   type(vector_R8P)              :: V1P                              !< `v1 - point`.
   real(R8P)                     :: nx, ny, nz, nn                   !< Normal e12 x e13, its squared norm.
   real(R8P)                     :: wx, wy, wz                       !< point - v1.
   real(R8P)                     :: cx, cy, cz                       !< Cross products.
   real(R8P)                     :: fx, fy, fz, ff, u                !< Edge 2-3, its squared length, edge parameter.
   real(R8P)                     :: s, t, d2                         !< Candidate coordinates, squared distance.

   ! Device-callable VecFor `_dev` API (issue #20 §Step 4): `vector - vector` lowers to
   ! TBP-generic operator dispatch, unsupported inside `!$acc routine seq`.
   call vector_sub_vector_R8P_dev(v1, point, V1P)
   wx = -V1P%x ; wy = -V1P%y ; wz = -V1P%z
   nx = e12%y * e13%z - e12%z * e13%y
   ny = e12%z * e13%x - e12%x * e13%z
   nz = e12%x * e13%y - e12%y * e13%x
   nn = nx * nx + ny * ny + nz * nz
   if (nn > 0._R8P) then ! interior foot by triple products
      cx = wy * e13%z - wz * e13%y ; cy = wz * e13%x - wx * e13%z ; cz = wx * e13%y - wy * e13%x
      s = (cx * nx + cy * ny + cz * nz) / nn
      cx = e12%y * wz - e12%z * wy ; cy = e12%z * wx - e12%x * wz ; cz = e12%x * wy - e12%y * wx
      t = (cx * nx + cy * ny + cz * nz) / nn
      if (s >= 0._R8P .and. t >= 0._R8P .and. s + t <= 1._R8P) then
         sq = s ; tq = t
         distance = residual2(V1P, e12, e13, sq, tq)
         return
      endif
   endif
   ! outside (or degenerate): nearest clamped projection on the three edges, first wins ties
   u = 0._R8P ; if (a > 0._R8P) u = min(max((wx * e12%x + wy * e12%y + wz * e12%z) / a, 0._R8P), 1._R8P)
   sq = u ; tq = 0._R8P ; distance = residual2(V1P, e12, e13, sq, tq)
   u = 0._R8P ; if (c > 0._R8P) u = min(max((wx * e13%x + wy * e13%y + wz * e13%z) / c, 0._R8P), 1._R8P)
   d2 = residual2(V1P, e12, e13, 0._R8P, u)
   if (d2 < distance) then
      sq = 0._R8P ; tq = u ; distance = d2
   endif
   fx = e13%x - e12%x ; fy = e13%y - e12%y ; fz = e13%z - e12%z
   ff = fx * fx + fy * fy + fz * fz
   u = 0._R8P
   if (ff > 0._R8P) u = min(max(((wx - e12%x) * fx + (wy - e12%y) * fy + (wz - e12%z) * fz) / ff, 0._R8P), 1._R8P)
   d2 = residual2(V1P, e12, e13, 1._R8P - u, u)
   if (d2 < distance) then
      sq = 1._R8P - u ; tq = u ; distance = d2
   endif
   endsubroutine robust_st

   pure function residual2(v1p, e12, e13, ss, tt) result(r2)
   !$acc routine seq
   !< Squared norm of the residual (v1 - point) + ss e12 + tt e13 (issue #39): accurate to ~eps L in the
   !< distance, unlike the expanded quadratic form (noise below ~1e-8 L).
   type(vector_R8P), intent(in) :: v1p        !< `v1 - point`.
   type(vector_R8P), intent(in) :: e12, e13   !< Edges.
   real(R8P),        intent(in) :: ss, tt     !< Coordinates.
   real(R8P)                    :: r2         !< Squared distance.
   real(R8P)                    :: rx, ry, rz !< Residual.

   rx = v1p%x + ss * e12%x + tt * e13%x
   ry = v1p%y + ss * e12%y + tt * e13%y
   rz = v1p%z + ss * e12%z + tt * e13%z
   r2 = rx * rx + ry * ry + rz * rz
   endfunction residual2
endmodule fossil_triangle_closest
