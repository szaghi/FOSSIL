!< FOSSIL, the exact crossing of an axis-aligned line with a triangle (issue #40), shared by the line queries and the
!< ray queries (issue #82).
module fossil_line_kernel
!< FOSSIL, the exact crossing of an axis-aligned line with a triangle (issue #40), shared by the line queries and the
!< ray queries (issue #82).
!<
!< The rule (one global symbolic perturbation of the line, issue #69) is described in the header of
!< `fossil_line_query`, which re-exports `tri_line_hit_oac`. It lives here, below the AABB tree and the ray queries, so
!< that `intersect_ray_*` decide an axis-aligned ray with the very same kernel as `line_crossings` and the cut cells.
use fossil_predicates, only : orient2d_xy
use penf,              only : I4P, R8P

implicit none
private
public :: tri_line_hit_oac, first_normal_sign

contains
   pure subroutine tri_line_hit_oac(tv, axis, a, b, hit, t, enter)
   !< Crossing of the line along `axis` at (a, b) with the triangle tv = [A, B, C] (9 coordinates).
   !$acc routine seq
   real(R8P),    intent(in)  :: tv(9)  !< Triangle vertices (x, y, z of A, then B, then C).
   integer(I4P), intent(in)  :: axis   !< Line axis (1, 2, 3).
   real(R8P),    intent(in)  :: a, b   !< Line position (the other two coordinates, cyclic).
   logical,      intent(out) :: hit    !< The line crosses the triangle.
   real(R8P),    intent(out) :: t      !< Crossing coordinate along the axis.
   logical,      intent(out) :: enter  !< The line (in +axis) enters the solid there.
   integer(I4P)              :: i1, i2 !< Projected coordinate indices.
   real(R8P)                 :: wa, wb, wc, sa, sb, sc !< Weights and perturbed signs.
   logical                   :: bfirst !< Perturb b before a.

   i1 = mod(axis, 3_I4P) + 1_I4P ; i2 = mod(axis + 1_I4P, 3_I4P) + 1_I4P
   hit = .false. ; t = 0._R8P ; enter = .false.
   bfirst = axis == 2_I4P ! global perturbation x + eps1, y + eps2, z + eps3: on a y line (a, b) = (z, x), x first
   call psign(tv(3 + i1), tv(3 + i2), tv(6 + i1), tv(6 + i2), a, b, bfirst, wa, sa) ! edge B -> C: weight of A
   call psign(tv(6 + i1), tv(6 + i2), tv(i1), tv(i2), a, b, bfirst, wb, sb)         ! edge C -> A: weight of B
   call psign(tv(i1), tv(i2), tv(3 + i1), tv(3 + i2), a, b, bfirst, wc, sc)         ! edge A -> B: weight of C
   if (sa == 0._R8P .or. sa /= sb .or. sa /= sc) return
   hit = .true.
   enter = sa < 0._R8P
   if (wa + wb + wc /= 0._R8P) then
      t = (wa * tv(axis) + wb * tv(3 + axis) + wc * tv(6 + axis)) / (wa + wb + wc)
   else
      t = (tv(axis) + tv(3 + axis) + tv(6 + axis)) / 3._R8P
   endif
   endsubroutine tri_line_hit_oac

   pure subroutine psign(p1, p2, q1, q2, a, b, bfirst, w, s)
   !< Weight orient2d(p, q, (a, b)) (approximate magnitude, exact sign) and its perturbed sign: (a, b) + (eps, eps^2),
   !< or (a, b) + (eps^2, eps) with `bfirst` (issue #69: the two orders make one global perturbation of the grid).
   !$acc routine seq
   real(R8P), intent(in)  :: p1, p2, q1, q2 !< Edge end points (projected).
   real(R8P), intent(in)  :: a, b           !< Line position.
   logical,   intent(in)  :: bfirst         !< Perturb b first.
   real(R8P), intent(out) :: w, s           !< Weight, perturbed sign (-1, 0, +1).
   real(R8P)              :: da, db         !< d/da = p2 - q2, d/db = q1 - p1 of orient2d(p, q, P).

   w = orient2d_xy(p1, p2, q1, q2, a, b)
   da = 0._R8P ; db = 0._R8P
   if (p2 > q2) da = 1._R8P
   if (p2 < q2) da = -1._R8P
   if (q1 > p1) db = 1._R8P
   if (q1 < p1) db = -1._R8P
   if (w > 0._R8P) then
      s = 1._R8P
   elseif (w < 0._R8P) then
      s = -1._R8P
   elseif (bfirst) then
      s = db ; if (s == 0._R8P) s = da
   else
      s = da ; if (s == 0._R8P) s = db
   endif
   endsubroutine psign

   pure function first_normal_sign(tri) result(sg)
   !< Sign (-1, 0, +1) of the first nonzero component of n = (B - A) x (C - A) in the order x, y, z, exactly (the
   !< components are the orientations of the projections on (y, z), (z, x), (x, y)). A device routine (issue #100).
   !$acc routine seq
   real(R8P), intent(in) :: tri(3, 3) !< Triangle vertices (columns).
   real(R8P)             :: sg
   integer(I4P)          :: k, i1, i2

   sg = 0._R8P
   do k = 1, 3
      i1 = mod(k, 3) + 1 ; i2 = mod(k + 1, 3) + 1
      sg = orient2d_xy(tri(i1, 1), tri(i2, 1), tri(i1, 2), tri(i2, 2), tri(i1, 3), tri(i2, 3))
      if (sg /= 0._R8P) then
         sg = sign(1._R8P, sg) ; return
      endif
   enddo
   endfunction first_normal_sign
endmodule fossil_line_kernel
