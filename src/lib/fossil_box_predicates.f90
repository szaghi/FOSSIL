!< FOSSIL, box predicates for AMR tagging: box-triangle (3-D) and box-segment (2-D) (issue #38).
module fossil_box_predicates
!< FOSSIL, box predicates for AMR tagging: box-triangle (3-D) and box-segment (2-D) (issue #38).
!<
!< Plain-real, pure procedures (device-callable, `!$acc routine seq`).
!<
!< 3-D, `box_triangle_overlap`: the separating-axis test of Akenine-Moeller (2001) on the closed
!< box [lo, hi]: the 3 box normals, the triangle normal and the 9 cross products of box axes with
!< triangle edges. It is evaluated in floating point on the box inflated outward by `inflate`
!< (callers pass ~1e-12 of the problem size), so it is **conservative**: it never misses a box that
!< touches the triangle, and may report a box closer than `inflate` as touching. For AMR tagging a
!< false positive costs one refined cell, a false negative loses a cut cell.
!<
!< 2-D, `box_segment_overlap_2d`: **exact**. A segment meets the closed box iff an endpoint lies in
!< it or the segment shares a point with one of the four box edges, decided with `orient2d_xy`.
!< `box_segment_distance2_2d` is the squared distance (0 exactly when they overlap).
use fossil_predicates, only : orient2d_xy
use penf,              only : I4P, R8P

implicit none
private
public :: box_triangle_overlap
public :: box_segment_overlap_2d, box_segment_distance2_2d, segments_intersect_2d

contains
   pure function box_triangle_overlap(lo, hi, v1, v2, v3, inflate) result(overlap)
   !< True if the closed box [lo - inflate, hi + inflate] and the triangle (v1, v2, v3) share a point
   !< (separating-axis test, conservative by the inflation).
   !$acc routine seq
   real(R8P), intent(in) :: lo(3), hi(3)     !< Box.
   real(R8P), intent(in) :: v1(3), v2(3), v3(3) !< Triangle.
   real(R8P), intent(in) :: inflate          !< Outward inflation (>= 0).
   logical               :: overlap          !< Result.
   real(R8P)             :: c(3), h(3), a(3), b(3), d(3), e(3, 3), n(3), ax(3), p0, p1, p2, r, u(3) !< Work.
   integer(I4P)          :: i, j             !< Counters.

   c = 0.5_R8P * (lo + hi)
   h = 0.5_R8P * (hi - lo) + inflate
   a = v1 - c ; b = v2 - c ; d = v3 - c
   overlap = .false.
   ! box face normals: triangle AABB against the box
   do i = 1, 3
      if (min(a(i), b(i), d(i)) > h(i) .or. max(a(i), b(i), d(i)) < -h(i)) return
   enddo
   ! triangle normal: plane against the box
   e(:, 1) = b - a ; e(:, 2) = d - b ; e(:, 3) = a - d
   n = cross(e(:, 1), e(:, 2))
   r = h(1) * abs(n(1)) + h(2) * abs(n(2)) + h(3) * abs(n(3))
   if (abs(dot_product(n, a)) > r) return
   ! the 9 cross products of the box axes with the triangle edges
   do i = 1, 3       ! box axis
      u = 0._R8P ; u(i) = 1._R8P
      do j = 1, 3    ! triangle edge
         ax = cross(u, e(:, j))
         if (all(ax == 0._R8P)) cycle ! edge parallel to the axis: covered by the face tests
         p0 = dot_product(ax, a) ; p1 = dot_product(ax, b) ; p2 = dot_product(ax, d)
         r = h(1) * abs(ax(1)) + h(2) * abs(ax(2)) + h(3) * abs(ax(3))
         if (min(p0, p1, p2) > r .or. max(p0, p1, p2) < -r) return
      enddo
   enddo
   overlap = .true.
   endfunction box_triangle_overlap

   pure function box_segment_overlap_2d(lo, hi, a, b) result(overlap)
   !< True if the segment ab and the closed box [lo, hi] share a point (exact).
   !$acc routine seq
   real(R8P), intent(in) :: lo(2), hi(2) !< Box.
   real(R8P), intent(in) :: a(2), b(2)   !< Segment.
   logical               :: overlap      !< Result.

   overlap = in_box(a) .or. in_box(b)
   if (overlap) return
   ! quick reject: segment box disjoint from the box
   if (max(a(1), b(1)) < lo(1) .or. min(a(1), b(1)) > hi(1) .or. max(a(2), b(2)) < lo(2) .or. min(a(2), b(2)) > hi(2)) return
   overlap = segments_intersect_2d(a(1), a(2), b(1), b(2), lo(1), lo(2), hi(1), lo(2)) .or. &
             segments_intersect_2d(a(1), a(2), b(1), b(2), hi(1), lo(2), hi(1), hi(2)) .or. &
             segments_intersect_2d(a(1), a(2), b(1), b(2), hi(1), hi(2), lo(1), hi(2)) .or. &
             segments_intersect_2d(a(1), a(2), b(1), b(2), lo(1), hi(2), lo(1), lo(2))
   contains
      pure function in_box(p) result(yes)
      !< p in the closed box.
      real(R8P), intent(in) :: p(2) !< Point.
      logical               :: yes  !< Result.

      yes = p(1) >= lo(1) .and. p(1) <= hi(1) .and. p(2) >= lo(2) .and. p(2) <= hi(2)
      endfunction in_box
   endfunction box_segment_overlap_2d

   pure function box_segment_distance2_2d(lo, hi, a, b) result(d2)
   !< Squared distance between the segment ab and the closed box [lo, hi]; exactly 0 when they
   !< overlap (exact test). Otherwise the minimum over the segment ends against the box and the box
   !< corners against the segment (the closest pair of a convex polygon and a segment that do not
   !< meet involves a vertex of one of them).
   !$acc routine seq
   real(R8P), intent(in) :: lo(2), hi(2) !< Box.
   real(R8P), intent(in) :: a(2), b(2)   !< Segment.
   real(R8P)             :: d2           !< Squared distance.
   real(R8P)             :: corner(2, 4) !< Box corners.
   integer(I4P)          :: k            !< Counter.

   d2 = 0._R8P
   if (box_segment_overlap_2d(lo, hi, a, b)) return
   d2 = min(point_box_d2(a), point_box_d2(b))
   corner = reshape([lo(1), lo(2), hi(1), lo(2), hi(1), hi(2), lo(1), hi(2)], [2, 4])
   do k = 1, 4
      d2 = min(d2, point_segment_d2(corner(:, k)))
   enddo
   contains
      pure function point_box_d2(p) result(dd)
      !< Squared distance from p to the box.
      real(R8P), intent(in) :: p(2) !< Point.
      real(R8P)             :: dd   !< Squared distance.
      real(R8P)             :: e(2) !< Excess.

      e = max(lo - p, 0._R8P, p - hi)
      dd = e(1)**2 + e(2)**2
      endfunction point_box_d2

      pure function point_segment_d2(p) result(dd)
      !< Squared distance from p to the segment ab.
      real(R8P), intent(in) :: p(2) !< Point.
      real(R8P)             :: dd   !< Squared distance.
      real(R8P)             :: w(2), s(2), t, l2 !< Work.

      s = b - a ; w = p - a ; l2 = s(1)**2 + s(2)**2
      t = 0._R8P
      if (l2 > 0._R8P) t = max(0._R8P, min(1._R8P, (w(1) * s(1) + w(2) * s(2)) / l2))
      w = p - (a + t * s)
      dd = w(1)**2 + w(2)**2
      endfunction point_segment_d2
   endfunction box_segment_distance2_2d

   pure function segments_intersect_2d(p1u, p1v, p2u, p2v, q1u, q1v, q2u, q2v) result(hit)
   !< True if the closed segments p1p2 and q1q2 share a point (exact, `orient2d_xy`).
   !$acc routine seq
   real(R8P), intent(in) :: p1u, p1v, p2u, p2v, q1u, q1v, q2u, q2v !< Endpoints.
   logical               :: hit                                     !< Result.
   real(R8P)             :: o1, o2, o3, o4                           !< Orientations.

   o1 = orient2d_xy(p1u, p1v, p2u, p2v, q1u, q1v)
   o2 = orient2d_xy(p1u, p1v, p2u, p2v, q2u, q2v)
   o3 = orient2d_xy(q1u, q1v, q2u, q2v, p1u, p1v)
   o4 = orient2d_xy(q1u, q1v, q2u, q2v, p2u, p2v)
   hit = opposite(o1, o2) .and. opposite(o3, o4)
   if (hit) return
   hit = (o1 == 0._R8P .and. on_box(q1u, q1v, p1u, p1v, p2u, p2v)) .or. &
         (o2 == 0._R8P .and. on_box(q2u, q2v, p1u, p1v, p2u, p2v)) .or. &
         (o3 == 0._R8P .and. on_box(p1u, p1v, q1u, q1v, q2u, q2v)) .or. &
         (o4 == 0._R8P .and. on_box(p2u, p2v, q1u, q1v, q2u, q2v))
   contains
      pure function opposite(x, y) result(yes)
      !< x and y nonzero with opposite signs.
      real(R8P), intent(in) :: x, y !< Values.
      logical               :: yes  !< Result.

      yes = (x > 0._R8P .and. y < 0._R8P) .or. (x < 0._R8P .and. y > 0._R8P)
      endfunction opposite

      pure function on_box(pu, pv, au, av, bu, bv) result(yes)
      !< p in the bounding box of ab (p collinear with ab: then on ab).
      real(R8P), intent(in) :: pu, pv, au, av, bu, bv !< Points.
      logical               :: yes                    !< Result.

      yes = pu >= min(au, bu) .and. pu <= max(au, bu) .and. pv >= min(av, bv) .and. pv <= max(av, bv)
      endfunction on_box
   endfunction segments_intersect_2d

   pure function cross(x, y) result(z)
   !< Cross product.
   !$acc routine seq
   real(R8P), intent(in) :: x(3), y(3) !< Vectors.
   real(R8P)             :: z(3)       !< x cross y.

   z = [x(2) * y(3) - x(3) * y(2), x(3) * y(1) - x(1) * y(3), x(1) * y(2) - x(2) * y(1)]
   endfunction cross
endmodule fossil_box_predicates
