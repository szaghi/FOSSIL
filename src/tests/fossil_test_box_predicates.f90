!< FOSSIL test: box predicates (issue #38 step 1).
program fossil_test_box_predicates
!< FOSSIL test: box predicates (issue #38 step 1).
!<
!< 3-D box_triangle_overlap (separating axes, conservative):
!< - touching cases on the unit box: a vertex on a face, a triangle in a face plane, an edge along a
!<   box edge, a triangle through the interior with no vertex inside, a corner touched by the plane
!<   x + y + z = 3; separated cases: the same vertex moved off the face by 1e-7, the plane x + y + z
!<   = 3.2;
!< - 20000 random box/triangle pairs against an independent reference (Sutherland-Hodgman clipping of
!<   the triangle by the six box planes): no missed overlap; disagreements only on contacts.
!< 2-D box_segment_overlap_2d (exact) and box_segment_distance2_2d:
!< - touching: an end on an edge, a segment along an edge, through a corner; separated by one ulp;
!< - 20000 random pairs against Liang-Barsky clipping: identical answers away from contacts; the
!<   distance is 0 exactly on overlap and matches dense sampling otherwise.

use fossil_box_predicates, only : box_triangle_overlap, box_segment_overlap_2d, box_segment_distance2_2d
use penf,                  only : I4P, R8P

implicit none

real(R8P), parameter :: LO(3) = 0._R8P, HI(3) = 1._R8P
real(R8P)            :: r(9), blo(3), bhi(3), t1(3), t2(3), t3(3), b2lo(2), b2hi(2), a(2), b(2), d2, ds
integer(I4P)         :: k, miss, extra, ss, m, miss2, diff2, dbad
logical              :: ok(4), sat, ref

! 3-D touching / separated
ok(1) = box_triangle_overlap(LO, HI, [1._R8P, 0.5_R8P, 0.5_R8P], [2._R8P, 0.5_R8P, 0.5_R8P], [2._R8P, 0.6_R8P, 0.5_R8P], 0._R8P) &
  .and. .not. box_triangle_overlap(LO, HI, [1._R8P + 1.e-7_R8P, 0.5_R8P, 0.5_R8P], [2._R8P, 0.5_R8P, 0.5_R8P], &
                                   [2._R8P, 0.6_R8P, 0.5_R8P], 0._R8P) &
  .and. box_triangle_overlap(LO, HI, [0.2_R8P, 0.2_R8P, 1._R8P], [0.8_R8P, 0.2_R8P, 1._R8P], [0.5_R8P, 0.8_R8P, 1._R8P], 0._R8P) &
  .and. box_triangle_overlap(LO, HI, [1._R8P, 1._R8P, 0.2_R8P], [1._R8P, 1._R8P, 0.8_R8P], [2._R8P, 2._R8P, 0.5_R8P], 0._R8P) &
  .and. box_triangle_overlap(LO, HI, [-5._R8P, -5._R8P, 0.5_R8P], [5._R8P, -5._R8P, 0.5_R8P], [0._R8P, 5._R8P, 0.5_R8P], 0._R8P) &
  .and. box_triangle_overlap(LO, HI, [3._R8P, 0._R8P, 0._R8P], [0._R8P, 3._R8P, 0._R8P], [0._R8P, 0._R8P, 3._R8P], 0._R8P) &
  .and. .not. box_triangle_overlap(LO, HI, [3.2_R8P, 0._R8P, 0._R8P], [0._R8P, 3.2_R8P, 0._R8P], [0._R8P, 0._R8P, 3.2_R8P], 0._R8P)
print '(A,L1)', '3-D touching and separated cases: ', ok(1)

! 3-D random vs clipping
call random_seed(size=ss)
call random_seed(put=[(2024 + k, k = 1, ss)])
miss = 0 ; extra = 0
do k = 1, 20000
   call random_number(r)
   blo = r(1:3) * 0.7_R8P ; call random_number(r(1:3)) ; bhi = blo + 0.05_R8P + 0.25_R8P * r(1:3)
   call random_number(r)
   t1 = -0.2_R8P + 1.4_R8P * r(1:3) ; t2 = -0.2_R8P + 1.4_R8P * r(4:6) ; t3 = -0.2_R8P + 1.4_R8P * r(7:9)
   sat = box_triangle_overlap(blo, bhi, t1, t2, t3, 1.e-12_R8P)
   ref = clip_nonempty(blo, bhi, t1, t2, t3)
   if (ref .and. .not. sat) miss = miss + 1
   if (sat .and. .not. ref) extra = extra + 1
enddo
ok(2) = miss == 0 .and. extra == 0
print '(A,I0,A,I0)', '3-D random vs clipping: missed ', miss, '  extra ', extra

! 2-D touching / separated
b2lo = [0._R8P, 0._R8P] ; b2hi = [1._R8P, 1._R8P]
ok(3) = box_segment_overlap_2d(b2lo, b2hi, [2._R8P, 0.5_R8P], [1._R8P, 0.5_R8P]) .and. &
        box_segment_overlap_2d(b2lo, b2hi, [0.2_R8P, 1._R8P], [0.8_R8P, 1._R8P]) .and. &
        box_segment_overlap_2d(b2lo, b2hi, [2._R8P, 0._R8P], [0._R8P, 2._R8P]) .and. &
        .not. box_segment_overlap_2d(b2lo, b2hi, [2._R8P, 0.5_R8P], [1._R8P + spacing(1._R8P), 0.5_R8P]) .and. &
        .not. box_segment_overlap_2d(b2lo, b2hi, [2._R8P, 0._R8P], [0._R8P, 2._R8P + 2._R8P * spacing(2._R8P)]) .and. &
        box_segment_distance2_2d(b2lo, b2hi, [2._R8P, 0.5_R8P], [1._R8P, 0.5_R8P]) == 0._R8P .and. &
        abs(box_segment_distance2_2d(b2lo, b2hi, [3._R8P, 0.5_R8P], [2._R8P, 0.5_R8P]) - 1._R8P) < 1.e-15_R8P
print '(A,L1)', '2-D touching and separated cases: ', ok(3)

! 2-D random vs Liang-Barsky, distance vs sampling
miss2 = 0 ; diff2 = 0 ; dbad = 0
do k = 1, 20000
   call random_number(r(1:8))
   b2lo = r(1:2) * 0.7_R8P ; b2hi = b2lo + 0.05_R8P + 0.25_R8P * r(3:4)
   a = -0.2_R8P + 1.4_R8P * r(5:6) ; b = -0.2_R8P + 1.4_R8P * r(7:8)
   sat = box_segment_overlap_2d(b2lo, b2hi, a, b)
   ref = liang_barsky(b2lo, b2hi, a, b)
   if (sat .neqv. ref) diff2 = diff2 + 1
   d2 = box_segment_distance2_2d(b2lo, b2hi, a, b)
   if (sat .neqv. (d2 == 0._R8P)) dbad = dbad + 1
   ds = huge(1._R8P)
   do m = 0, 2000
      ds = min(ds, pbox(b2lo, b2hi, a + (b - a) * (m / 2000._R8P)))
   enddo
   if (d2 > ds + 1.e-14_R8P .or. sqrt(ds) - sqrt(d2) > 1.5_R8P * norm2(b - a) / 2000._R8P + 1.e-14_R8P) dbad = dbad + 1
enddo
ok(4) = diff2 == 0 .and. dbad == 0 .and. miss2 == 0
print '(A,I0,A,I0)', '2-D random vs Liang-Barsky: differences ', diff2, '  distance errors ', dbad

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function clip_nonempty(lo, hi, p1, p2, p3) result(nonempty)
   !< Sutherland-Hodgman: clip the triangle by the six planes of the box; true if anything remains.
   real(R8P), intent(in) :: lo(3), hi(3), p1(3), p2(3), p3(3) !< Box, triangle.
   logical               :: nonempty !< Result.
   real(R8P)             :: poly(3, 16), outp(3, 16), s, e, t !< Polygons.
   integer(I4P)          :: n, no, ax, sd, i, j !< Work.

   poly(:, 1) = p1 ; poly(:, 2) = p2 ; poly(:, 3) = p3 ; n = 3
   do ax = 1, 3
      do sd = 1, 2
         no = 0
         do i = 1, n
            j = mod(i, n) + 1
            s = side(poly(:, i), lo, hi, ax, sd) ; e = side(poly(:, j), lo, hi, ax, sd)
            if (s >= 0._R8P) then
               no = no + 1 ; outp(:, no) = poly(:, i)
            endif
            if ((s >= 0._R8P) .neqv. (e >= 0._R8P)) then
               t = s / (s - e)
               no = no + 1 ; outp(:, no) = poly(:, i) + t * (poly(:, j) - poly(:, i))
            endif
         enddo
         n = no ; poly(:, 1:n) = outp(:, 1:n)
         if (n == 0) then
            nonempty = .false.
            return
         endif
      enddo
   enddo
   nonempty = .true.
   endfunction clip_nonempty

   pure function side(q, lo, hi, ax, sd) result(v)
   !< Inside measure of q for the plane (axis ax, side sd: 1 low, 2 high) of the box (>= 0 inside).
   real(R8P),    intent(in) :: q(3), lo(3), hi(3) !< Point, box.
   integer(I4P), intent(in) :: ax, sd             !< Plane.
   real(R8P)                :: v                  !< Measure.

   if (sd == 1) then
      v = q(ax) - lo(ax)
   else
      v = hi(ax) - q(ax)
   endif
   endfunction side

   pure function liang_barsky(lo, hi, p, q) result(hit)
   !< Liang-Barsky clipping of the segment pq by the box: true if a part remains.
   real(R8P), intent(in) :: lo(2), hi(2), p(2), q(2) !< Box, segment.
   logical               :: hit !< Result.
   real(R8P)             :: t0, t1, dd(2), pp(4), qq(4), rr !< Work.
   integer(I4P)          :: i !< Counter.

   dd = q - p ; t0 = 0._R8P ; t1 = 1._R8P
   pp = [-dd(1), dd(1), -dd(2), dd(2)]
   qq = [p(1) - lo(1), hi(1) - p(1), p(2) - lo(2), hi(2) - p(2)]
   hit = .false.
   do i = 1, 4
      if (pp(i) == 0._R8P) then
         if (qq(i) < 0._R8P) return
      else
         rr = qq(i) / pp(i)
         if (pp(i) < 0._R8P) then
            t0 = max(t0, rr)
         else
            t1 = min(t1, rr)
         endif
      endif
   enddo
   hit = t0 <= t1
   endfunction liang_barsky

   pure function pbox(lo, hi, p) result(d)
   !< Squared distance from p to the box.
   real(R8P), intent(in) :: lo(2), hi(2), p(2) !< Box, point.
   real(R8P)             :: d    !< Squared distance.
   real(R8P)             :: e(2) !< Excess.

   e = max(lo - p, 0._R8P, p - hi)
   d = e(1)**2 + e(2)**2
   endfunction pbox
endprogram fossil_test_box_predicates
