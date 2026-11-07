!< FOSSIL test: accuracy of the closest-point kernels of the body-body distance (issue #69, step 10).
!<
!< 1. `closest_segment_segment` on 4000 fixed-seed pairs: near-parallel (angles 1e-2 .. 1e-15), offsets 1e-3 ..
!<    1e-12, overlapping, staggered and disjoint along their lines, plus degenerate (zero-length) segments. The
!<    distance |c1 - c2| against a quadruple-precision reference (real128, test side only: the four endpoint-segment
!<    distances, and the line-line distance when the closest points of the lines fall inside both segments, from the
!<    exactly computed cross product). Bound: |d - d_ref| <= 16 eps L, L the extent of the pair (measured: 0.62).
!<    Before issue #69 (Ericson 5.1.9) the error grew as eps L / angle: 139 of 4000 over the bound, by more than 1e5.
!< 2. `closest_on_triangle` on 4000 points against needle and regular triangles (aspect up to 1e8): |p - c| against
!<    the quadruple-precision distance (plane distance when the projection falls inside, else the three edges). Bound
!<    64 eps L: the point kernel `triangle_closest_st` (robust_st on needles, issue #39) errs ~eps/|n| on needles of
!<    aspect 1e8 (measured: 33). Before issue #69 (Ericson 5.1.5) 113 of 4000 were over it, by more than 1e5 eps L.
program fossil_test_segment_distance_accuracy

use, intrinsic :: iso_fortran_env, only : real128
use fossil_body_distance, only : closest_segment_segment, closest_on_triangle
use penf,                 only : I4P, I8P, R8P

implicit none

integer,   parameter :: QP = real128
real(R8P), parameter :: BOUND_SEG = 16._R8P, BOUND_TRI = 64._R8P
real(R8P)            :: p1(3), q1(3), p2(3), q2(3), c1(3), c2(3), t(3, 3), p(3), cc(3), d, dref, l, ratio1, ratio2
real(R8P)            :: ang, h, u(3), w(3), e(3)
integer(I4P)         :: k, nbad1, nbad2
integer(I8P)         :: seed

seed = 1069_I8P ; nbad1 = 0 ; ratio1 = 0._R8P
do k = 1, 4000
   u = unit3([rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P])
   w = unit3(cross3(u, [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P]))
   ang = 10._R8P**(-2._R8P - 13._R8P * rnd()) ; h = 10._R8P**(-3._R8P - 9._R8P * rnd())
   p1 = [rnd(), rnd(), rnd()] ; q1 = p1 + (0.2_R8P + rnd()) * u
   e = unit3(u + ang * w)
   p2 = p1 + h * unit3(cross3(u, w)) + (rnd() - 0.5_R8P) * 1.5_R8P * u
   q2 = p2 + (0.2_R8P + rnd()) * e
   if (mod(k, 50) == 0) q2 = p2 ! degenerate
   call closest_segment_segment(p1, q1, p2, q2, c1, c2)
   d = norm2(c1 - c2) ; dref = real(segseg_ref(p1, q1, p2, q2), R8P)
   l = max(maxval(abs([p1, q1, p2, q2])), 1._R8P)
   ratio1 = max(ratio1, abs(d - dref) / (epsilon(1._R8P) * l))
   if (abs(d - dref) > BOUND_SEG * epsilon(1._R8P) * l) nbad1 = nbad1 + 1
enddo
print '(A,I0,A,F10.2)', '1. segment-segment: over the bound ', nbad1, '  max |d - ref| / (eps L) ', min(ratio1, 99999._R8P)

nbad2 = 0 ; ratio2 = 0._R8P
do k = 1, 4000
   u = unit3([rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P])
   w = unit3(cross3(u, [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P]))
   t(:, 1) = [rnd(), rnd(), rnd()]
   t(:, 2) = t(:, 1) + (0.5_R8P + rnd()) * u
   t(:, 3) = t(:, 1) + (0.5_R8P * rnd()) * u + 10._R8P**(-8._R8P * rnd()) * w ! needles down to 1e-8
   p = 0.5_R8P * (t(:, 1) + t(:, 2)) + [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P] * 10._R8P**(-6._R8P * rnd())
   cc = closest_on_triangle(p, t)
   d = norm2(p - cc) ; dref = real(pointtri_ref(p, t), R8P)
   l = max(maxval(abs([p, t(:, 1), t(:, 2), t(:, 3)])), 1._R8P)
   ratio2 = max(ratio2, abs(d - dref) / (epsilon(1._R8P) * l))
   if (abs(d - dref) > BOUND_TRI * epsilon(1._R8P) * l) nbad2 = nbad2 + 1
enddo
print '(A,I0,A,F10.2)', '2. point-triangle: over the bound ', nbad2, '  max |d - ref| / (eps L) ', min(ratio2, 99999._R8P)
print '(A,L1)', 'Are all tests passed? ', nbad1 == 0 .and. nbad2 == 0
contains
   function segseg_ref(a1, b1, a2, b2) result(dd)
   !< Distance between segments a1-b1 and a2-b2 in quadruple precision.
   real(R8P), intent(in) :: a1(3), b1(3), a2(3), b2(3)
   real(QP)              :: dd, x1(3), y1(3), x2(3), y2(3), d1(3), d2(3), r(3), n(3), nn, s, tt

   x1 = real(a1, QP) ; y1 = real(b1, QP) ; x2 = real(a2, QP) ; y2 = real(b2, QP)
   dd = min(ps(x1, x2, y2), ps(y1, x2, y2), ps(x2, x1, y1), ps(y2, x1, y1))
   d1 = y1 - x1 ; d2 = y2 - x2 ; r = x2 - x1
   n = [d1(2) * d2(3) - d1(3) * d2(2), d1(3) * d2(1) - d1(1) * d2(3), d1(1) * d2(2) - d1(2) * d2(1)]
   nn = dot_product(n, n)
   if (nn > 0._QP) then
      s = dot_product([r(2) * d2(3) - r(3) * d2(2), r(3) * d2(1) - r(1) * d2(3), r(1) * d2(2) - r(2) * d2(1)], n) / nn
      tt = dot_product([r(2) * d1(3) - r(3) * d1(2), r(3) * d1(1) - r(1) * d1(3), r(1) * d1(2) - r(2) * d1(1)], n) / nn
      if (s > 0._QP .and. s < 1._QP .and. tt > 0._QP .and. tt < 1._QP) dd = min(dd, abs(dot_product(r, n)) / sqrt(nn))
   endif
   endfunction segseg_ref

   function ps(x, a, b) result(dd)
   !< Distance from x to the segment a-b, quadruple precision.
   real(QP), intent(in) :: x(3), a(3), b(3)
   real(QP)             :: dd, ab(3), s

   ab = b - a
   s = 0._QP
   if (dot_product(ab, ab) > 0._QP) s = min(1._QP, max(0._QP, dot_product(x - a, ab) / dot_product(ab, ab)))
   dd = sqrt(dot_product(x - a - s * ab, x - a - s * ab))
   endfunction ps

   function pointtri_ref(pp, tr) result(dd)
   !< Distance from pp to the triangle tr, quadruple precision.
   real(R8P), intent(in) :: pp(3), tr(3, 3)
   real(QP)              :: dd, x(3), a(3), b(3), cq(3), n(3), nn, o1, o2, o3, s

   x = real(pp, QP) ; a = real(tr(:, 1), QP) ; b = real(tr(:, 2), QP) ; cq = real(tr(:, 3), QP)
   dd = min(ps(x, a, b), ps(x, b, cq), ps(x, cq, a))
   n = cr(b - a, cq - a) ; nn = dot_product(n, n)
   if (nn > 0._QP) then
      o1 = dot_product(cr(b - a, x - a), n) ; o2 = dot_product(cr(cq - b, x - b), n) ; o3 = dot_product(cr(a - cq, x - cq), n)
      if (o1 >= 0._QP .and. o2 >= 0._QP .and. o3 >= 0._QP) then
         s = abs(dot_product(x - a, n)) / sqrt(nn)
         dd = min(dd, s)
      endif
   endif
   endfunction pointtri_ref

   pure function cr(x, y) result(z)
   real(QP), intent(in) :: x(3), y(3)
   real(QP)             :: z(3)

   z = [x(2) * y(3) - x(3) * y(2), x(3) * y(1) - x(1) * y(3), x(1) * y(2) - x(2) * y(1)]
   endfunction cr

   pure function cross3(x, y) result(z)
   real(R8P), intent(in) :: x(3), y(3)
   real(R8P)             :: z(3)

   z = [x(2) * y(3) - x(3) * y(2), x(3) * y(1) - x(1) * y(3), x(1) * y(2) - x(2) * y(1)]
   endfunction cross3

   pure function unit3(x) result(y)
   real(R8P), intent(in) :: x(3)
   real(R8P)             :: y(3)

   y = x / norm2(x)
   endfunction unit3

   function rnd() result(r)
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_segment_distance_accuracy
