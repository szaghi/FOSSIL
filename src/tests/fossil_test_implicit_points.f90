!< FOSSIL test: implicit points and exact predicates on them (issue #65, step 1).
!<
!< 1. Golden signs (src/tests/implicit_points_golden.txt, from an exact rational evaluation): orient3d of a line-plane
!<    intersection against a plane through it with one vertex moved by 0..50 ulps, and orient2d of two such points
!<    on nearly the same line with an endpoint of that line. Every sign must match; the same predicates evaluated on
!<    the points rounded to doubles must not (negative control: the rounding changes some signs).
!< 2. Exact zeros by construction, on random inputs: a line-plane point lies on its plane (orient3d = 0) and on its
!<    line (orient2d with the line's endpoints = 0 in every projection); a segment-segment point lies on both
!<    segments in its projection, and on its first (canonical) segment in 3D; reordering the inputs gives the same
!<    identity and the same rounded bits.
!< 3. Explicit points agree with fossil_predicates.
program fossil_test_implicit_points

use fossil_implicit_points, only : implicit_point_t, ip_explicit, ip_lpi, ip_ssi, ip_same, ip_approx, &
                                   ip_orient3d, ip_orient2d, ip_compare
use fossil_predicates,      only : orient3d_xyz, orient2d_xy
use penf,                   only : I4P, I8P, R8P

implicit none

type(implicit_point_t) :: x1, x2, e
character(2)           :: tag
real(R8P)              :: v(30), p(3), q(3), r(3), s(3), t(3), a(3), b(3), c(3), d(3)
integer(I4P)           :: u, ios, ex, nv, drop, ngold, bad_gold, naive_off, nzero_bad, nagree_bad, k, i
integer(I8P)           :: seed
logical                :: ok

ngold = 0 ; bad_gold = 0 ; naive_off = 0
open(newunit=u, file='src/tests/implicit_points_golden.txt', status='old', action='read')
read(u, *) ! header
do
   read(u, '(A2)', advance='no', iostat=ios) tag
   if (ios /= 0) exit
   if (tag == 'O3') then
      read(u, *) ex, nv, v(1:24)
      x1 = ip_lpi(v(1:3), v(4:6), v(7:9), v(10:12), v(13:15))
      k = ip_orient3d(x1, v(16:18), v(19:21), v(22:24))
   else
      read(u, *) ex, drop, nv, v(1:30)
      x1 = ip_lpi(v(1:3), v(4:6), v(7:9), v(10:12), v(13:15))
      x2 = ip_lpi(v(16:18), v(19:21), v(22:24), v(25:27), v(28:30))
      k = ip_orient2d(x1, x2, ip_explicit(v(1:3)), drop)
   endif
   ngold = ngold + 1
   if (k /= ex) bad_gold = bad_gold + 1
   if (nv /= ex) naive_off = naive_off + 1
enddo
close(u)
print '(A,I0,A,I0,A,I0)', 'golden cases ', ngold, ', exact disagreements ', bad_gold, &
                          ', disagreements of the rounded evaluation (negative control) ', naive_off

nzero_bad = 0 ; seed = 20261106_I8P
do i = 1, 500
   p = rv() ; q = rv() ; r = rv() ; s = rv() ; t = rv()
   x1 = ip_lpi(p, q, r, s, t)
   if (ip_orient3d(x1, r, s, t) /= 0) nzero_bad = nzero_bad + 1
   if (ip_orient3d(x1, t, r, s) /= 0) nzero_bad = nzero_bad + 1
   do drop = 1, 3
      if (ip_orient2d(x1, ip_explicit(p), ip_explicit(q), drop) /= 0) nzero_bad = nzero_bad + 1
   enddo
   x2 = ip_lpi(q, p, s, t, r) ! same construction, inputs reordered
   if (.not. ip_same(x1, x2)) nzero_bad = nzero_bad + 1
   if (any(ip_approx(x1) /= ip_approx(x2))) nzero_bad = nzero_bad + 1
   if (ip_compare(x1, x2, 1) /= 0) nzero_bad = nzero_bad + 1
   ! two points on the line (p, q) and p: collinear in every projection
   x2 = ip_lpi(p, q, rv(), rv(), rv())
   do drop = 1, 3
      if (ip_orient2d(x1, x2, ip_explicit(p), drop) /= 0) nzero_bad = nzero_bad + 1
   enddo
   ! segment-segment point: on both lines in its projection
   a = rv() ; b = rv() ; c = rv() ; d = rv() ; drop = int(mod(rnd(), 3_I8P), I4P) + 1
   x1 = ip_ssi(a, b, c, d, drop)
   if (ip_orient2d(x1, ip_explicit(a), ip_explicit(b), drop) /= 0) nzero_bad = nzero_bad + 1
   if (ip_orient2d(x1, ip_explicit(c), ip_explicit(d), drop) /= 0) nzero_bad = nzero_bad + 1
   do k = 1, 3 ! in 3D it lies on its first canonical segment, in every projection
      if (ip_orient2d(x1, ip_explicit(x1%v(:, 1)), ip_explicit(x1%v(:, 2)), k) /= 0) nzero_bad = nzero_bad + 1
   enddo
   if (.not. ip_same(x1, ip_ssi(d, c, b, a, drop))) nzero_bad = nzero_bad + 1
enddo
print '(A,I0)', 'exact-zero and identity checks failing: ', nzero_bad

nagree_bad = 0
do i = 1, 500
   a = rv() ; b = rv() ; c = rv() ; d = rv()
   e = ip_explicit(d)
   if (ip_orient3d(e, a, b, c) /= isgn(orient3d_xyz(a, b, c, d))) nagree_bad = nagree_bad + 1
   if (ip_orient2d(ip_explicit(a), ip_explicit(b), ip_explicit(c), 1) /= &
       isgn(orient2d_xy(a(2), a(3), b(2), b(3), c(2), c(3)))) nagree_bad = nagree_bad + 1
enddo
print '(A,I0)', 'explicit points disagreeing with fossil_predicates: ', nagree_bad

ok = ngold > 0 .and. bad_gold == 0 .and. naive_off > 0 .and. nzero_bad == 0 .and. nagree_bad == 0
print '(A,L1)', 'Are all tests passed? ', ok

contains
   function rnd() result(r)
   integer(I8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P) ; r = seed
   endfunction rnd

   function rv() result(x)
   real(R8P) :: x(3)
   integer(I4P) :: k

   do k = 1, 3
      x(k) = 2._R8P * real(rnd(), R8P) / 2147483647._R8P - 1._R8P
   enddo
   endfunction rv

   pure function isgn(x) result(s)
   real(R8P), intent(in) :: x
   integer(I4P)          :: s

   s = 0 ; if (x > 0._R8P) s = 1 ; if (x < 0._R8P) s = -1
   endfunction isgn
endprogram fossil_test_implicit_points
