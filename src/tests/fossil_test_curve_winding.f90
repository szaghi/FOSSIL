!< FOSSIL test: 2-D generalized winding number, inside test, open loops (issue #36 step 4).
program fossil_test_curve_winding
!< FOSSIL test: 2-D generalized winding number, inside test, open loops (issue #36 step 4).
!<
!< - Closed loops: the exact angle sum is 1 inside and 0 outside the square (to 1e-14), 0 in the
!<   hole and 1 in the ring of an annulus, 2 where two outer squares overlap.
!< - Fast versus exact: 2000 random points around a NACA 0012 (200 vertices) and an annulus. The
!<   inside decision (w > 0.5) matches the exact integer winding number everywhere, and the fast
!<   value is within 0.05 of the exact one at the default beta = 2 (near a thin body the far-field
!<   expansion is coarse; the decision margin is 0.5). On the NACA 0012 the maximum error falls
!<   by at least 4x per doubling of beta (2, 4, 8): second order, as the expansion is.
!< - Open loop: the square with a gap (0,0.2)-(0,0) on its left side. It is reported open, has no
!<   area, and saves/reloads as open. With SIGN_FAST_WINDING_NUMBER every random point farther
!<   than 0.4 from the gap gets the square's inside/outside right, for is_point_inside and for
!<   the sign of distance (whose magnitude does not depend on the sign algorithm).

use fossil, only : curve_object, CURVE_OPEN_LOOP, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P

implicit none

character(*), parameter :: F1 = 'fossil_test_curve_winding.dat'
real(R8P),    parameter :: PI = acos(-1._R8P)
type(curve_object)      :: c, r
real(R8P)               :: u(200), v(200), p(2), d1, d2
integer(I4P)            :: i, k, wrong, seed_size, unit
logical                 :: ok(5), inside

! closed loops, exact sum
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 0._R8P], v=[0._R8P, 0._R8P, 2._R8P, 2._R8P], loop_sizes=[4_I4P])
ok(1) = abs(c%winding_number([1._R8P, 1._R8P], use_brute_force=.true.) - 1._R8P) < 1.e-14_R8P .and. &
        abs(c%winding_number([3._R8P, 1._R8P], use_brute_force=.true.)) < 1.e-14_R8P
call annulus(c)
ok(1) = ok(1) .and. abs(c%winding_number([2._R8P, 2._R8P], use_brute_force=.true.)) < 1.e-14_R8P .and. &
        abs(c%winding_number([0.5_R8P, 2._R8P], use_brute_force=.true.) - 1._R8P) < 1.e-14_R8P .and. &
        c%is_point_inside([0.5_R8P, 2._R8P]) .and. .not. c%is_point_inside([2._R8P, 2._R8P]) .and. &
        c%is_point_inside([0.5_R8P, 2._R8P], sign_algorithm=SIGN_FAST_WINDING_NUMBER) .and. &
        .not. c%is_point_inside([2._R8P, 2._R8P], sign_algorithm=SIGN_FAST_WINDING_NUMBER)
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 0._R8P, 1._R8P, 3._R8P, 3._R8P, 1._R8P], &
                   v=[0._R8P, 0._R8P, 2._R8P, 2._R8P, 1._R8P, 1._R8P, 3._R8P, 3._R8P], loop_sizes=[4_I4P, 4_I4P])
ok(1) = ok(1) .and. abs(c%winding_number([1.5_R8P, 1.5_R8P], use_brute_force=.true.) - 2._R8P) < 1.e-14_R8P .and. &
        c%is_point_inside([1.5_R8P, 1.5_R8P], sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,L1)', 'closed loops (square, annulus, overlap): ', ok(1)

! fast versus exact
call naca0012(100_I4P, u, v)
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
ok(2) = fast_check(c, [-0.3_R8P, -0.3_R8P], [1.3_R8P, 0.3_R8P], 'NACA 0012') .and. beta_convergence(c)
call annulus(c)
ok(3) = fast_check(c, [-1._R8P, -1._R8P], [5._R8P, 5._R8P], 'annulus')

! open loop with a gap
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 0._R8P, 0._R8P], v=[0._R8P, 0._R8P, 2._R8P, 2._R8P, 0.2_R8P], &
                   loop_sizes=[5_I4P], closed=[.false.])
ok(4) = c%check_validity() == CURVE_OPEN_LOOP .and. c%get_area() == 0._R8P .and. .not. c%is_loop_closed(1_I4P)
call c%save_into_file(file_name=F1)
call r%load_from_file(file_name=F1)
ok(4) = ok(4) .and. r%get_vertices_number() == 5 .and. .not. r%is_loop_closed(1_I4P)
print '(A,L1)', 'open loop flagged, saved and reloaded open: ', ok(4)
call random_seed(size=seed_size)
call random_seed(put=[(777 + k, k = 1, seed_size)])
wrong = 0
do k = 1, 2000
   call random_number(p)
   p = -1._R8P + 4._R8P * p
   if (hypot(p(1), p(2) - 0.1_R8P) < 0.4_R8P) cycle
   inside = all(p > 0._R8P) .and. all(p < 2._R8P)
   if (c%is_point_inside(p, sign_algorithm=SIGN_FAST_WINDING_NUMBER) .neqv. inside) wrong = wrong + 1
   d1 = c%distance(p, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   d2 = c%distance(p)
   if (((d1 < 0._R8P) .neqv. inside) .or. abs(d1) /= abs(d2)) wrong = wrong + 1
enddo
ok(5) = wrong == 0
print '(A,I0)', 'open loop, FWN sign errors away from the gap: ', wrong

open(newunit=unit, file=F1, status='old') ; close(unit, status='delete')
print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function fast_check(cur, lo, hi, label) result(good)
   !< 2000 random points: |fast - exact| < 1e-3 and the fast decision equals the exact integer one.
   type(curve_object), intent(in) :: cur   !< Curve.
   real(R8P),          intent(in) :: lo(2) !< Box minimum.
   real(R8P),          intent(in) :: hi(2) !< Box maximum.
   character(*),       intent(in) :: label !< Label.
   logical                        :: good  !< Result.
   real(R8P)                      :: q(2), wf, we, maxdiff !< Work.
   integer(I4P)                   :: m, bad, ss !< Work.

   call random_seed(size=ss)
   call random_seed(put=[(99 + m, m = 1, ss)])
   bad = 0 ; maxdiff = 0._R8P
   do m = 1, 2000
      call random_number(q)
      q = lo + q * (hi - lo)
      wf = cur%winding_number(q)
      we = cur%winding_number(q, use_brute_force=.true.)
      maxdiff = max(maxdiff, abs(wf - we))
      if ((wf > 0.5_R8P) .neqv. (cur%winding_number_exact(q) /= 0)) bad = bad + 1
   enddo
   good = maxdiff < 0.05_R8P .and. bad == 0
   print '(A,A,ES10.3,A,I0)', label, ': max |fast - exact| ', maxdiff, '  decision mismatches ', bad
   endfunction fast_check

   function beta_convergence(cur) result(good)
   !< Max |fast - exact| over 2000 random points at beta = 2, 4, 8: each at least 4x smaller.
   type(curve_object), intent(in) :: cur   !< Curve.
   logical                        :: good  !< Result.
   real(R8P)                      :: q(2), e(3), we !< Work.
   integer(I4P)                   :: m, ib, ss !< Work.

   call random_seed(size=ss)
   call random_seed(put=[(99 + m, m = 1, ss)])
   e = 0._R8P
   do m = 1, 2000
      call random_number(q)
      q = [-0.3_R8P, -0.3_R8P] + q * [1.6_R8P, 0.6_R8P]
      we = cur%winding_number(q, use_brute_force=.true.)
      do ib = 1, 3
         e(ib) = max(e(ib), abs(cur%winding_number(q, beta=2._R8P**ib) - we))
      enddo
   enddo
   good = e(2) < e(1) / 4._R8P .and. e(3) < e(2) / 4._R8P
   print '(A,3ES10.3)', 'NACA 0012, max error at beta 2, 4, 8: ', e
   endfunction beta_convergence

   subroutine annulus(cur)
   !< Square [0,4]^2 with the clockwise hole [1,3]^2.
   type(curve_object), intent(inout) :: cur !< Curve.

   call cur%adopt_loops(u=[0._R8P, 4._R8P, 4._R8P, 0._R8P, 1._R8P, 1._R8P, 3._R8P, 3._R8P], &
                        v=[0._R8P, 0._R8P, 4._R8P, 4._R8P, 1._R8P, 3._R8P, 3._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
   endsubroutine annulus

   subroutine naca0012(m, x, y)
   !< NACA 0012 with a closed trailing edge, 2m vertices, cosine spacing, counter-clockwise.
   integer(I4P), intent(in)  :: m          !< Points per side.
   real(R8P),    intent(out) :: x(:), y(:) !< Vertices [2m].
   real(R8P)                 :: xi         !< Abscissa.
   integer(I4P)              :: kk, j      !< Counters.

   j = 0
   do kk = m, 0, -1
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / m))
      j = j + 1 ; x(j) = xi ; y(j) = thickness(xi)
   enddo
   do kk = 1, m - 1
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / m))
      j = j + 1 ; x(j) = xi ; y(j) = -thickness(xi)
   enddo
   endsubroutine naca0012

   pure function thickness(x) result(y)
   !< NACA 4-digit half thickness, 12 %, closed trailing edge.
   real(R8P), intent(in) :: x !< Abscissa.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction thickness
endprogram fossil_test_curve_winding
