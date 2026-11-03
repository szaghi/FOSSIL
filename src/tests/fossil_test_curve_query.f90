!< FOSSIL test: batched rich queries on curves (issue #37 step 4).
program fossil_test_curve_query
!< FOSSIL test: batched rich queries on curves (issue #37 step 4).
!<
!< - Unit square, analytic: below the bottom side (closest point (0.5, 0), segment 1, normal
!<   (0, -1)); beyond a corner (vertex 3, normal (1, 1)/sqrt2).
!< - A NACA 0012, an annulus with an island in its hole (bodies 1 and 2) and an open square;
!<   3000 random points each: every output of query_many equals compute_distance bit for bit, with
!<   the pseudo-normal and the fast-winding-number sign; bodies equal loop_body of the segment's loop;
!<   query_many_device (after enter_device; on the host in tests-gnu) equals query_many bit for bit.
!< - Size mismatch and device use before enter_device are refused.

use fossil, only : curve_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P)
type(curve_object)   :: c
real(R8P)            :: u(200), v(200), d(2), cp(2, 2), nn(2, 2)
integer(I4P)         :: sg(2), vx(2), st, st2, k
logical              :: ok(5)

call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call c%query_many(points=reshape([0.5_R8P, -1._R8P, 2._R8P, 2._R8P], [2, 2]), distances=d, closest_points=cp, &
                  segments=sg, vertices=vx, normals=nn)
ok(1) = d(1) == 1._R8P .and. all(cp(:, 1) == [0.5_R8P, 0._R8P]) .and. sg(1) == 1 .and. vx(1) == 0 .and. &
        all(nn(:, 1) == [0._R8P, -1._R8P]) .and. vx(2) == 3 .and. all(abs(nn(:, 2) - [1._R8P, 1._R8P] / sqrt(2._R8P)) < 1.e-15_R8P)
print '(A,L1)', 'square, segment and vertex analytic: ', ok(1)

do k = 0, 99
   u(k + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - k) / 100)) ; v(k + 1) = thickness(u(k + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do k = 1, 99
   u(101 + k) = 0.5_R8P * (1._R8P - cos(PI * k / 100)) ; v(101 + k) = -thickness(u(101 + k))
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
ok(2) = check(c, [-0.3_R8P, -0.3_R8P], [1.3_R8P, 0.3_R8P], 'NACA 0012')
call c%adopt_loops(u=[0._R8P, 6._R8P, 6._R8P, 0._R8P,  1._R8P, 1._R8P, 5._R8P, 5._R8P,  2._R8P, 4._R8P, 4._R8P, 2._R8P], &
                   v=[0._R8P, 0._R8P, 6._R8P, 6._R8P,  1._R8P, 5._R8P, 5._R8P, 1._R8P,  2._R8P, 2._R8P, 4._R8P, 4._R8P], &
                   loop_sizes=[4_I4P, 4_I4P, 4_I4P])
ok(3) = check(c, [-1._R8P, -1._R8P], [7._R8P, 7._R8P], 'annulus with an island')
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 0._R8P, 0._R8P], v=[0._R8P, 0._R8P, 2._R8P, 2._R8P, 0.2_R8P], &
                   loop_sizes=[5_I4P], closed=[.false.])
ok(4) = check(c, [-1._R8P, -1._R8P], [3._R8P, 3._R8P], 'open square')

call c%query_many(points=reshape([0.5_R8P, -1._R8P], [2, 1]), distances=d(1:1), closest_points=cp, status=st)
call c%query_many_device(points=reshape([0.5_R8P, -1._R8P], [2, 1]), distances=d(1:1), status=st2)
ok(5) = st == STATUS_INVALID_INPUT .and. st2 == STATUS_INVALID_INPUT
print '(A,L1)', 'refusals: ', ok(5)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check(cur, lo, hi, label) result(good)
   !< query_many vs compute_distance, and the device path vs query_many, every output, both signs.
   type(curve_object), intent(inout) :: cur   !< Curve.
   real(R8P),          intent(in)    :: lo(2) !< Box minimum.
   real(R8P),          intent(in)    :: hi(2) !< Box maximum.
   character(*),       intent(in)    :: label !< Label.
   logical                           :: good  !< Result.
   integer(I4P), parameter :: N = 3000 !< Points.
   integer(I4P), parameter :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.
   real(R8P)    :: p(2, N), d1(N), d2(N), c1(2, N), c2(2, N), n1(2, N), n2(2, N), ds, cs(2), ns(2) !< Work.
   integer(I4P) :: s1(N), s2(N), x1(N), x2(N), b1(N), b2(N), ss, xs, bs, m, a, mism, mdev, nseed, l, f, la !< Work.

   call random_seed(size=nseed)
   call random_seed(put=[(17 + m, m = 1, nseed)])
   call random_number(p)
   do m = 1, N
      p(:, m) = lo + p(:, m) * (hi - lo)
   enddo
   mism = 0 ; mdev = 0
   call cur%enter_device
   do a = 1, 2
      call cur%query_many(points=p, distances=d1, closest_points=c1, segments=s1, vertices=x1, normals=n1, bodies=b1, &
                          sign_algorithm=ALGO(a))
      call cur%query_many_device(points=p, distances=d2, closest_points=c2, segments=s2, vertices=x2, normals=n2, &
                                 bodies=b2, sign_algorithm=ALGO(a))
      do m = 1, N
         call cur%compute_distance(p=p(:, m), distance=ds, closest_point=cs, closest_segment=ss, closest_vertex=xs, &
                                   sign_algorithm=ALGO(a), normal=ns, body=bs)
         if (d1(m) /= ds .or. any(c1(:, m) /= cs) .or. s1(m) /= ss .or. x1(m) /= xs .or. any(n1(:, m) /= ns) .or. &
             b1(m) /= bs) mism = mism + 1
         do l = 1, cur%get_loops_number() ! body of the segment's loop
            call cur%get_loop_range(l, f, la)
            if (ss >= f .and. ss <= la .and. bs /= cur%loop_body(l)) mism = mism + 1
         enddo
         if (d2(m) /= d1(m) .or. any(c2(:, m) /= c1(:, m)) .or. s2(m) /= s1(m) .or. x2(m) /= x1(m) .or. &
             any(n2(:, m) /= n1(:, m)) .or. b2(m) /= b1(m)) mdev = mdev + 1
      enddo
   enddo
   call cur%exit_device
   good = mism == 0 .and. mdev == 0
   print '(A,A,I0,A,I0)', label, ': query_many vs compute_distance mismatches ', mism, '  device path vs CPU ', mdev
   endfunction check

   pure function thickness(x) result(y)
   !< NACA 4-digit half thickness, 12 %, closed trailing edge.
   real(R8P), intent(in) :: x !< Abscissa.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction thickness
endprogram fossil_test_curve_query
