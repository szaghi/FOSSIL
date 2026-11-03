!< FOSSIL test: batched curve distances, CPU and device code path (issue #36 step 7).
program fossil_test_curve_many
!< FOSSIL test: batched curve distances, CPU and device code path (issue #36 step 7).
!<
!< NACA 0012 (200 vertices), an annulus and an open square with a gap; 4000 random points each:
!< - distance_many equals compute_distance bit for bit: signed with the pseudo-normal, signed with
!<   the fast winding number, and unsigned;
!< - distance_many_device after enter_device equals distance_many bit for bit. Built without
!<   OpenACC (tests-gnu), the device routine runs its loop on the host through the same kernel:
!<   this checks the device code path's logic; the GPU check is the manual bench (pending, see #36);
!< - distance_many_device before enter_device, mismatched sizes and an unsupported sign algorithm
!<   are refused with STATUS_INVALID_INPUT.

use fossil, only : curve_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, SIGN_SOLID_ANGLE, STATUS_OK, &
                   STATUS_INVALID_INPUT
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter :: PI = acos(-1._R8P)
type(curve_object)   :: c
real(R8P)            :: u(200), v(200)
real(R8P)            :: dist(10)
integer(I4P)         :: st, st2, st3
logical              :: ok(4)

call naca0012(100_I4P, u, v)
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
ok(1) = check(c, [-0.3_R8P, -0.3_R8P], [1.3_R8P, 0.3_R8P], 'NACA 0012')
call c%adopt_loops(u=[0._R8P, 4._R8P, 4._R8P, 0._R8P, 1._R8P, 1._R8P, 3._R8P, 3._R8P], &
                   v=[0._R8P, 0._R8P, 4._R8P, 4._R8P, 1._R8P, 3._R8P, 3._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
ok(2) = check(c, [-1._R8P, -1._R8P], [5._R8P, 5._R8P], 'annulus')
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 0._R8P, 0._R8P], v=[0._R8P, 0._R8P, 2._R8P, 2._R8P, 0.2_R8P], &
                   loop_sizes=[5_I4P], closed=[.false.])
ok(3) = check(c, [-1._R8P, -1._R8P], [3._R8P, 3._R8P], 'open square')

! refusals
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
call c%distance_many_device(points=reshape([0.5_R8P, 0._R8P], [2, 1]), distances=dist(1:1), status=st)
call c%distance_many(points=reshape([0.5_R8P, 0._R8P], [2, 1]), distances=dist(1:2), status=st2)
call c%distance_many(points=reshape([0.5_R8P, 0._R8P], [2, 1]), distances=dist(1:1), sign_algorithm=SIGN_SOLID_ANGLE, &
                     status=st3)
ok(4) = st == STATUS_INVALID_INPUT .and. st2 == STATUS_INVALID_INPUT .and. st3 == STATUS_INVALID_INPUT
print '(A,L1)', 'refusals: ', ok(4)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check(cur, lo, hi, label) result(good)
   !< distance_many == compute_distance and distance_many_device == distance_many, bit for bit.
   type(curve_object), intent(inout) :: cur   !< Curve.
   real(R8P),          intent(in)    :: lo(2) !< Box minimum.
   real(R8P),          intent(in)    :: hi(2) !< Box maximum.
   character(*),       intent(in)    :: label !< Label.
   logical                           :: good  !< Result.
   integer(I4P), parameter           :: N = 4000 !< Points.
   real(R8P)                         :: p(2, N), dm(N), dd(N), dref !< Work.
   integer(I4P)                      :: k, a, mism, mdev, ss, status !< Work.
   integer(I4P), parameter           :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Sign algorithms.

   call random_seed(size=ss)
   call random_seed(put=[(31 + k, k = 1, ss)])
   call random_number(p)
   do k = 1, N
      p(:, k) = lo + p(:, k) * (hi - lo)
   enddo
   mism = 0 ; mdev = 0
   call cur%enter_device
   do a = 1, 3
      if (a <= 2) then
         call cur%distance_many(points=p, distances=dm, sign_algorithm=ALGO(a), status=status)
         call cur%distance_many_device(points=p, distances=dd, sign_algorithm=ALGO(a))
      else
         call cur%distance_many(points=p, distances=dm, is_signed=.false., status=status)
         call cur%distance_many_device(points=p, distances=dd, is_signed=.false.)
      endif
      do k = 1, N
         if (a <= 2) then
            call cur%compute_distance(p=p(:, k), distance=dref, sign_algorithm=ALGO(a))
         else
            call cur%compute_distance(p=p(:, k), distance=dref, is_signed=.false.)
         endif
         if (dm(k) /= dref) mism = mism + 1
         if (dd(k) /= dm(k)) mdev = mdev + 1
      enddo
   enddo
   call cur%exit_device
   good = mism == 0 .and. mdev == 0 .and. status == STATUS_OK
   print '(A,A,I0,A,I0)', label, ': distance_many vs compute_distance mismatches ', mism, '  device path vs CPU ', mdev
   endfunction check

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
endprogram fossil_test_curve_many
