!< FOSSIL test/bench: 2-D batched distances, device against serial (issue #36 step 7).
program fossil_test_curve_bench
!< FOSSIL test/bench: 2-D batched distances, device against serial (issue #36 step 7).
!<
!< Usage: fossil_test_curve_bench [n]   (n x n query grid, default 64)
!<
!< A procedural NACA 0012 (2000 vertices) and a 1000-gon circle with a clockwise hole, queried on
!< an n x n grid around them: serial `compute_distance`, `distance_many` (CPU) and
!< `distance_many_device` (after `enter_device`), with the pseudo-normal and the fast-winding-
!< number signs. Prints the timings and `device vs serial max |d| err` (the manual GPU gate of
!< #20 / #36 requires 0.000E+00). Built by `tests-gnu` it runs the device routine's loop on the
!< host (no OpenACC); built by the `tests-nvidia-curve` mode it runs on the GPU.

use fossil, only : curve_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, I8P, R8P

implicit none

real(R8P), parameter      :: PI = acos(-1._R8P)
type(curve_object)        :: c
real(R8P), allocatable    :: u(:), v(:)
character(32)             :: arg
integer(I4P)              :: n, k
logical                   :: ok

n = 64
if (command_argument_count() >= 1) then
   call get_command_argument(1, arg)
   read(arg, *) n
endif
ok = .true.

allocate(u(2000), v(2000))
call naca0012(1000_I4P, u, v)
call c%adopt_loops(u=u, v=v, loop_sizes=[2000_I4P])
call bench(c, [-0.5_R8P, -0.5_R8P], [1.5_R8P, 0.5_R8P], 'NACA 0012 (2000 vertices)')

deallocate(u, v) ; allocate(u(1100), v(1100))
do k = 1, 1000 ! outer circle, counter-clockwise
   u(k) = cos(2._R8P * PI * (k - 1) / 1000) ; v(k) = sin(2._R8P * PI * (k - 1) / 1000)
enddo
do k = 1, 100  ! hole, clockwise
   u(1000 + k) = 0.3_R8P * cos(-2._R8P * PI * (k - 1) / 100) ; v(1000 + k) = 0.3_R8P * sin(-2._R8P * PI * (k - 1) / 100)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[1000_I4P, 100_I4P])
call bench(c, [-1.5_R8P, -1.5_R8P], [1.5_R8P, 1.5_R8P], 'disc with a hole (1100 vertices)')

print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine bench(cur, lo, hi, label)
   !< Serial, CPU-batched and device-path distances on an n x n grid; errors and timings.
   type(curve_object), intent(inout) :: cur   !< Curve.
   real(R8P),          intent(in)    :: lo(2) !< Grid minimum.
   real(R8P),          intent(in)    :: hi(2) !< Grid maximum.
   character(*),       intent(in)    :: label !< Label.
   real(R8P), allocatable            :: p(:,:), ds(:), dm(:), dd(:) !< Points, distances.
   real(R8P)                         :: err_cpu, err_dev, t_ser, t_cpu, t_dev !< Errors, timings.
   integer(I8P)                      :: t0, t1, rate !< Clock.
   integer(I4P)                      :: i, j, m, a   !< Counters.
   integer(I4P), parameter           :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.
   character(4), parameter           :: NAME(2) = ['pn  ', 'fwn '] !< Sign names.

   allocate(p(2, n * n), ds(n * n), dm(n * n), dd(n * n))
   m = 0
   do j = 1, n
      do i = 1, n
         m = m + 1
         p(:, m) = lo + (hi - lo) * [real(i - 1, R8P), real(j - 1, R8P)] / real(max(n - 1, 1), R8P)
      enddo
   enddo
   print '(A,A,I0,A,I0)', label, ', grid ', n, ' x ', n
   call cur%enter_device
   do a = 1, 2
      call system_clock(t0, rate)
      do m = 1, n * n
         call cur%compute_distance(p=p(:, m), distance=ds(m), sign_algorithm=ALGO(a))
      enddo
      call system_clock(t1) ; t_ser = real(t1 - t0, R8P) / rate
      call system_clock(t0)
      call cur%distance_many(points=p, distances=dm, sign_algorithm=ALGO(a))
      call system_clock(t1) ; t_cpu = real(t1 - t0, R8P) / rate
      call system_clock(t0)
      call cur%distance_many_device(points=p, distances=dd, sign_algorithm=ALGO(a))
      call system_clock(t1) ; t_dev = real(t1 - t0, R8P) / rate
      err_cpu = maxval(abs(dm - ds)) ; err_dev = maxval(abs(dd - ds))
      print '(A,A,A,3F10.4)', '   ', NAME(a), 'time serial / cpu batch / device [s]: ', t_ser, t_cpu, t_dev
      print '(A,A,A,ES10.3)', '   ', NAME(a), 'cpu batch vs serial max |d| err: ', err_cpu
      print '(A,A,A,ES10.3)', '   ', NAME(a), 'device vs serial max |d| err: ', err_dev
      ok = ok .and. err_cpu == 0._R8P .and. err_dev == 0._R8P
      call rich(cur, p, ALGO(a), NAME(a))
   enddo
   call cur%exit_device
   endsubroutine bench

   subroutine rich(cur, p, algo, nm)
   !< query_many_device vs query_many (issue #37): distance, closest point, normal, ids.
   type(curve_object), intent(inout) :: cur     !< Curve (on device).
   real(R8P),          intent(in)    :: p(:,:)  !< Points [2, n].
   integer(I4P),       intent(in)    :: algo    !< Sign algorithm.
   character(*),       intent(in)    :: nm      !< Sign name.
   real(R8P),    allocatable :: d1(:), d2(:), c1(:,:), c2(:,:), n1(:,:), n2(:,:) !< Outputs.
   integer(I4P), allocatable :: s1(:), s2(:), x1(:), x2(:), b1(:), b2(:) !< Outputs.
   real(R8P)                 :: err   !< Max error.
   integer(I4P)              :: nq, ids !< Work.

   nq = size(p, 2)
   allocate(d1(nq), d2(nq), c1(2, nq), c2(2, nq), n1(2, nq), n2(2, nq), s1(nq), s2(nq), x1(nq), x2(nq), b1(nq), b2(nq))
   call cur%query_many(points=p, distances=d1, closest_points=c1, segments=s1, vertices=x1, normals=n1, bodies=b1, &
                       sign_algorithm=algo)
   call cur%query_many_device(points=p, distances=d2, closest_points=c2, segments=s2, vertices=x2, normals=n2, bodies=b2, &
                              sign_algorithm=algo)
   err = max(maxval(abs(d1 - d2)), maxval(abs(c1 - c2)), maxval(abs(n1 - n2)))
   ids = count(s1 /= s2 .or. x1 /= x2 .or. b1 /= b2)
   print '(A,A,A,ES10.3,A,I0)', '   ', nm, 'rich device vs serial max err: ', err, '  id mismatches ', ids
   ok = ok .and. err == 0._R8P .and. ids == 0
   endsubroutine rich

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
endprogram fossil_test_curve_bench
