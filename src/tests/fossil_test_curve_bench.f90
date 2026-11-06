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
!< Also `rich ... id mismatches` (query_many_device, issue #37) and `box device vs serial mismatches`
!< (box_classify_many_device on a box around every grid point, issue #38), and `sign bench device vs
!< serial mismatches` (sections of the issue #39 benchmark bodies), `deviation device vs serial mismatches` (issue #60:
!< the last curve rotated by 2 degrees against itself, on the device and on the CPU) and `line crossings device vs serial
!< mismatches` (issue #40) and `thin device vs serial mismatches` (issue #62: box_thin_many_device on the cells of a
!< 64^2 block, medial thickness and gap): all must be 0.

use fossil, only : curve_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, THICKNESS_MEDIAL
use penf,   only : I4P, I8P, R8P
use fossil_sign_bench_bodies, only : body_t, bench_bodies, device_sign_mismatches_2d, device_line_mismatches_2d
use fossil, only : rigid_pose_2d_t, pose_2d, to_world

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

call sign_bench
call deviation_case(c)
call feature_case(c)
call thin_case(c)
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
      call box(cur, p, (hi - lo) / real(max(n, 1), R8P), ALGO(a), NAME(a))
      call posed(cur, p, ALGO(a), NAME(a))
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
   ! patches (issue #66): a third of the segments in a second patch, the device's patches against the CPU's
   if (cur%get_patches_number() == 1_I4P) call cur%set_patch([(ids, ids = 1, cur%get_vertices_number(), 3)], 'third')
   call cur%query_many(points=p, distances=d1, segments=s1, patches=x1, sign_algorithm=algo)
   call cur%query_many_device(points=p, distances=d2, patches=x2, sign_algorithm=algo)
   ids = count(x1 /= x2)
   print '(A,A,A,I0)', '   ', nm, 'patch device vs serial mismatches: ', ids
   ok = ok .and. ids == 0
   endsubroutine rich

   subroutine sign_bench
   !< Device vs serial signed distances near the sections of the issue #39 benchmark bodies.
   type(body_t), allocatable :: bb(:) !< Bodies.
   integer                   :: nbad  !< Mismatches.

   call bench_bodies(bb)
   nbad = device_sign_mismatches_2d(bb)
   print '(A,I0)', 'sign bench device vs serial mismatches: ', nbad
   ok = ok .and. nbad == 0
   nbad = device_line_mismatches_2d(bb)
   print '(A,I0)', 'line crossings device vs serial mismatches: ', nbad
   ok = ok .and. nbad == 0
   endsubroutine sign_bench

   subroutine posed(cur, p, algo, nm)
   !< query_many_device vs query_many with an in-plane pose (issue #41), the grid points moved by it.
   type(curve_object), intent(inout) :: cur    !< Curve (on device).
   real(R8P),          intent(in)    :: p(:,:) !< Points [2, n].
   integer(I4P),       intent(in)    :: algo   !< Sign algorithm.
   character(*),       intent(in)    :: nm     !< Sign name.
   type(rigid_pose_2d_t)     :: pose !< Pose.
   real(R8P),    allocatable :: q(:,:), d1(:), d2(:), c1(:,:), c2(:,:) !< Work.
   integer(I4P)              :: nq, m, bad !< Work.

   pose = pose_2d(0.9_R8P, [5._R8P, -2._R8P])
   nq = size(p, 2)
   allocate(q(2, nq), d1(nq), d2(nq), c1(2, nq), c2(2, nq))
   do m = 1, nq
      q(:, m) = to_world(pose, p(:, m))
   enddo
   call cur%query_many(points=q, distances=d1, closest_points=c1, sign_algorithm=algo, pose=pose)
   call cur%query_many_device(points=q, distances=d2, closest_points=c2, sign_algorithm=algo, pose=pose)
   bad = count(d1 /= d2) + count(c1 /= c2)
   print '(A,A,A,I0)', '   ', nm, 'posed device vs serial mismatches: ', bad
   ok = ok .and. bad == 0
   endsubroutine posed

   subroutine deviation_case(cur)
   !< curve%deviation with device=.true. vs the CPU path (issue #60): the curve rotated by 2 degrees about the origin
   !< against the curve (entered), certified at 1e-4, lower, upper, status and per_vertex compared bit for bit.
   type(curve_object), intent(inout) :: cur        !< Curve B (entered here).
   type(curve_object)                :: moved      !< The rotated copy.
   real(R8P), allocatable            :: uu(:), vv(:), pv1(:), pv2(:) !< Coordinates, per-vertex distances.
   real(R8P)                         :: l1, u1, l2, u2, w(2), cs, sn !< Bounds, work.
   integer(I4P)                      :: s1, s2, i, nl, first, last, bad !< Statuses, counters.
   integer(I4P), allocatable         :: sizes(:)   !< Loop sizes.

   cs = cos(2._R8P * PI / 180._R8P) ; sn = sin(2._R8P * PI / 180._R8P)
   nl = cur%get_loops_number()
   allocate(uu(cur%get_vertices_number()), vv(cur%get_vertices_number()), sizes(nl))
   do i = 1, nl
      call cur%get_loop_range(i, first, last)
      sizes(i) = last - first + 1
   enddo
   do i = 1, cur%get_vertices_number()
      w = cur%get_vertex(i)
      uu(i) = cs * w(1) - sn * w(2) ; vv(i) = sn * w(1) + cs * w(2)
   enddo
   call moved%adopt_loops(u=uu, v=vv, loop_sizes=sizes)
   call moved%deviation(cur, lower=l1, upper=u1, per_vertex=pv1, tolerance=1.e-4_R8P, status=s1)
   call cur%enter_device
   call moved%deviation(cur, lower=l2, upper=u2, per_vertex=pv2, tolerance=1.e-4_R8P, device=.true., status=s2)
   call cur%exit_device
   bad = merge(0, 1, l1 == l2) + merge(0, 1, u1 == u2) + merge(0, 1, s1 == s2) + count(pv1 /= pv2)
   print '(A,2ES12.4,A,I0)', 'deviation of a 2-degree rotated copy: ', l1, u1, ', status ', s1
   print '(A,I0)', 'deviation device vs serial mismatches: ', bad
   ok = ok .and. bad == 0
   endsubroutine deviation_case

   subroutine feature_case(cur)
   !< Corner distances (issue #61): query_many_device vs query_many on a 64^2 grid over the curve, at 10 degrees.
   type(curve_object), intent(inout) :: cur  !< Curve (entered here).
   real(R8P), allocatable            :: p(:,:), d(:), f1(:), f2(:) !< Points, outputs.
   real(R8P)                         :: lo(2), hi(2) !< Box.
   integer(I4P)                      :: i, j, m, bad !< Counters.

   lo = cur%get_bmin() - 0.1_R8P ; hi = cur%get_bmax() + 0.1_R8P
   allocate(p(2, 64 * 64), d(64 * 64), f1(64 * 64), f2(64 * 64))
   m = 0
   do j = 1, 64 ; do i = 1, 64
      m = m + 1 ; p(:, m) = lo + (hi - lo) * [(i - 0.5_R8P) / 64, (j - 0.5_R8P) / 64]
   enddo ; enddo
   call cur%query_many(points=p, distances=d, feature_distances=f1, feature_angle=10._R8P)
   call cur%enter_device
   call cur%query_many_device(points=p, distances=d, feature_distances=f2, feature_angle=10._R8P)
   call cur%exit_device
   bad = count(f1 /= f2)
   print '(A,I0)', 'feature device vs serial mismatches: ', bad
   ok = ok .and. bad == 0
   endsubroutine feature_case

   subroutine thin_case(cur)
   !< Thin cells (issue #62): box_thin_many_device vs box_thin_many on the 64^2 cells of a block over the curve, the local
   !< sizes from the medial thickness and gap, factor 25 (the hole and the ring are 0.6 and 0.7 across, the cells 0.034).
   type(curve_object), intent(inout) :: cur  !< Curve (entered here).
   real(R8P), allocatable            :: lo(:,:), hi(:,:), th(:), gp(:) !< Cells, local sizes.
   logical,   allocatable            :: t1(:), t2(:) !< Touched.
   real(R8P)                         :: b0(2), b1(2) !< Box.
   integer(I4P)                      :: i, j, m, bad !< Counters.

   b0 = cur%get_bmin() - 0.1_R8P ; b1 = cur%get_bmax() + 0.1_R8P
   allocate(lo(2, 64 * 64), hi(2, 64 * 64), t1(64 * 64), t2(64 * 64))
   m = 0
   do j = 1, 64 ; do i = 1, 64
      m = m + 1 ; lo(:, m) = b0 + (b1 - b0) * [i - 1, j - 1] / 64._R8P ; hi(:, m) = b0 + (b1 - b0) * [i, j] / 64._R8P
   enddo ; enddo
   call cur%thickness(thickness=th, gap=gp, method=THICKNESS_MEDIAL)
   call cur%build_thin(thickness=th, gap=gp)
   call cur%box_thin_many(lo, hi, 25._R8P, t1)
   call cur%enter_device
   call cur%box_thin_many_device(lo, hi, 25._R8P, t2)
   call cur%exit_device
   bad = count(t1 .neqv. t2)
   print '(A,I0,A,I0)', 'thin cells touched: ', count(t1), ' of ', 64 * 64
   print '(A,I0)', 'thin device vs serial mismatches: ', bad
   ok = ok .and. bad == 0 .and. count(t1) > 0
   endsubroutine thin_case

   subroutine box(cur, p, h, algo, nm)
   !< box_classify_many_device vs box_classify_many (issue #38): a box of half-size h around every
   !< grid point, near radius |h|.
   type(curve_object), intent(inout) :: cur    !< Curve (on device).
   real(R8P),          intent(in)    :: p(:,:) !< Points [2, n].
   real(R8P),          intent(in)    :: h(2)   !< Half-sizes.
   integer(I4P),       intent(in)    :: algo   !< Sign algorithm.
   character(*),       intent(in)    :: nm     !< Sign name.
   real(R8P),    allocatable :: lo(:,:), hi(:,:) !< Boxes.
   integer(I4P), allocatable :: c1(:), c2(:) !< Codes.
   logical,      allocatable :: i1(:), i2(:) !< Inside flags.
   integer(I4P)              :: nq, m, bad !< Work.

   nq = size(p, 2)
   allocate(lo(2, nq), hi(2, nq), c1(nq), c2(nq), i1(nq), i2(nq))
   do m = 1, nq
      lo(:, m) = p(:, m) - h ; hi(:, m) = p(:, m) + h
   enddo
   call cur%box_classify_many(lo=lo, hi=hi, codes=c1, radius=norm2(h), sign_algorithm=algo, inside=i1)
   call cur%box_classify_many_device(lo=lo, hi=hi, codes=c2, radius=norm2(h), sign_algorithm=algo, inside=i2)
   bad = count(c1 /= c2 .or. (i1 .neqv. i2))
   print '(A,A,A,I0)', '   ', nm, 'box device vs serial mismatches: ', bad
   ok = ok .and. bad == 0
   endsubroutine box

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
