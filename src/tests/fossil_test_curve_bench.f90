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
use fossil_sign_bench_bodies, only : body_t, bench_bodies, device_sign_mismatches_2d, device_line_mismatches_2d, to_surface
use fossil, only : surface_stl_object, slice_surface, cut_cells_2d_t
use fossil_cut_cells_compare, only : cut_cells_twin_mismatches
use fossil, only : rigid_pose_2d_t, pose_2d, to_world
use fossil, only : curvilinear_block_2d_t, detect_topology_2d, topo_record_t, TOPO_BOUNDARY, wall_cutter_2d
use fossil_overset, only : hole_flags_2d, block_logical_t, donor_set_t, donor_candidates_2d

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
call block_case
call cut_cells_case
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

   subroutine block_case
   !< Point location in a curvilinear block on the device against the CPU (issue #71): a 2-D O-grid of 129 x 33
   !< vertices around the unit cylinder, 20000 points in its box; cells, xi (bit for bit) and statuses.
   type(curvilinear_block_2d_t) :: blk(1)                 !< Block.
   real(R8P),    allocatable    :: bu(:,:), bv(:,:), pu(:), pv(:), x1(:,:), x2(:,:) !< Vertices, points, xi.
   integer(I4P), allocatable    :: c1(:,:), c2(:,:), s1(:), s2(:) !< Results.
   integer(I4P)                 :: i, j, n, bad          !< Counters.
   real(R8P)                    :: r, rr(2)              !< Work.

   allocate(bu(129, 33), bv(129, 33))
   do j = 1, 33
      r = 1._R8P + 2._R8P * ((j - 1._R8P) / 32) ** 1.4_R8P
      do i = 1, 129
         bu(i, j) = r * cos(2 * PI * (i - 1) / 128) ; bv(i, j) = r * sin(2 * PI * (i - 1) / 128)
      enddo
   enddo
   bu(129,:) = bu(1,:) ; bv(129,:) = bv(1,:)
   call blk(1)%initialize(bu, bv) ; call detect_topology_2d(blk) ; call blk(1)%build_search
   allocate(pu(20000), pv(20000))
   do n = 1, 20000
      call random_number(rr) ; pu(n) = 6.2_R8P * rr(1) - 3.1_R8P ; pv(n) = 6.2_R8P * rr(2) - 3.1_R8P
   enddo
   call blk(1)%locate_many(pu, pv, c1, x1, s1)
   call blk(1)%locate_many_device(pu, pv, c2, x2, s2)
   bad = count(any(c1 /= c2, dim=1) .or. any(x1 /= x2, dim=1) .or. s1 /= s2)
   print '(A,I0,A,I0)', 'block location: points inside ', count(s1 == 1), ' of ', size(s1)
   print '(A,I0)', 'block location device vs serial mismatches: ', bad
   ok = ok .and. bad == 0 .and. count(s1 == 1) > 10000
   call hole_case(blk(1))
   endsubroutine block_case

   subroutine hole_case(og)
   !< Hole flags on the device against the CPU (issue #72): the wall of the O-grid cuts a 129^2 background block.
   type(curvilinear_block_2d_t), intent(in) :: og        !< O-grid.
   type(curvilinear_block_2d_t)             :: blk(2)    !< O-grid with its wall, background.
   type(curve_object)                       :: cut       !< Cutter.
   type(topo_record_t)                      :: rec       !< Wall record.
   type(block_logical_t), allocatable       :: h1(:), h2(:) !< Flags.
   type(donor_set_t),     allocatable       :: d1(:), d2(:) !< Donor candidates.
   real(R8P)                                :: bu(129, 129), bv(129, 129) !< Background.
   integer(I4P)                             :: i, j, bad !< Counters.

   blk(1) = og
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [og%ni, 1, 1] ; rec%name = 'wall'
   call blk(1)%add_record(rec)
   do j = 1, 129
      do i = 1, 129
         bu(i, j) = -3.1_R8P + 6.2_R8P * (i - 1) / 128 ; bv(i, j) = -3.1_R8P + 6.2_R8P * (j - 1) / 128
      enddo
   enddo
   call blk(2)%initialize(bu, bv)
   call wall_cutter_2d(blk(1:1), cut)
   call hole_flags_2d(blk, cut, h1)
   call hole_flags_2d(blk, cut, h2, device=.true.)
   bad = count(h1(2)%v .neqv. h2(2)%v) + count(h1(1)%v .neqv. h2(1)%v)
   print '(A,I0,A,I0)', 'hole flags: holes ', count(h1(2)%v), ' of ', size(h1(2)%v)
   print '(A,I0)', 'hole flags device vs serial mismatches: ', bad
   ok = ok .and. bad == 0 .and. count(h1(2)%v) > 1000
   call blk(1)%build_search ; call blk(2)%build_search
   call donor_candidates_2d(blk, d1)
   call blk(1)%enter_device ; call blk(2)%enter_device
   call donor_candidates_2d(blk, d2, device=.true.)
   call blk(1)%exit_device ; call blk(2)%exit_device
   bad = 0
   do i = 1, 2
      if (any(d1(i)%start /= d2(i)%start)) then
         bad = bad + 1 ; cycle
      endif
      bad = bad + count(any(d1(i)%cell /= d2(i)%cell, dim=1) .or. any(d1(i)%xi /= d2(i)%xi, dim=1) .or. &
                        d1(i)%block /= d2(i)%block .or. d1(i)%status /= d2(i)%status)
   enddo
   print '(A,I0)', 'donor candidates: ', size(d1(1)%block) + size(d1(2)%block)
   print '(A,I0)', 'donor candidates device vs serial mismatches: ', bad
   ok = ok .and. bad == 0 .and. size(d1(2)%block) > 1000
   endsubroutine hole_case

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

   subroutine cut_cells_case
   !< The 2-D cut cells on the device (issue #100 step 2): `cut_cells_device` against `cut_cells`, every output element
   !< bitwise, with the fluid edge intervals and without (`cut_cells_twin_mismatches`), on the mid-height sections of
   !< the benchmark bodies (n^2, n the bench size) and on four curves (n^2; with n >= 256 also (4n)^2, at most 1024^2, and
   !< a timing table: best of 5 calls, transfers included, intervals on). The regression run (n = 64) stays within its
   !< budget; the device validation runs `fossil_test_curve_bench 256`. `cut cells 2-D device vs serial mismatches` = 0.
   type(body_t), allocatable :: bb(:)                  !< Bodies.
   type(surface_stl_object)  :: s                      !< Surface.
   type(curve_object)        :: cc(4)                  !< Ellipse, NACA 0012, disc with a hole, 50 circles.
   type(cut_cells_2d_t)      :: cut_h, cut_d, cut_r    !< Cut cells (host, device, resident: each call refills its own).
   character(16)             :: cname(4)               !< Curve names.
   real(R8P),    allocatable :: uu(:), vv(:), x(:), y(:) !< Work.
   real(R8P)                 :: lo(2), hi(2), zmin, zmax, th, t(3, 5) !< Work.
   integer(I4P)              :: bad, st, sts, i, k, m, g, nn, r, ngrid(3), nl(50) !< Work.
   integer(I8P)              :: c0, c1, rate       !< Clock.

   ngrid = [n, min(4_I4P * n, 1024_I4P), 0_I4P]
   bad = 0_I4P ; sts = 0_I4P
   call bench_bodies(bb)
   do i = 1, size(bb)
      call to_surface(bb(i)%soup, s)
      zmin = minval(bb(i)%soup%v(3, :, :)) ; zmax = maxval(bb(i)%soup%v(3, :, :))
      call slice_surface(surface=s, normal_axis=3_I4P, level=0.5_R8P * (zmin + zmax), curve=cc(1))
      if (cc(1)%get_vertices_number() == 0) cycle
      lo = cc(1)%get_bmin() ; hi = cc(1)%get_bmax()
      call grid_edges(lo, hi, ngrid(1), x, y)
      bad = bad + cut_cells_twin_mismatches(cc(1), x, y, st) ; sts = sts + st
   enddo
   cname = [character(16) :: 'ellipse', 'NACA 0012', 'disc with a hole', '50 circles']
   allocate(uu(1024), vv(1024))
   do k = 1, 1024
      th = 2._R8P * PI * (k - 1) / 1024
      uu(k) = 0.5_R8P + 0.4_R8P * cos(th) ; vv(k) = 0.5_R8P + 0.25_R8P * sin(th)
   enddo
   call cc(1)%adopt_loops(u=uu, v=vv, loop_sizes=[1024_I4P])
   deallocate(uu, vv) ; allocate(uu(1022), vv(1022))
   call naca0012(511_I4P, uu, vv)
   uu = 0.15_R8P + 0.7_R8P * uu ; vv = 0.5_R8P + 0.7_R8P * vv
   call cc(2)%adopt_loops(u=uu, v=vv, loop_sizes=[1022_I4P])
   deallocate(uu, vv) ; allocate(uu(1100), vv(1100))
   do k = 1, 1000
      uu(k) = 0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (k - 1) / 1000) ; vv(k) = 0.5_R8P + 0.4_R8P * sin(2._R8P * PI * (k - 1) / 1000)
   enddo
   do k = 1, 100
      uu(1000 + k) = 0.5_R8P + 0.12_R8P * cos(-2._R8P * PI * (k - 1) / 100)
      vv(1000 + k) = 0.5_R8P + 0.12_R8P * sin(-2._R8P * PI * (k - 1) / 100)
   enddo
   call cc(3)%adopt_loops(u=uu, v=vv, loop_sizes=[1000_I4P, 100_I4P])
   deallocate(uu, vv) ; allocate(uu(3200), vv(3200))
   do m = 1, 50
      do k = 1, 64
         uu(64 * (m - 1) + k) = 0.1_R8P + 0.8_R8P * mod(m - 1, 10) / 9._R8P + 0.03_R8P * cos(2._R8P * PI * (k - 1) / 64)
         vv(64 * (m - 1) + k) = 0.1_R8P + 0.8_R8P * ((m - 1) / 10) / 4._R8P + 0.03_R8P * sin(2._R8P * PI * (k - 1) / 64)
      enddo
   enddo
   nl = 64_I4P
   call cc(4)%adopt_loops(u=uu, v=vv, loop_sizes=nl)
   lo = [0.0031_R8P, 0.0031_R8P] ; hi = [0.9969_R8P, 0.9969_R8P]
   do i = 1, 4
      do g = 1, merge(2, 1, n >= 256_I4P)
         call grid_edges(lo, hi, ngrid(g), x, y)
         bad = bad + cut_cells_twin_mismatches(cc(i), x, y, st) ; sts = sts + st
      enddo
   enddo
   print '(A,I0)', 'cut cells 2-D device vs serial mismatches: ', bad
   ok = ok .and. bad == 0 .and. sts == 0
   if (n < 256_I4P) return
   print '(A)', 'cut cells 2-D timing (5 calls, intervals on, transfers included) [s]: curve, cells, host, device, '// &
                'device resident (best, median); speed-ups of the best'
   ! warm-up (issue #100 step 6a): after the host-heavy twin checks the first device calls copied back at ~3.6 GB/s
   ! instead of ~13.5 GB/s for about a second (RTX 4070, WSL), whatever the code path; ten untimed calls end that window
   call grid_edges(lo, hi, ngrid(2), x, y)
   do r = 1, 10
      call cc(1)%cut_cells_device(x=x, y=y, cut=cut_d, status=st)
   enddo
   call system_clock(count_rate=rate)
   do i = 1, 4
      if (i == 3) cycle
      do g = 1, 2
         call grid_edges(lo, hi, ngrid(g), x, y)
         do r = 1, 5
            call system_clock(c0) ; call cc(i)%cut_cells(x=x, y=y, cut=cut_h, status=st) ; call system_clock(c1)
            t(1, r) = real(c1 - c0, R8P) / rate
            call system_clock(c0) ; call cc(i)%cut_cells_device(x=x, y=y, cut=cut_d, status=st) ; call system_clock(c1)
            t(2, r) = real(c1 - c0, R8P) / rate
            call system_clock(c0) ; call cc(i)%cut_cells_device(x=x, y=y, cut=cut_r, status=st, resident=.true.)
            call system_clock(c1)
            t(3, r) = real(c1 - c0, R8P) / rate
         enddo
         call cut_r%exit_device
         nn = ngrid(g)
         print '(3X,A16,1X,I4,A,3(2X,ES9.2,1X,ES9.2),2(A,F6.2),A)', cname(i), nn, '^2', &
               (minval(t(k, :)), median5(t(k, :)), k = 1, 3), '  (', minval(t(1, :)) / minval(t(2, :)), 'x, ', &
               minval(t(1, :)) / minval(t(3, :)), 'x)'
      enddo
   enddo
   endsubroutine cut_cells_case

   pure function median5(t) result(m)
   !< Median of five timings.
   real(R8P), intent(in) :: t(5) !< Timings.
   real(R8P)             :: m    !< Median.
   real(R8P)             :: w(5), z !< Work.
   integer(I4P)          :: a, b !< Counters.

   w = t
   do a = 2, 5
      z = w(a) ; b = a - 1
      do while (b >= 1)
         if (w(b) <= z) exit
         w(b + 1) = w(b) ; b = b - 1
      enddo
      w(b + 1) = z
   enddo
   m = w(3)
   endfunction median5

   subroutine grid_edges(lo, hi, n, x, y)
   !< n x n cells over the box [lo, hi] grown by 2 % (cell edges x(n+1), y(n+1)).
   real(R8P),              intent(in)  :: lo(2), hi(2) !< Box.
   integer(I4P),           intent(in)  :: n            !< Cells per axis.
   real(R8P), allocatable, intent(out) :: x(:), y(:)   !< Cell edges.
   real(R8P)                           :: a(2), b(2)   !< Grown box.
   integer(I4P)                        :: k            !< Counter.

   a = lo - 0.02_R8P * (hi - lo) ; b = hi + 0.02_R8P * (hi - lo)
   allocate(x(n + 1), y(n + 1))
   do k = 1, n + 1
      x(k) = a(1) + (b(1) - a(1)) * (k - 1) / real(n, R8P) ; y(k) = a(2) + (b(2) - a(2)) * (k - 1) / real(n, R8P)
   enddo
   endsubroutine grid_edges

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
