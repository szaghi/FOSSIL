!< FOSSIL benchmark: signed/unsigned distance-query throughput on a regular
!< grid of probe points — the Adam/Xall embedded-boundary access pattern
!< (one `surface%distance` call per Cartesian grid point against a static
!< surface). Establishes the baseline for the issue #19 distance hot-path
!< memory-layout overhaul; every optimization step in that issue is measured
!< against the numbers this harness prints.
!<
!< What it does:
!<   1. Loads the STL (default `src/tests/dragon-fine.stl`, ~24k facets).
!<   2. Builds a deterministic `N x N x N` probe grid spanning 1.5x the
!<      surface bounding box, so probes fall inside, outside, and near the
!<      surface. `N` is the second optional argument (default 40 -> 64000
!<      queries).
!<   3. Times the `SIGN_PSEUDO_NORMAL` signed query (the Xall path) and the
!<      unsigned query, reporting total wall-clock, queries/sec, and mean
!<      us/query. `cpu_time` is used for consistency with the other FOSSIL
!<      bench tests (`fossil_test_refactor_bench`, `fossil_test_pool_speedup`).
!<   4. Asserts a correctness invariant: on a small deterministic sub-sample
!<      the BVH-accelerated signed distance must be bit-exact against the
!<      brute-force scan (`aabb%set_use_index(.false.)`). This makes the
!<      harness double as a regression guard while the issue #19 layout
!<      changes land — any step that perturbs a distance value fails here.
!<
!<   5. Device lines (issues #20, #37, #38): `device`, `rich` and `box device vs serial mismatches`
!<      (box_classify_many_device on a box of half-size diagonal/64 around every probe) and
!<      `sign bench device vs serial mismatches` (the issue #39 benchmark points) and
!<      `line crossings device vs serial mismatches` (issue #40) must be 0.
!<
!< Usage: fossil_test_distance_bench [<stl_path> [<grid_n>]]
!<        default <stl_path> = src/tests/dragon-fine.stl
!<        default <grid_n>   = 40
!<
!< `perf stat` companion (baseline cache-miss / IPC capture for issue #19):
!<   perf stat -e cache-references,cache-misses,instructions,cycles \
!<       ./exe/fossil_test_distance_bench src/tests/dragon-fine.stl 40

program fossil_test_distance_bench

use fossil,              only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, &
                               SIGN_SOLID_ANGLE, STATUS_OK, STATUS_INVALID_INPUT
use fossil_facet_object, only : facet_object
use penf,                only : I4P, I8P, R8P
use vecfor,              only : vector_R8P
use fossil_distance_field, only : distance_field_object
use fossil_sign_bench_bodies, only : body_t, bench_bodies, device_sign_mismatches_3d, device_line_mismatches_3d
use fossil,              only : rigid_pose_t, pose_from_axis_angle, to_world

implicit none

type(surface_stl_object)      :: surface
type(vector_R8P)              :: bmin, bmax, span, probe
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_signed(:), d_unsigned(:), d_batch(:), d_device(:)
real(R8P),        allocatable :: d_fwn(:), d_fwn_device(:)   !< FWN-signed serial / device results (issue #30).
real(R8P)                     :: t0, t1, t_signed, t_unsigned
real(R8P)                     :: lo(3), hi(3), frac
real(R8P)                     :: d_tree, d_brute, max_abs_err, batch_abs_err
real(R8P)                     :: t_serial_wall, t_batch_wall, t_device_wall
real(R8P)                     :: device_abs_err
real(R8P)                     :: t_fwn_serial_wall, t_fwn_device_wall, fwn_device_abs_err
integer(I4P)                  :: fwn_sign_mismatch, fwn_inside, refuse_status
real(R8P)                     :: d_probe(1)                !< Output of the refused call (must stay untouched).
integer(I8P)                  :: clk0, clk1, clk_rate
integer(I4P)                  :: grid_n, n_queries, nf
integer(I4P)                  :: i, j, k, q
integer(I4P)                  :: argc, n_sample, s, stride
integer(I4P)                  :: n_threads
character(len=512)            :: stl_path, arg2
logical                       :: correctness_ok, batch_ok, device_ok, fwn_ok, refuse_ok, holed_ok, field_ok
logical                       :: rich_ok                   !< Rich queries: device equals CPU (issue #37).
real(R8P)                     :: rich_err                  !< Rich queries: max |device - CPU| over distance, closest point, normal.
integer(I4P)                  :: rich_id_mismatch          !< Rich queries: facet/feature/index/body mismatches.
integer(I4P)                  :: patch_mismatch            !< Rich queries: patch mismatches (issue #66).
integer(I4P)                  :: box_mismatch              !< Box classification: device vs CPU code/inside mismatches (issue #38).
integer                       :: sign_mismatch             !< Sign benchmark: device vs CPU signed-distance mismatches (issue #39).
integer                       :: line_mismatch             !< Grid-line crossings: device vs CPU mismatches (issue #40).
integer                       :: posed_mismatch            !< Posed queries: device vs CPU mismatches (issue #41).
integer                       :: deviation_mismatch        !< Surface deviation: device vs CPU mismatches (issue #60).
integer                       :: feature_mismatch          !< Feature-edge boxes: device vs CPU mismatches (issue #61).
integer                       :: thin_mismatch             !< Thin-facet boxes: device vs CPU mismatches (issue #62).
type(body_t), allocatable     :: sign_bodies(:)            !< Sign-benchmark bodies (issue #39).
real(R8P), parameter          :: BBOX_INFLATE = 1.5_R8P    !< Probe grid spans 1.5x the surface bbox.
integer(I4P), parameter       :: N_SAMPLE_MAX = 200_I4P    !< Sub-sample size for the brute-force correctness check.
!$ integer, external          :: omp_get_max_threads       !< OpenMP runtime thread-count query (only declared under -fopenmp).

! ---- arguments ---------------------------------------------------------------
argc = command_argument_count()
if (argc >= 1) then
   call get_command_argument(1, stl_path)
else
   stl_path = 'src/tests/dragon-fine.stl'
endif
grid_n = 40_I4P
if (argc >= 2) then
   call get_command_argument(2, arg2)
   read(arg2, *) grid_n
endif
if (grid_n < 2_I4P) grid_n = 2_I4P

! ---- load --------------------------------------------------------------------
call surface%load_from_file(file_name=trim(stl_path), guess_format=.true.)
nf = surface%get_facets_number()

! ---- build the deterministic probe grid -------------------------------------
bmin = surface%get_bmin()
bmax = surface%get_bmax()
span = bmax - bmin
! Inflate about the bbox centre so probes straddle the surface.
lo(1) = bmin%x - 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%x
lo(2) = bmin%y - 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%y
lo(3) = bmin%z - 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%z
hi(1) = bmax%x + 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%x
hi(2) = bmax%y + 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%y
hi(3) = bmax%z + 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%z

n_queries = grid_n * grid_n * grid_n
allocate(points(n_queries))
allocate(d_signed(n_queries))
allocate(d_unsigned(n_queries))
allocate(d_batch(n_queries))
allocate(d_device(n_queries))
allocate(d_fwn(n_queries))
allocate(d_fwn_device(n_queries))

q = 0_I4P
do k = 1_I4P, grid_n
   do j = 1_I4P, grid_n
      do i = 1_I4P, grid_n
         q = q + 1_I4P
         frac = real(i - 1_I4P, R8P) / real(grid_n - 1_I4P, R8P)
         probe%x = lo(1) + frac * (hi(1) - lo(1))
         frac = real(j - 1_I4P, R8P) / real(grid_n - 1_I4P, R8P)
         probe%y = lo(2) + frac * (hi(2) - lo(2))
         frac = real(k - 1_I4P, R8P) / real(grid_n - 1_I4P, R8P)
         probe%z = lo(3) + frac * (hi(3) - lo(3))
         points(q) = probe
      enddo
   enddo
enddo

! ---- timed: signed distance (SIGN_PSEUDO_NORMAL — the Xall path) -------------
call cpu_time(t0)
do q = 1_I4P, n_queries
   d_signed(q) = surface%distance(point=points(q), is_signed=.true., &
                                  sign_algorithm=SIGN_PSEUDO_NORMAL, is_square_root=.true.)
enddo
call cpu_time(t1)
t_signed = t1 - t0

! ---- timed: unsigned distance -----------------------------------------------
call cpu_time(t0)
do q = 1_I4P, n_queries
   d_unsigned(q) = surface%distance(point=points(q), is_square_root=.true.)
enddo
call cpu_time(t1)
t_unsigned = t1 - t0

! ---- timed: serial loop vs distance_many batch (issue #19 §B3) --------------
! Both timed with system_clock (wall time): cpu_time would sum across OpenMP
! threads and hide any parallel speedup. The serial loop here is the apples-to-
! apples baseline for the batch call below — same query, same options.
n_threads = 1_I4P
!$ n_threads = omp_get_max_threads()

call system_clock(count=clk0, count_rate=clk_rate)
do q = 1_I4P, n_queries
   d_signed(q) = surface%distance(point=points(q), is_signed=.true., &
                                  sign_algorithm=SIGN_PSEUDO_NORMAL, is_square_root=.true.)
enddo
call system_clock(count=clk1)
t_serial_wall = real(clk1 - clk0, R8P) / real(clk_rate, R8P)

call system_clock(count=clk0, count_rate=clk_rate)
call surface%distance_many(points=points, distances=d_batch, is_signed=.true., &
                           sign_algorithm=SIGN_PSEUDO_NORMAL, is_square_root=.true.)
call system_clock(count=clk1)
t_batch_wall = real(clk1 - clk0, R8P) / real(clk_rate, R8P)

! distance_many must be bit-exact against the serial loop — it is semantically a
! loop of the same `distance` calls, just OpenMP-parallel over disjoint outputs.
batch_abs_err = 0._R8P
do q = 1_I4P, n_queries
   batch_abs_err = max(batch_abs_err, abs(d_batch(q) - d_signed(q)))
enddo
batch_ok = (batch_abs_err == 0._R8P)

! ---- timed: distance_many_device (issue #20 §Step 5) -------------------------
! GPU-offloadable path: enter_device once, distance_many_device queries the BVH
! on-accelerator via `!$acc parallel loop`, exit_device releases. Under
! `tests-gnu` (no OpenACC compiler), the directives are inert comments and the
! loop runs sequentially in the host process — that gives us a CPU-side
! regression on the new code path. Under `tests-nvidia`, the same source runs
! on the GPU. Correctness gate: bit-exact vs the serial signed result, exactly
! like distance_many.
call surface%enter_device
call system_clock(count=clk0, count_rate=clk_rate)
call surface%distance_many_device(points=points, distances=d_device, is_signed=.true., &
                                  is_square_root=.true.)
call system_clock(count=clk1)
t_device_wall = real(clk1 - clk0, R8P) / real(clk_rate, R8P)

! ---- FWN sign on the device (issue #30) --------------------------------------
! Serial CPU reference first (one `distance` call per point, the FWN sign of
! `is_point_inside`), then the device run with SIGN_FAST_WINDING_NUMBER. Gate:
! bit-exact, like the pseudo-normal run above. The sign mismatch count is
! reported separately so a magnitude-only difference is told apart from a
! wrong inside/outside decision.
call system_clock(count=clk0, count_rate=clk_rate)
do q = 1_I4P, n_queries
   d_fwn(q) = surface%distance(point=points(q), is_signed=.true., &
                               sign_algorithm=SIGN_FAST_WINDING_NUMBER, is_square_root=.true.)
enddo
call system_clock(count=clk1)
t_fwn_serial_wall = real(clk1 - clk0, R8P) / real(clk_rate, R8P)

call system_clock(count=clk0, count_rate=clk_rate)
call surface%distance_many_device(points=points, distances=d_fwn_device, is_signed=.true., &
                                  sign_algorithm=SIGN_FAST_WINDING_NUMBER, is_square_root=.true.)
call system_clock(count=clk1)
t_fwn_device_wall = real(clk1 - clk0, R8P) / real(clk_rate, R8P)

! An unsupported sign algorithm is refused, not run: STATUS_INVALID_INPUT and
! the output left untouched (sentinel check).
d_probe(1) = -12345._R8P
call surface%distance_many_device(points=points(1:1), distances=d_probe, is_signed=.true., &
                                  sign_algorithm=SIGN_SOLID_ANGLE, status=refuse_status)
refuse_ok = (refuse_status == STATUS_INVALID_INPUT .and. d_probe(1) == -12345._R8P)
call rich_device_case(err=rich_err, id_mismatch=rich_id_mismatch, patch_mismatch=patch_mismatch)
rich_ok = rich_err == 0._R8P .and. rich_id_mismatch == 0_I4P .and. patch_mismatch == 0_I4P
call box_device_case(mismatch=box_mismatch)
rich_ok = rich_ok .and. box_mismatch == 0_I4P
call surface%exit_device
call bench_bodies(sign_bodies)
sign_mismatch = device_sign_mismatches_3d(sign_bodies)
rich_ok = rich_ok .and. sign_mismatch == 0
line_mismatch = device_line_mismatches_3d(sign_bodies)
rich_ok = rich_ok .and. line_mismatch == 0
call posed_case(mismatch=posed_mismatch)
rich_ok = rich_ok .and. posed_mismatch == 0
call deviation_case(mismatch=deviation_mismatch)
rich_ok = rich_ok .and. deviation_mismatch == 0
call feature_case(mismatch=feature_mismatch)
rich_ok = rich_ok .and. feature_mismatch == 0
call thin_case(mismatch=thin_mismatch)
rich_ok = rich_ok .and. thin_mismatch == 0

fwn_device_abs_err = 0._R8P
fwn_sign_mismatch  = 0_I4P
fwn_inside         = 0_I4P
do q = 1_I4P, n_queries
   fwn_device_abs_err = max(fwn_device_abs_err, abs(d_fwn_device(q) - d_fwn(q)))
   if ((d_fwn_device(q) < 0._R8P) .neqv. (d_fwn(q) < 0._R8P)) fwn_sign_mismatch = fwn_sign_mismatch + 1_I4P
   if (d_fwn(q) < 0._R8P) fwn_inside = fwn_inside + 1_I4P
enddo
fwn_ok = (fwn_device_abs_err == 0._R8P)

device_abs_err = 0._R8P
do q = 1_I4P, n_queries
   device_abs_err = max(device_abs_err, abs(d_device(q) - d_signed(q)))
enddo
device_ok = (device_abs_err == 0._R8P)

! ---- correctness invariant: BVH vs brute force on a sub-sample --------------
! Force the brute-force scan, recompute signed distance on a strided
! sub-sample, and require bit-exact agreement with the BVH result. Both tree
! kinds are documented to produce identical output; this guards that property
! through the issue #19 layout changes.
n_sample = min(N_SAMPLE_MAX, n_queries)
stride   = max(1_I4P, n_queries / n_sample)
max_abs_err = 0._R8P
call surface%aabb%set_use_index(.false.)   ! brute-force path
s = 0_I4P
do q = 1_I4P, n_queries, stride
   s = s + 1_I4P
   d_brute = surface%distance(point=points(q), is_signed=.true., &
                              sign_algorithm=SIGN_PSEUDO_NORMAL, is_square_root=.true.)
   d_tree  = d_signed(q)
   max_abs_err = max(max_abs_err, abs(d_tree - d_brute))
enddo
call surface%aabb%set_use_index(.true.)    ! restore
correctness_ok = (max_abs_err == 0._R8P)

! ---- holed mesh on the device (issue #30) -----------------------------------
call holed_device_case(ok=holed_ok)

! ---- distance-field band on the device (issue #32) ---------------------------
call field_band_case(ok=field_ok)

! ---- report ------------------------------------------------------------------
print '(A,A)',       '--- FOSSIL distance benchmark: ', trim(stl_path)
print '(A,I0)',      'facets:                       ', nf
print '(A,I0,A,I0)', 'probe grid:                   ', grid_n, ' ^3  ->  queries = ', n_queries
print '(A)',         ''
print '(A)',         'signed distance (SIGN_PSEUDO_NORMAL — the Xall path):'
print '(A,F12.4,A)', '   total wall-clock:        ', 1000._R8P * t_signed, ' ms'
print '(A,F12.2,A)', '   throughput:             ', real(n_queries, R8P) / max(t_signed, tiny(1._R8P)), ' queries/s'
print '(A,F12.4,A)', '   mean per query:         ', 1.0e6_R8P * t_signed / real(n_queries, R8P), ' us'
print '(A)',         ''
print '(A)',         'unsigned distance:'
print '(A,F12.4,A)', '   total wall-clock:        ', 1000._R8P * t_unsigned, ' ms'
print '(A,F12.2,A)', '   throughput:             ', real(n_queries, R8P) / max(t_unsigned, tiny(1._R8P)), ' queries/s'
print '(A,F12.4,A)', '   mean per query:         ', 1.0e6_R8P * t_unsigned / real(n_queries, R8P), ' us'
print '(A)',         ''
print '(A,I0,A)',    'distance_many batch (signed, OpenMP — ', n_threads, ' thread(s)):'
print '(A,F12.4,A)', '   serial loop wall-clock: ', 1000._R8P * t_serial_wall, ' ms'
print '(A,F12.4,A)', '   distance_many wall-clock:', 1000._R8P * t_batch_wall, ' ms'
print '(A,F12.4,A)', '   batch mean per query:   ', 1.0e6_R8P * t_batch_wall / real(n_queries, R8P), ' us'
print '(A,F12.2,A)', '   speedup vs serial:      ', t_serial_wall / max(t_batch_wall, tiny(1._R8P)), ' x'
print '(A,ES10.3)',  '   batch vs serial max |d| err: ', batch_abs_err
print '(A)',         ''
print '(A)',         'distance_many_device (signed, !$acc parallel loop -- inert under tests-gnu):'
print '(A,F12.4,A)', '   wall-clock:             ', 1000._R8P * t_device_wall, ' ms'
print '(A,F12.4,A)', '   mean per query:         ', 1.0e6_R8P * t_device_wall / real(n_queries, R8P), ' us'
print '(A,F12.2,A)', '   speedup vs serial:      ', t_serial_wall / max(t_device_wall, tiny(1._R8P)), ' x'
print '(A,ES10.3)',  '   device vs serial max |d| err: ', device_abs_err
print '(A)',         ''
print '(A)',         'distance_many_device (SIGN_FAST_WINDING_NUMBER, issue #30):'
print '(A,F12.4,A)', '   serial FWN wall-clock:  ', 1000._R8P * t_fwn_serial_wall, ' ms'
print '(A,F12.4,A)', '   device FWN wall-clock:  ', 1000._R8P * t_fwn_device_wall, ' ms'
print '(A,F12.4,A)', '   mean per query:         ', 1.0e6_R8P * t_fwn_device_wall / real(n_queries, R8P), ' us'
print '(A,F12.2,A)', '   speedup vs serial FWN:  ', t_fwn_serial_wall / max(t_fwn_device_wall, tiny(1._R8P)), ' x'
print '(A,I0,A,I0)', '   inside (serial FWN):    ', fwn_inside, ' of ', n_queries
print '(A,I0)',      '   sign mismatches:        ', fwn_sign_mismatch
print '(A,ES10.3)',  '   fwn device vs serial max |d| err: ', fwn_device_abs_err
print '(A,L1)',      '   unsupported sign refused:  ', refuse_ok
print '(A)',         ''
print '(A)',         'query_many_device vs query_many (issue #37: distance, closest point, normal, ids; both signs):'
print '(A,ES10.3)',  '   rich device vs serial max err: ', rich_err
print '(A,I0)',      '   rich id mismatches:      ', rich_id_mismatch
print '(A,I0)',      '   patch device vs serial mismatches: ', patch_mismatch
print '(A,I0)',      '   box device vs serial mismatches: ', box_mismatch
print '(A,I0)',      '   sign bench device vs serial mismatches: ', sign_mismatch
print '(A,I0)',      '   line crossings device vs serial mismatches: ', line_mismatch
print '(A,I0)',      '   posed device vs serial mismatches: ', posed_mismatch
print '(A,I0)',      '   deviation device vs serial mismatches: ', deviation_mismatch
print '(A,I0)',      '   feature device vs serial mismatches: ', feature_mismatch
print '(A,I0)',      '   thin device vs serial mismatches: ', thin_mismatch
print '(A)',         ''
print '(A,I0,A,ES10.3)', 'correctness (BVH vs brute force, ', s, ' samples): max |d_tree - d_brute| = ', max_abs_err
print '(A)',         ''
print '(A,L1)',      'Are all tests passed? ', (correctness_ok .and. batch_ok .and. device_ok .and. fwn_ok .and. refuse_ok &
                                                .and. holed_ok .and. field_ok .and. rich_ok)
if (.not. (correctness_ok .and. batch_ok .and. device_ok .and. fwn_ok .and. refuse_ok .and. holed_ok .and. field_ok &
           .and. rich_ok)) &
   error stop 1

contains
   subroutine rich_device_case(err, id_mismatch, patch_mismatch)
   !< query_many_device against query_many on the bench points (surface entered), both signs.
   real(R8P),    intent(out) :: err         !< Max |device - CPU| over distance, closest point, normal.
   integer(I4P), intent(out) :: id_mismatch !< Facet, feature, index, body mismatches.
   integer(I4P), intent(out) :: patch_mismatch !< Patch mismatches (issue #66; a third of the facets in a 2nd patch).
   integer(I4P), allocatable :: p1(:), p2(:) !< Patches.
   integer(I4P)              :: k !< Counter.
   real(R8P),    allocatable :: d1(:), d2(:), c1(:,:), c2(:,:), n1(:,:), n2(:,:) !< Outputs.
   integer(I4P), allocatable :: f1(:), f2(:), t1(:), t2(:), i1(:), i2(:), b1(:), b2(:) !< Outputs.
   integer(I4P)              :: a, nq !< Work.
   integer(I4P), parameter   :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.

   nq = size(points, kind=I4P)
   allocate(d1(nq), d2(nq), c1(3, nq), c2(3, nq), n1(3, nq), n2(3, nq))
   allocate(f1(nq), f2(nq), t1(nq), t2(nq), i1(nq), i2(nq), b1(nq), b2(nq))
   err = 0._R8P ; id_mismatch = 0_I4P ; patch_mismatch = 0_I4P
   allocate(p1(nq), p2(nq))
   call surface%set_patch([(k, k = 1_I4P, surface%get_facets_number(), 3_I4P)], 'third')
   do a = 1, 2
      call surface%query_many(points=points, distances=d1, closest_points=c1, facets=f1, features=t1, feature_indexes=i1, &
                              normals=n1, bodies=b1, sign_algorithm=ALGO(a))
      call surface%query_many_device(points=points, distances=d2, closest_points=c2, facets=f2, features=t2, &
                                     feature_indexes=i2, normals=n2, bodies=b2, sign_algorithm=ALGO(a))
      err = max(err, maxval(abs(d1 - d2)), maxval(abs(c1 - c2)), maxval(abs(n1 - n2)))
      id_mismatch = id_mismatch + count(f1 /= f2 .or. t1 /= t2 .or. i1 /= i2 .or. b1 /= b2, kind=I4P)
      call surface%query_many(points=points, distances=d1, patches=p1, sign_algorithm=ALGO(a))
      call surface%query_many_device(points=points, distances=d2, patches=p2, sign_algorithm=ALGO(a))
      patch_mismatch = patch_mismatch + count(p1 /= p2, kind=I4P)
   enddo
   endsubroutine rich_device_case

   subroutine posed_case(mismatch)
   !< query_many_device vs query_many with a pose (issue #41) on the bench points moved by that pose.
   integer, intent(out)          :: mismatch !< Mismatches over distance, closest point, normal, facet.
   type(rigid_pose_t)            :: pose     !< Pose.
   type(vector_R8P), allocatable :: q(:)     !< Moved points.
   real(R8P),        allocatable :: d1(:), d2(:), c1(:,:), c2(:,:), n1(:,:), n2(:,:) !< Outputs.
   integer(I4P),     allocatable :: f1(:), f2(:) !< Outputs.
   real(R8P)                     :: w(3)     !< Work.
   integer(I4P)                  :: m, nq    !< Work.

   pose = pose_from_axis_angle([1._R8P, -2._R8P, 0.5_R8P], 1.1_R8P, [20._R8P, 3._R8P, -7._R8P])
   nq = size(points, kind=I4P)
   allocate(q(nq), d1(nq), d2(nq), c1(3, nq), c2(3, nq), n1(3, nq), n2(3, nq), f1(nq), f2(nq))
   do m = 1, nq
      w = to_world(pose, [points(m)%x, points(m)%y, points(m)%z]) ; q(m) = vector_R8P(w(1), w(2), w(3))
   enddo
   call surface%enter_device
   call surface%query_many(points=q, distances=d1, closest_points=c1, normals=n1, facets=f1, pose=pose)
   call surface%query_many_device(points=q, distances=d2, closest_points=c2, normals=n2, facets=f2, pose=pose)
   call surface%exit_device
   mismatch = count(d1 /= d2) + count(c1 /= c2) + count(n1 /= n2) + count(f1 /= f2)
   endsubroutine posed_case

   subroutine deviation_case(mismatch)
   !< surface%deviation with device=.true. vs the CPU path (issue #60): a copy of the surface rotated by 2 degrees about
   !< its centre against the surface (entered), certified bounds at 1e-2 of the diagonal, work capped at 200 000
   !< triangles; lower, upper, status and per_vertex compared bit for bit.
   integer, intent(out)     :: mismatch     !< Mismatches.
   type(surface_stl_object) :: moved        !< The rotated copy.
   type(rigid_pose_t)       :: pose         !< Rotation about the centre.
   type(vector_R8P)         :: lo, hi       !< Bounding box.
   real(R8P), allocatable   :: pv1(:), pv2(:) !< Per-vertex distances.
   real(R8P)                :: l1, u1, l2, u2, c(3), tol !< Bounds, centre, tolerance.
   integer(I4P)             :: s1, s2       !< Statuses.

   lo = surface%get_bmin() ; hi = surface%get_bmax()
   c = 0.5_R8P * [lo%x + hi%x, lo%y + hi%y, lo%z + hi%z]
   tol = 1.e-2_R8P * norm2([hi%x - lo%x, hi%y - lo%y, hi%z - lo%z])
   pose = pose_from_axis_angle([1._R8P, 1._R8P, 1._R8P], 2._R8P * acos(-1._R8P) / 180._R8P)
   pose%t = c - matmul(pose%r, c)
   moved = surface
   call moved%move(pose)
   call moved%deviation(surface, lower=l1, upper=u1, per_vertex=pv1, tolerance=tol, max_triangles=200000_I4P, status=s1)
   call surface%enter_device
   call moved%deviation(surface, lower=l2, upper=u2, per_vertex=pv2, tolerance=tol, max_triangles=200000_I4P, &
                        device=.true., status=s2)
   call surface%exit_device
   mismatch = merge(0, 1, l1 == l2) + merge(0, 1, u1 == u2) + merge(0, 1, s1 == s2) + count(pv1 /= pv2)
   print '(A,2ES12.4,A,I0)', '   deviation of a 2-degree rotated copy: ', l1, u1, ', status ', s1
   endsubroutine deviation_case

   subroutine feature_case(mismatch)
   !< box_feature_many_device vs box_feature_many and query_many_device vs query_many feature distances (issue #61): the
   !< bench surface's feature edges (default 30 degrees), a 24^3 block of cells over its bounding box and their centres.
   integer, intent(out)     :: mismatch     !< Mismatches.
   integer(I4P), parameter  :: NC = 24      !< Cells per direction.
   type(vector_R8P)         :: lo_, hi_     !< Bounding box.
   real(R8P), allocatable   :: blo(:,:), bhi(:,:) !< Cell boxes.
   logical,   allocatable   :: t1(:), t2(:) !< Touched.
   type(vector_R8P), allocatable :: cq(:)   !< Cell centres.
   real(R8P), allocatable   :: fd1(:), fd2(:), dq(:) !< Feature distances, distances.
   real(R8P)                :: h(3)         !< Cell size.
   integer(I4P)             :: i, j, k, c, st !< Counters.

   lo_ = surface%get_bmin() ; hi_ = surface%get_bmax()
   h = [hi_%x - lo_%x, hi_%y - lo_%y, hi_%z - lo_%z] / NC
   allocate(blo(3, NC**3), bhi(3, NC**3), t1(NC**3), t2(NC**3))
   c = 0
   do k = 1, NC ; do j = 1, NC ; do i = 1, NC
      c = c + 1 ; blo(:, c) = [lo_%x, lo_%y, lo_%z] + h * [i - 1, j - 1, k - 1] ; bhi(:, c) = blo(:, c) + h
   enddo ; enddo ; enddo
   call surface%build_features(status=st)
   call surface%box_feature_many(blo, bhi, t1)
   allocate(cq(NC**3), fd1(NC**3), fd2(NC**3), dq(NC**3))
   do c = 1, NC**3                                   ! feature distances at the cell centres (step 3)
      cq(c) = vector_R8P(0.5_R8P * (blo(1, c) + bhi(1, c)), 0.5_R8P * (blo(2, c) + bhi(2, c)), 0.5_R8P * (blo(3, c) + bhi(3, c)))
   enddo
   call surface%query_many(points=cq, distances=dq, feature_distances=fd1)
   call surface%enter_device
   call surface%box_feature_many_device(blo, bhi, t2)
   call surface%query_many_device(points=cq, distances=dq, feature_distances=fd2)
   call surface%exit_device
   mismatch = count(t1 .neqv. t2) + count(fd1 /= fd2)
   print '(A,I0,A,I0)', '   feature cells touched: ', count(t1), ' of ', NC**3
   endsubroutine feature_case

   subroutine thin_case(mismatch)
   !< box_thin_many_device vs box_thin_many (issue #62): the axial-ray thickness of the bench surface (rays=1; open
   !< surfaces leave huge values, never thin), build_thin, a 24^3 block of cells over its bounding box, the factor set so
   !< the threshold is the mean finite thickness on the largest cell side.
   integer, intent(out)     :: mismatch     !< Mismatches.
   integer(I4P), parameter  :: NC = 24      !< Cells per direction.
   type(vector_R8P)         :: lo_, hi_     !< Bounding box.
   real(R8P), allocatable   :: blo(:,:), bhi(:,:), th(:) !< Cell boxes, thickness.
   logical,   allocatable   :: t1(:), t2(:) !< Touched.
   real(R8P)                :: h(3), factor !< Cell size, threshold factor.
   integer(I4P)             :: i, j, k, c, st !< Counters.

   lo_ = surface%get_bmin() ; hi_ = surface%get_bmax()
   h = [hi_%x - lo_%x, hi_%y - lo_%y, hi_%z - lo_%z] / NC
   allocate(blo(3, NC**3), bhi(3, NC**3), t1(NC**3), t2(NC**3))
   c = 0
   do k = 1, NC ; do j = 1, NC ; do i = 1, NC
      c = c + 1 ; blo(:, c) = [lo_%x, lo_%y, lo_%z] + h * [i - 1, j - 1, k - 1] ; bhi(:, c) = blo(:, c) + h
   enddo ; enddo ; enddo
   call surface%thickness(thickness=th, rays=1_I4P, status=st)
   factor = sum(th, th < huge(1._R8P)) / max(1, count(th < huge(1._R8P))) / maxval(h)
   call surface%build_thin(thickness=th, status=st)
   call surface%box_thin_many(blo, bhi, factor, t1)
   call surface%enter_device
   call surface%box_thin_many_device(blo, bhi, factor, t2)
   call surface%exit_device
   mismatch = count(t1 .neqv. t2)
   print '(A,I0,A,I0)', '   thin cells touched: ', count(t1), ' of ', NC**3
   endsubroutine thin_case

   subroutine box_device_case(mismatch)
   !< box_classify_many_device against box_classify_many (issue #38): a box of half-size 1/64 of the
   !< surface diagonal around every bench point, near radius of the same size, both signs.
   integer(I4P), intent(out) :: mismatch !< Code or inside-flag mismatches.
   real(R8P),    allocatable :: lo(:,:), hi(:,:) !< Boxes.
   integer(I4P), allocatable :: c1(:), c2(:) !< Codes.
   logical,      allocatable :: i1(:), i2(:) !< Inside flags.
   type(vector_R8P)          :: bmin, bmax   !< Surface box.
   real(R8P)                 :: h            !< Half-size.
   integer(I4P)              :: a, nq, m     !< Work.
   integer(I4P), parameter   :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.

   nq = size(points, kind=I4P)
   allocate(lo(3, nq), hi(3, nq), c1(nq), c2(nq), i1(nq), i2(nq))
   bmin = surface%get_bmin() ; bmax = surface%get_bmax()
   h = sqrt((bmax%x - bmin%x)**2 + (bmax%y - bmin%y)**2 + (bmax%z - bmin%z)**2) / 64._R8P
   do m = 1, nq
      lo(:, m) = [points(m)%x, points(m)%y, points(m)%z] - h
      hi(:, m) = [points(m)%x, points(m)%y, points(m)%z] + h
   enddo
   mismatch = 0_I4P
   do a = 1, 2
      call surface%box_classify_many(lo=lo, hi=hi, codes=c1, radius=h, sign_algorithm=ALGO(a), inside=i1)
      call surface%box_classify_many_device(lo=lo, hi=hi, codes=c2, radius=h, sign_algorithm=ALGO(a), inside=i2)
      mismatch = mismatch + count(c1 /= c2 .or. (i1 .neqv. i2), kind=I4P)
   enddo
   endsubroutine box_device_case
   subroutine field_band_case(ok)
   !< Narrow band and signs of a distance field over the probe box (grid_n^3 points, PN sign):
   !< band and sign point queries through distance_many_device must equal the host ones bit for
   !< bit (issue #32). Timings printed (host: distance_many, OpenMP if enabled).
   logical, intent(out)        :: ok              !< Pass flag.
   type(distance_field_object) :: host, dev       !< Host and device fields.
   real(R8P),    allocatable   :: gx(:), gy(:), gz(:)
   integer(I4P)                :: m, st_h, st_d, fb_h, fb_d
   integer(I8P)                :: c0, c1, rate
   real(R8P)                   :: t_host, t_dev, ts_host, ts_dev
   logical                     :: band_ok

   allocate(gx(grid_n), gy(grid_n), gz(grid_n))
   do m = 1_I4P, grid_n
      frac  = real(m - 1_I4P, R8P) / real(grid_n - 1_I4P, R8P)
      gx(m) = lo(1) + frac * (hi(1) - lo(1))
      gy(m) = lo(2) + frac * (hi(2) - lo(2))
      gz(m) = lo(3) + frac * (hi(3) - lo(3))
   enddo
   call host%initialize(x=gx, y=gy, z=gz)
   call dev%initialize(x=gx, y=gy, z=gz)
   call system_clock(c0, rate)
   call host%compute_band(surface=surface, sign_algorithm=SIGN_PSEUDO_NORMAL, status=st_h)
   call system_clock(c1)
   t_host = real(c1 - c0, R8P) / real(rate, R8P)
   call system_clock(c0, rate)
   call dev%compute_band(surface=surface, sign_algorithm=SIGN_PSEUDO_NORMAL, use_device=.true., status=st_d)
   call system_clock(c1)
   t_dev = real(c1 - c0, R8P) / real(rate, R8P)
   band_ok = st_h == STATUS_OK .and. st_d == STATUS_OK .and. host%get_band_cells() > 0_I4P .and. &
        all(host%phi == dev%phi) .and. all(host%in_band .eqv. dev%in_band)
   ! signs: the point queries follow the band's device flag (open meshes query most cells)
   call system_clock(c0, rate)
   call host%compute_signs(surface=surface, fallback_cells=fb_h, status=st_h)
   call system_clock(c1)
   ts_host = real(c1 - c0, R8P) / real(rate, R8P)
   call system_clock(c0, rate)
   call dev%compute_signs(surface=surface, fallback_cells=fb_d, status=st_d)
   call system_clock(c1)
   ts_dev = real(c1 - c0, R8P) / real(rate, R8P)
   ok = band_ok .and. st_h == STATUS_OK .and. st_d == STATUS_OK .and. fb_h == fb_d .and. all(host%phi == dev%phi)
   print '(A)',          'distance-field band (issue #32, PN sign, probe box):'
   print '(A,I0,A,ES10.3)', '   band cells:             ', host%get_band_cells(), '   band ', host%get_band_width()
   print '(A,F12.4,A)',  '   host band wall-clock:   ', 1000._R8P * t_host, ' ms'
   print '(A,F12.4,A)',  '   device band wall-clock: ', 1000._R8P * t_dev, ' ms (incl. enter/exit_device)'
   print '(A,L1)',       '   device band == host band: ', band_ok
   print '(A,I0)',       '   sign fallback cells:    ', fb_h
   print '(A,F12.4,A)',  '   host signs wall-clock:  ', 1000._R8P * ts_host, ' ms'
   print '(A,F12.4,A)',  '   device signs wall-clock:', 1000._R8P * ts_dev, ' ms (incl. enter/exit_device)'
   print '(A,L1)',       '   device field == host field (band + signs): ', ok
   print '(A)',          ''
   endsubroutine field_band_case

   subroutine holed_device_case(ok)
   !< FWN sign on the device for a surface with a hole (issue #30): cube-02 with
   !< its first facet removed, probed on a 16^3 grid spanning 1.5x its bbox.
   !< Pass when the device FWN result is bit-exact against the serial CPU FWN
   !< result and the case is not vacuous (some probes inside, the centroid
   !< among them: seen from the centre the missing facet removes far less than
   !< half of the 4 pi solid angle). The pseudo-normal sign mismatch count
   !< against FWN is printed only as context, since PN is undefined on a hole.
   logical, intent(out)            :: ok                  !< Pass flag.
   type(surface_stl_object)        :: cube, holed         !< Source and holed surfaces.
   type(facet_object), allocatable :: facets(:)           !< Facets minus the first.
   type(facet_object), pointer     :: f                   !< Facet cursor.
   type(vector_R8P), allocatable   :: pts(:)              !< Probe points (centroid first).
   real(R8P),        allocatable   :: d_cpu(:), d_dev(:)  !< FWN serial / device.
   real(R8P),        allocatable   :: d_pn(:)             !< PN device (context only).
   type(vector_R8P)                :: c, lo_v, hi_v       !< Centroid and grid corners.
   real(R8P)                       :: err                 !< Max |device - serial|.
   integer(I4P), parameter         :: NG = 16_I4P         !< Grid size per axis.
   integer(I4P)                    :: ii, jj, kk, m, n_in, n_pn_mis, st

   call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
   allocate(facets(cube%get_facets_number() - 1_I4P))
   do ii = 2_I4P, cube%get_facets_number()
      f => cube%facet_at(ii)
      facets(ii - 1_I4P) = f
   enddo
   call holed%adopt_facets(facets=facets)

   lo_v = holed%get_bmin() ; hi_v = holed%get_bmax()
   c%x = 0.5_R8P * (lo_v%x + hi_v%x) ; c%y = 0.5_R8P * (lo_v%y + hi_v%y) ; c%z = 0.5_R8P * (lo_v%z + hi_v%z)
   allocate(pts(1_I4P + NG * NG * NG), d_cpu(1_I4P + NG * NG * NG), d_dev(1_I4P + NG * NG * NG), d_pn(1_I4P + NG * NG * NG))
   pts(1) = c
   m = 1_I4P
   do kk = 0_I4P, NG - 1_I4P
      do jj = 0_I4P, NG - 1_I4P
         do ii = 0_I4P, NG - 1_I4P
            m = m + 1_I4P
            pts(m)%x = c%x + 1.5_R8P * (hi_v%x - lo_v%x) * (real(ii, R8P) / real(NG - 1_I4P, R8P) - 0.5_R8P)
            pts(m)%y = c%y + 1.5_R8P * (hi_v%y - lo_v%y) * (real(jj, R8P) / real(NG - 1_I4P, R8P) - 0.5_R8P)
            pts(m)%z = c%z + 1.5_R8P * (hi_v%z - lo_v%z) * (real(kk, R8P) / real(NG - 1_I4P, R8P) - 0.5_R8P)
         enddo
      enddo
   enddo

   do m = 1_I4P, size(pts, kind=I4P)
      d_cpu(m) = holed%distance(point=pts(m), is_signed=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER, &
                                is_square_root=.true.)
   enddo
   call holed%enter_device
   call holed%distance_many_device(points=pts, distances=d_dev, is_signed=.true., &
                                   sign_algorithm=SIGN_FAST_WINDING_NUMBER, is_square_root=.true., status=st)
   call holed%distance_many_device(points=pts, distances=d_pn, is_signed=.true., is_square_root=.true.)
   call holed%exit_device

   err = maxval(abs(d_dev - d_cpu))
   n_in = count(d_cpu < 0._R8P)
   n_pn_mis = count((d_pn < 0._R8P) .neqv. (d_cpu < 0._R8P))
   ok = (st == STATUS_OK .and. err == 0._R8P .and. n_in > 0_I4P .and. d_cpu(1) < 0._R8P .and. d_dev(1) < 0._R8P)

   print '(A)',         'holed mesh on the device (cube-02 minus one facet, FWN sign, issue #30):'
   print '(A,L1)',      '   watertight:             ', holed%is_watertight()
   print '(A,I0,A,I0)', '   inside (serial FWN):    ', n_in, ' of ', size(pts)
   print '(A,L1)',      '   centroid inside:        ', d_dev(1) < 0._R8P
   print '(A,ES10.3)',  '   fwn device vs serial max |d| err: ', err
   print '(A,I0)',      '   PN sign mismatches vs FWN (context): ', n_pn_mis
   print '(A,L1)',      '   holed case ok:          ', ok
   print '(A)',         ''
   endsubroutine holed_device_case
endprogram fossil_test_distance_bench
