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
!<      `line crossings device vs serial mismatches` (issue #40) must be 0, and so must every later device line
!<      (`orient3d device vs serial mismatches`, issue #100 step 3c, among them).
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
use fossil_sign_bench_bodies, only : body_t, bench_bodies, device_sign_mismatches_3d, device_line_mismatches_3d, to_surface
use fossil,              only : rigid_pose_t, pose_from_axis_angle, to_world
use fossil,              only : curvilinear_block_t, detect_topology, topo_record_t, TOPO_BOUNDARY, wall_cutter
use fossil_overset,      only : hole_flags, block_logical_t, donor_set_t, donor_candidates
use fossil,              only : curve_object, slice_surface, cut_cells_3d_t, cut_cells_3d, cut_cells_3d_device
use fossil_boolean_bodies, only : BODY_SPHERE, mc_body
use fossil_cut_cells_3d, only : cut_cells_3d_sections, cut_cells_3d_nodes_device
use fossil_cut_cells_compare, only : cut_cells_many_mismatches, cut_cells_3d_face_mismatches, cut_cells_3d_mismatches
use fossil_orient3d_cases, only : orient3d_quadruples, lcb_cases, orient3d_many, lcb_many
use fossil_predicates,   only : orient3d_xyz
use fossil_line_query,   only : line_crossing_before

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
integer                       :: block_mismatch            !< Block location: device vs CPU mismatches (issue #71).
integer                       :: hole_mismatch             !< Hole flags: device vs CPU mismatches (issue #72).
integer                       :: donor_mismatch            !< Donor candidates: device vs CPU mismatches (issue #72).
integer                       :: sections_mismatch         !< 3-D cut-cell sections: device vs CPU mismatches (issue #100).
integer                       :: orient3d_mismatch         !< orient3d, line_crossing_before: device vs CPU (#100 step 3c).
integer                       :: nodes_mismatch            !< 3-D cut-cell node states, edge fractions: device vs CPU (3c).
integer                       :: cells_mismatch            !< 3-D cut cells: device vs CPU mismatches (#100 step 3c).
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
call block_case(mismatch=block_mismatch, hole_mismatch=hole_mismatch, donor_mismatch=donor_mismatch)
rich_ok = rich_ok .and. block_mismatch == 0 .and. hole_mismatch == 0 .and. donor_mismatch == 0
call sections_case(mismatch=sections_mismatch)
rich_ok = rich_ok .and. sections_mismatch == 0
call orient3d_case(mismatch=orient3d_mismatch)
rich_ok = rich_ok .and. orient3d_mismatch == 0
call nodes_case(mismatch=nodes_mismatch)
rich_ok = rich_ok .and. nodes_mismatch == 0
call cells_case(mismatch=cells_mismatch)
rich_ok = rich_ok .and. cells_mismatch == 0

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
print '(A,I0)',      '   block location device vs serial mismatches: ', block_mismatch
print '(A,I0)',      '   hole flags device vs serial mismatches: ', hole_mismatch
print '(A,I0)',      '   donor candidates device vs serial mismatches: ', donor_mismatch
print '(A,I0)',      '   cut cells 3-D sections device vs serial mismatches: ', sections_mismatch
print '(A,I0)',      '   orient3d device vs serial mismatches: ', orient3d_mismatch
print '(A,I0)',      '   cut cells 3-D nodes device vs serial mismatches: ', nodes_mismatch
print '(A,I0)',      '   cut cells 3-D device vs serial mismatches: ', cells_mismatch
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

   subroutine block_case(mismatch, hole_mismatch, donor_mismatch)
   !< Point location in a curvilinear block on the device against the CPU (issue #71): a wavy 3-D O-grid of
   !< 65 x 17 x 17 vertices with its seam, 50000 points in its box; cells, xi (bit for bit) and statuses. Timed.
   integer, intent(out)      :: mismatch               !< Mismatches.
   integer, intent(out)      :: hole_mismatch          !< Hole flag mismatches (issue #72).
   integer, intent(out)      :: donor_mismatch         !< Donor candidate mismatches (issue #72).
   type(donor_set_t), allocatable :: d1(:), d2(:)      !< Donor candidates.
   type(curvilinear_block_t) :: blk(1)                 !< Block.
   type(curvilinear_block_t) :: hb(2)                  !< O-grid with its inner wall, background.
   type(surface_stl_object)  :: cut                    !< Cutter.
   type(topo_record_t)       :: rec                    !< Wall record.
   type(block_logical_t), allocatable :: h1(:), h2(:)  !< Hole flags.
   real(R8P),    allocatable :: gx(:,:,:), gy(:,:,:), gz(:,:,:) !< Background.
   real(R8P),    allocatable :: bx(:,:,:), by(:,:,:), bz(:,:,:), pu(:), pv(:), pw(:), x1(:,:), x2(:,:) !< Data.
   integer(I4P), allocatable :: c1(:,:), c2(:,:), s1(:), s2(:) !< Results.
   integer(I4P)              :: i, j, k, n             !< Counters.
   integer(I8P)              :: t0, t1, t2, rate       !< Clock.
   real(R8P)                 :: r, rr(3), pi2          !< Work.

   pi2 = 2 * acos(-1._R8P)
   allocate(bx(65, 17, 17), by(65, 17, 17), bz(65, 17, 17))
   do k = 1, 17
      do j = 1, 17
         do i = 1, 65
            r = (1._R8P + 0.12_R8P * (j - 1)) * (1._R8P + 0.05_R8P * sin(3 * pi2 * (i - 1) / 64 + 0.4_R8P * k))
            bx(i, j, k) = r * cos(pi2 * (i - 1) / 64) ; by(i, j, k) = r * sin(pi2 * (i - 1) / 64)
            bz(i, j, k) = 0.1_R8P * (k - 1) + 0.03_R8P * cos(pi2 * (i - 1) / 64) * j
         enddo
      enddo
   enddo
   bx(65,:,:) = bx(1,:,:) ; by(65,:,:) = by(1,:,:) ; bz(65,:,:) = bz(1,:,:)
   call blk(1)%initialize(bx, by, bz) ; call detect_topology(blk) ; call blk(1)%build_search
   allocate(pu(50000), pv(50000), pw(50000))
   do n = 1, 50000
      call random_number(rr) ; pu(n) = 6.4_R8P * rr(1) - 3.2_R8P ; pv(n) = 6.4_R8P * rr(2) - 3.2_R8P
      pw(n) = 2.1_R8P * rr(3)
   enddo
   call system_clock(t0, rate)
   call blk(1)%locate_many(pu, pv, pw, c1, x1, s1)
   call system_clock(t1)
   call blk(1)%locate_many_device(pu, pv, pw, c2, x2, s2)
   call system_clock(t2)
   mismatch = count(any(c1 /= c2, dim=1) .or. any(x1 /= x2, dim=1) .or. s1 /= s2)
   print '(A,I0,A,I0,A,F9.2,A,F9.2,A)', 'block location: inside ', count(s1 == 1), ' of ', size(s1), ', serial ', &
         1e3_R8P * (t1 - t0) / rate, ' ms, device ', 1e3_R8P * (t2 - t1) / rate, ' ms (with copies)'
   if (count(s1 == 1) < 10000) mismatch = mismatch + 1
   ! hole flags (issue #72): the inner wavy wall of the O-grid (an open tube) cuts a 41^3 background block
   hb(1) = blk(1)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [65, 1, 17] ; rec%name = 'wall'
   call hb(1)%add_record(rec)
   allocate(gx(41, 41, 41), gy(41, 41, 41), gz(41, 41, 41))
   do k = 1, 41
      do j = 1, 41
         do i = 1, 41
            gx(i, j, k) = -3.2_R8P + 0.16_R8P * (i - 1) ; gy(i, j, k) = -3.2_R8P + 0.16_R8P * (j - 1)
            gz(i, j, k) = 0.05_R8P * (k - 1)
         enddo
      enddo
   enddo
   call hb(2)%initialize(gx, gy, gz)
   call wall_cutter(hb(1:1), cut)
   call hole_flags(hb, cut, h1)
   call hole_flags(hb, cut, h2, device=.true.)
   hole_mismatch = count(h1(2)%v .neqv. h2(2)%v) + count(h1(1)%v .neqv. h2(1)%v)
   print '(A,I0,A,I0)', 'hole flags: holes ', count(h1(2)%v), ' of ', size(h1(2)%v)
   if (count(h1(2)%v) < 1000) hole_mismatch = hole_mismatch + 1
   call hb(1)%build_search ; call hb(2)%build_search
   call donor_candidates(hb, d1)
   call hb(1)%enter_device ; call hb(2)%enter_device
   call donor_candidates(hb, d2, device=.true.)
   call hb(1)%exit_device ; call hb(2)%exit_device
   donor_mismatch = 0
   do i = 1, 2
      if (any(d1(i)%start /= d2(i)%start)) then
         donor_mismatch = donor_mismatch + 1 ; cycle
      endif
      donor_mismatch = donor_mismatch + count(any(d1(i)%cell /= d2(i)%cell, dim=1) .or. any(d1(i)%xi /= d2(i)%xi, dim=1) &
                                              .or. d1(i)%block /= d2(i)%block .or. d1(i)%status /= d2(i)%status)
   enddo
   print '(A,I0)', 'donor candidates: ', size(d1(1)%block) + size(d1(2)%block)
   if (size(d1(2)%block) < 1000) donor_mismatch = donor_mismatch + 1
   endsubroutine block_case

   subroutine orient3d_case(mismatch)
   !< The exact orient3d and line_crossing_before on the device (issue #100 step 3c) against the host, bit for bit: the
   !< cases of `fossil_orient3d_cases` (near-coplanar grids at 2^0 and 2^40, exact zeros on shared coordinate planes,
   !< coincident points, inexact near-coplanar points far from the origin, random quadruples; the perturbation cases of
   !< line_crossing_before) in one device loop each against the scalar host calls.
   integer, intent(out)      :: mismatch               !< Mismatches.
   real(R8P),    allocatable :: q(:,:), od(:), ltri(:,:), lab(:,:), lc(:) !< Cases, device results.
   integer(I4P), allocatable :: expect(:), lax(:)      !< Expected signs, line axes.
   logical,      allocatable :: len_(:), lexp(:), lbd(:) !< Entering flags, exact answers, device results.
   real(R8P)                 :: o                      !< Host result.
   integer(I4P)              :: c                      !< Counter.

   mismatch = 0
   call orient3d_quadruples(q, expect)
   allocate(od(size(expect)))
   call orient3d_many(size(expect, kind=I4P), q, od)
   do c = 1, size(expect, kind=I4P)
      o = orient3d_xyz(q(1:3, c), q(4:6, c), q(7:9, c), q(10:12, c))
      if (transfer(o, 0_I8P) /= transfer(od(c), 0_I8P)) mismatch = mismatch + 1
   enddo
   call lcb_cases(ltri, lax, lab, lc, len_, lexp)
   allocate(lbd(size(lexp)))
   call lcb_many(size(lexp, kind=I4P), ltri, lax, lab, lc, len_, lbd)
   do c = 1, size(lexp, kind=I4P)
      if (lbd(c) .neqv. line_crossing_before(reshape(ltri(:, c), [3, 3]), lax(c), lab(1, c), lab(2, c), lc(c), len_(c))) &
         mismatch = mismatch + 1
   enddo
   endsubroutine orient3d_case

   subroutine nodes_case(mismatch)
   !< The node states and cell-edge fractions of the 3-D cut cells on the device (issue #100 step 3c) against
   !< `cut_cells_3d(..., edges=.true.)`, bitwise: the bench body in its probe box cut into grid_n^3 cells, and into
   !< (grid_n + 3)^3 cells of a box with corners at the body's bounding box (grid lines through its extreme vertices);
   !< node_in and edge_x/_y/_z, and node_in alone without the edges.
   integer, intent(out)      :: mismatch               !< Mismatches.
   type(cut_cells_3d_t)      :: ch                     !< Host cut cells.
   logical,   allocatable    :: nd(:,:,:)              !< Device node flags.
   real(R8P), allocatable    :: gx(:), gy(:), gz(:), ex(:,:,:), ey(:,:,:), ez(:,:,:) !< Cell edges, device fractions.
   type(vector_R8P)          :: b0, b1                 !< Bounding box.
   integer(I4P)              :: g, q, st, st2, nc      !< Work.

   mismatch = 0
   do g = 1, 2
      if (g == 1) then
         call section_grid(grid_n, gx, gy, gz)
      else
         b0 = surface%get_bmin() ; b1 = surface%get_bmax() ; nc = grid_n + 3_I4P
         allocate(gx(nc + 1), gy(nc + 1), gz(nc + 1))
         do q = 1, nc + 1
            gx(q) = b0%x + (b1%x - b0%x) * real(q - 1, R8P) / real(nc, R8P)
            gy(q) = b0%y + (b1%y - b0%y) * real(q - 1, R8P) / real(nc, R8P)
            gz(q) = b0%z + (b1%z - b0%z) * real(q - 1, R8P) / real(nc, R8P)
         enddo
      endif
      call cut_cells_3d(surface, gx, gy, gz, ch, status=st, edges=.true.)
      call cut_cells_3d_nodes_device(surface, gx, gy, gz, nd, ex, ey, ez, edges=.true., status=st2)
      mismatch = mismatch + abs(st - st2) + count(ch%node_in .neqv. nd)
      mismatch = mismatch + count(transfer(ch%edge_x, 0_I8P, size(ch%edge_x)) /= transfer(ex, 0_I8P, size(ex)))
      mismatch = mismatch + count(transfer(ch%edge_y, 0_I8P, size(ch%edge_y)) /= transfer(ey, 0_I8P, size(ey)))
      mismatch = mismatch + count(transfer(ch%edge_z, 0_I8P, size(ch%edge_z)) /= transfer(ez, 0_I8P, size(ez)))
      call cut_cells_3d_nodes_device(surface, gx, gy, gz, nd, ex, ey, ez, status=st2)
      mismatch = mismatch + abs(st - st2) + count(ch%node_in .neqv. nd)
      deallocate(gx, gy, gz)
   enddo
   endsubroutine nodes_case

   subroutine cells_case(mismatch)
   !< `cut_cells_3d_device` (issue #100 step 3c) against `cut_cells_3d` without links, every output bitwise
   !< (`cut_cells_3d_mismatches`), with and without the edges, copied back and resident (`update_host`): the bench body
   !< in its probe box at grid_n^3; with grid_n >= 32 also the bodies of step 0 (the marching-cubes sphere, the cone,
   !< cube-02.stl) on m^3 cells over their box grown by 10 % at m = 64 and 128, and a timing table: best and median of
   !< 5 calls after a warm-up, the CPU (`cut_cells_3d`) against the device copied back and resident (the surface mapped
   !< by `enter_device`; the resident call releases the previous one).
   integer, intent(out)      :: mismatch               !< Mismatches.
   type(surface_stl_object)  :: sb(3)                  !< Step-0 bodies.
   type(body_t), allocatable :: bb(:)                  !< Bench bodies.
   character(12), parameter  :: NAMES(3) = [character(12) :: 'mc sphere', 'cone', 'cube-02'] !< Their names.
   real(R8P),    allocatable :: gx(:), gy(:), gz(:)    !< Cell edges.
   integer(I4P)              :: b, g, m                !< Work.

   mismatch = 0
   call section_grid(grid_n, gx, gy, gz)
   mismatch = mismatch + cells_compare(surface, gx, gy, gz)
   if (grid_n < 32_I4P) return
   call mc_body(BODY_SPHERE, 24_I4P, sb(1))
   call bench_bodies(bb) ; call to_surface_bench(bb(7), sb(2))
   call sb(3)%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
   do b = 1, 3
      do g = 1, 2
         m = 64_I4P * g
         call box_grid(sb(b), m, gx, gy, gz)
         mismatch = mismatch + cells_compare(sb(b), gx, gy, gz)
      enddo
   enddo
   print '(A)', 'cut cells 3-D timing (no links; best / median of 5) [s]: body, cells, CPU, device copied back, '// &
                'device resident, speed-ups of the best'
   call section_grid(grid_n, gx, gy, gz)
   call cells_timing('bench body', surface, grid_n, gx, gy, gz)
   do b = 1, 3
      do g = 1, 2
         m = 64_I4P * g
         call box_grid(sb(b), m, gx, gy, gz)
         call cells_timing(NAMES(b), sb(b), m, gx, gy, gz)
      enddo
   enddo
   endsubroutine cells_case

   function cells_compare(s, gx, gy, gz) result(mismatch)
   !< `cut_cells_3d_device` against `cut_cells_3d` on s and the grid gx, gy, gz: with and without the edges, copied back
   !< and resident (then `update_host`; the resident cut overwritten by the next call).
   type(surface_stl_object), intent(in) :: s            !< Body.
   real(R8P),                intent(in) :: gx(:), gy(:), gz(:) !< Cell edges.
   integer                              :: mismatch     !< Mismatches.
   type(cut_cells_3d_t)                 :: ch, ce, cd   !< Host with and without edges, device.
   integer(I4P)                         :: st, st2      !< Status.

   call cut_cells_3d(s, gx, gy, gz, ch, status=st, edges=.true.)
   call cut_cells_3d(s, gx, gy, gz, ce, status=st)
   call cut_cells_3d_device(s, gx, gy, gz, cd, status=st2, edges=.true.)
   mismatch = abs(st - st2) + cut_cells_3d_mismatches(ch, cd)
   call cut_cells_3d_device(s, gx, gy, gz, cd, status=st2)
   mismatch = mismatch + abs(st - st2) + cut_cells_3d_mismatches(ce, cd)
   call cut_cells_3d_device(s, gx, gy, gz, cd, status=st2, resident=.true.)
   call cd%update_host
   mismatch = mismatch + abs(st - st2) + cut_cells_3d_mismatches(ce, cd)
   call cut_cells_3d_device(s, gx, gy, gz, cd, status=st2, edges=.true., resident=.true.)
   call cd%update_host
   mismatch = mismatch + abs(st - st2) + cut_cells_3d_mismatches(ch, cd)
   call cd%exit_device
   endfunction cells_compare

   subroutine cells_timing(label, s, m, gx, gy, gz)
   !< One row of the timing table of `cells_case`.
   character(*),             intent(in)    :: label        !< Body.
   type(surface_stl_object), intent(inout) :: s            !< Body (mapped for the device runs, then released).
   integer(I4P),             intent(in)    :: m            !< Cells per axis.
   real(R8P),                intent(in)    :: gx(:), gy(:), gz(:) !< Cell edges.
   type(cut_cells_3d_t)                    :: ch, cd       !< Host, device.
   real(R8P)                               :: t(3, 5)      !< Times.
   integer(I8P)                            :: c0, c1, rate !< Clock.
   integer(I4P)                            :: r            !< Counter.

   call system_clock(count_rate=rate)
   call s%enter_device
   do r = 1, 3 ! warm-up (issue #100 step 6a: slow pageable copies after a host-heavy phase)
      call cut_cells_3d_device(s, gx, gy, gz, cd)
   enddo
   do r = 1, 5
      call system_clock(c0) ; call cut_cells_3d(s, gx, gy, gz, ch) ; call system_clock(c1)
      t(1, r) = real(c1 - c0, R8P) / rate
      call system_clock(c0) ; call cut_cells_3d_device(s, gx, gy, gz, cd) ; call system_clock(c1)
      t(2, r) = real(c1 - c0, R8P) / rate
   enddo
   do r = 1, 5
      call system_clock(c0) ; call cut_cells_3d_device(s, gx, gy, gz, cd, resident=.true.) ; call system_clock(c1)
      t(3, r) = real(c1 - c0, R8P) / rate
   enddo
   call cd%exit_device
   call s%exit_device
   print '(3X,A12,I4,A,3(2X,ES9.2,A,ES9.2),A,F6.2,A,F6.2,A)', label, m, '^3', minval(t(1, :)), ' /', &
         median_of5(t(1, :)), minval(t(2, :)), ' /', median_of5(t(2, :)), minval(t(3, :)), ' /', median_of5(t(3, :)), &
         '  (', minval(t(1, :)) / minval(t(2, :)), 'x, ', minval(t(1, :)) / minval(t(3, :)), 'x)'
   endsubroutine cells_timing

   subroutine box_grid(s, m, gx, gy, gz)
   !< m^3 uniform cells over the box of s grown by 10 % (the grid of the perf bench's 3-D cut cells).
   type(surface_stl_object), intent(in)  :: s              !< Body.
   integer(I4P),             intent(in)  :: m              !< Cells per axis.
   real(R8P), allocatable,   intent(out) :: gx(:), gy(:), gz(:) !< Cell edges.
   type(vector_R8P)                      :: b0, b1         !< Box.
   integer(I4P)                          :: q              !< Counter.

   b0 = s%get_bmin() ; b1 = s%get_bmax()
   gx = [(b0%x - 0.1_R8P * (b1%x - b0%x) + 1.2_R8P * (b1%x - b0%x) * q / m, q = 0, m)]
   gy = [(b0%y - 0.1_R8P * (b1%y - b0%y) + 1.2_R8P * (b1%y - b0%y) * q / m, q = 0, m)]
   gz = [(b0%z - 0.1_R8P * (b1%z - b0%z) + 1.2_R8P * (b1%z - b0%z) * q / m, q = 0, m)]
   endsubroutine box_grid

   subroutine to_surface_bench(b, s)
   !< The bench body b as a surface, centred on its box (the perf bench's cone).
   type(body_t),             intent(in)    :: b      !< Body.
   type(surface_stl_object), intent(inout) :: s      !< Surface.
   type(vector_R8P)                        :: b0, b1 !< Box.

   call to_surface(b%soup, s)
   b0 = s%get_bmin() ; b1 = s%get_bmax()
   call s%translate(delta=-0.5_R8P * (b0 + b1)) ; call s%analyze
   endsubroutine to_surface_bench

   subroutine sections_case(mismatch)
   !< The plane sections of `cut_cells_3d` on the device (issue #100 step 3b) against the host path, bitwise, on the
   !< bench body and its probe box cut into grid_n^3 cells: `cut_cells_3d(..., device_sections=.true.)` against
   !< `cut_cells_3d` (face arrays, cell fractions, pieces), and the z sections as the blocks of `cut_cells_device_many`
   !< with pieces and intervals against `curve%cut_cells` per section (`cut_cells_many_mismatches`, copied back and
   !< resident). With grid_n >= 32 also the faces alone (`cut_cells_3d_sections`) at 64^3 and 128^3, and a timing table
   !< of the sections: best and median of 5 calls after a warm-up, host vs device (copies included).
   integer, intent(out)      :: mismatch               !< Mismatches.
   type(cut_cells_3d_t)      :: ch, cd                 !< Host, device sections.
   type(curve_object), allocatable :: secs(:)          !< z sections.
   real(R8P),    allocatable :: gx(:), gy(:), gz(:)    !< Cell edges.
   real(R8P)                 :: t(2, 5)                !< Times.
   integer(I4P)              :: ng(3), g, m, r, st, st2 !< Work.
   integer(I8P)              :: c0, c1, rate           !< Clock.

   mismatch = 0
   call section_grid(grid_n, gx, gy, gz)
   call cut_cells_3d(surface, gx, gy, gz, ch, status=st)
   call cut_cells_3d(surface, gx, gy, gz, cd, status=st2, device_sections=.true.)
   mismatch = mismatch + cut_cells_3d_face_mismatches(ch, cd) + abs(st - st2)
   mismatch = mismatch + count(transfer(ch%alpha, 0_I8P, size(ch%alpha)) /= transfer(cd%alpha, 0_I8P, size(cd%alpha)))
   mismatch = mismatch + count(ch%pieces /= cd%pieces)
   if (ch%np /= cd%np) mismatch = mismatch + 1
   allocate(secs(size(gz)))
   do g = 1, size(gz)
      call slice_surface(surface=surface, normal_axis=3_I4P, level=gz(g), curve=secs(g), on_plane_low=.true.)
   enddo
   mismatch = mismatch + cut_cells_many_mismatches(secs, gx, gy, st) + abs(st)
   deallocate(secs)
   if (grid_n < 32_I4P) return
   ng = [grid_n, 64_I4P, 128_I4P]
   do g = 2, 3
      call section_grid(ng(g), gx, gy, gz)
      call cut_cells_3d_sections(surface, gx, gy, gz, ch, status=st)
      call cut_cells_3d_sections(surface, gx, gy, gz, cd, device=.true., status=st2)
      mismatch = mismatch + cut_cells_3d_face_mismatches(ch, cd) + abs(st - st2)
   enddo
   print '(A)', 'cut cells 3-D sections timing (3(N+1) planes, faces copied back; best / median of 5) [s]: cells, host, '// &
                'device, speed-up of the best'
   call system_clock(count_rate=rate)
   do r = 1, 5 ! warm-up (issue #100 step 6a: slow pageable copies after a host-heavy phase)
      call cut_cells_3d_sections(surface, gx, gy, gz, cd, device=.true.)
   enddo
   do g = 1, 3
      if (g > 1 .and. ng(g) == ng(1)) cycle
      call section_grid(ng(g), gx, gy, gz)
      do r = 1, 5
         call system_clock(c0) ; call cut_cells_3d_sections(surface, gx, gy, gz, ch) ; call system_clock(c1)
         t(1, r) = real(c1 - c0, R8P) / rate
         call system_clock(c0) ; call cut_cells_3d_sections(surface, gx, gy, gz, cd, device=.true.) ; call system_clock(c1)
         t(2, r) = real(c1 - c0, R8P) / rate
      enddo
      m = ng(g)
      print '(3X,I4,A,2(2X,ES9.2,A,ES9.2),A,F6.2,A)', m, '^3', minval(t(1, :)), ' /', median_of5(t(1, :)), &
            minval(t(2, :)), ' /', median_of5(t(2, :)), '  (', minval(t(1, :)) / minval(t(2, :)), 'x)'
   enddo
   endsubroutine sections_case

   subroutine section_grid(nc, gx, gy, gz)
   !< The probe box [lo, hi] cut into nc^3 uniform cells.
   integer(I4P),           intent(in)  :: nc             !< Cells per axis.
   real(R8P), allocatable, intent(out) :: gx(:), gy(:), gz(:) !< Cell edges.
   integer(I4P)                        :: q              !< Counter.

   allocate(gx(nc + 1), gy(nc + 1), gz(nc + 1))
   do q = 1, nc + 1
      gx(q) = lo(1) + (hi(1) - lo(1)) * real(q - 1, R8P) / real(nc, R8P)
      gy(q) = lo(2) + (hi(2) - lo(2)) * real(q - 1, R8P) / real(nc, R8P)
      gz(q) = lo(3) + (hi(3) - lo(3)) * real(q - 1, R8P) / real(nc, R8P)
   enddo
   endsubroutine section_grid

   pure function median_of5(a) result(m)
   !< Median of five values.
   real(R8P), intent(in) :: a(5) !< Values.
   real(R8P)             :: m    !< Median.
   real(R8P)             :: b(5), w !< Sorted copy, swap.
   integer               :: i, j !< Counters.

   b = a
   do i = 2, 5
      w = b(i) ; j = i - 1
      do while (j >= 1)
         if (b(j) <= w) exit
         b(j + 1) = b(j) ; j = j - 1
      enddo
      b(j + 1) = w
   enddo
   m = b(3)
   endfunction median_of5

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
