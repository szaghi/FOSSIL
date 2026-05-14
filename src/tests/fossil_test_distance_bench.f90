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
!< Usage: fossil_test_distance_bench [<stl_path> [<grid_n>]]
!<        default <stl_path> = src/tests/dragon-fine.stl
!<        default <grid_n>   = 40
!<
!< `perf stat` companion (baseline cache-miss / IPC capture for issue #19):
!<   perf stat -e cache-references,cache-misses,instructions,cycles \
!<       ./exe/fossil_test_distance_bench src/tests/dragon-fine.stl 40

program fossil_test_distance_bench

use fossil,              only : surface_stl_object, SIGN_PSEUDO_NORMAL
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

type(surface_stl_object)      :: surface
type(vector_R8P)              :: bmin, bmax, span, probe
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_signed(:), d_unsigned(:)
real(R8P)                     :: t0, t1, t_signed, t_unsigned
real(R8P)                     :: lo(3), hi(3), frac
real(R8P)                     :: d_tree, d_brute, max_abs_err
integer(I4P)                  :: grid_n, n_queries, nf
integer(I4P)                  :: i, j, k, q
integer(I4P)                  :: argc, n_sample, s, stride
character(len=512)            :: stl_path, arg2
logical                       :: correctness_ok
real(R8P), parameter          :: BBOX_INFLATE = 1.5_R8P    !< Probe grid spans 1.5x the surface bbox.
integer(I4P), parameter       :: N_SAMPLE_MAX = 200_I4P    !< Sub-sample size for the brute-force correctness check.

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
print '(A,I0,A,ES10.3)', 'correctness (BVH vs brute force, ', s, ' samples): max |d_tree - d_brute| = ', max_abs_err
print '(A)',         ''
print '(A,L1)',      'Are all tests passed? ', correctness_ok

endprogram fossil_test_distance_bench
