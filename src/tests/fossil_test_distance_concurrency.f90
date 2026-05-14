!< FOSSIL test: concurrency-safety of the distance query path.
!<
!< Pins the contract that the ADAM immersed-boundary AMR-marker pass depends
!< on: many threads may call `surface%distance` / `surface%compute_distance`
!< concurrently on a single, **shared, read-only** `surface_stl_object` and
!< every call must return bit-exactly what the serial path returns.
!<
!< Why this test exists: `compute_distance` and its whole call tree are
!< `intent(in)` on `self` / `pure`, with no module `save` state on the query
!< path — so the parallelism is safe *by construction*. But "safe by
!< construction" is an invariant that a future change could silently break
!< (a lazy cache written into `self`, a module-level scratch buffer, a
!< `save` variable). A passing parallel *benchmark* is weak evidence —
!< races are nondeterministic and a benchmark only checks timing. This test
!< checks **correctness under contention** and is the regression guard that
!< catches the day someone sneaks mutable state into the distance path.
!<
!< Strategy — make a race, if one exists, *likely* to surface:
!<   1. Build one surface; compute a golden reference serially.
!<   2. Re-run the identical query set inside `!$omp parallel do
!<      schedule(dynamic, 1)` — chunk-1 dynamic scheduling maximises
!<      thread interleaving (threads grab single queries, not contiguous
!<      runs), so any shared-state corruption has the widest window to
!<      manifest.
!<   3. Repeat the parallel sweep `N_SWEEPS` times — races are
!<      nondeterministic; one clean run proves little.
!<   4. Exercise every sign path: unsigned, signed pseudo-normal (touches
!<      `pseudo_normal_for_region` + the closest-facet recompute), signed
!<      ray-intersections and signed solid-angle (both go through
!<      `is_point_inside`). The contract must cover the whole query
!<      surface, not just the default.
!<   5. Assert **bit-exact** equality parallel-vs-serial — any race in a
!<      read-only-claimed path would perturb a result.
!<
!< Build/run:
!<   - Under `tests-gnu-openmp` the `!$omp` directives activate and the
!<     test genuinely runs multi-threaded (set `OMP_NUM_THREADS`).
!<   - Under plain `tests-gnu` the directives are inert comments; the test
!<     still builds and runs (serially) and still passes — it then only
!<     proves the serial path is self-consistent, which is a weaker but
!<     still valid check. The real assertion happens in the OpenMP build.
!<
!< Usage: fossil_test_distance_concurrency [<stl_path>]
!<        default <stl_path> = src/tests/dragon.stl

program fossil_test_distance_concurrency

use fossil, only : surface_stl_object,                                  &
                   SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object)      :: surface
type(vector_R8P)              :: bmin, bmax, span, probe
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: ref_unsigned(:), ref_pseudo(:), ref_ray(:), ref_sa(:)
real(R8P),        allocatable :: par_unsigned(:), par_pseudo(:), par_ray(:), par_sa(:)
real(R8P)                     :: lo(3), hi(3), frac
integer(I4P)                  :: grid_n, n_queries, nf
integer(I4P)                  :: i, j, k, q, sweep
integer(I4P)                  :: argc
integer(I4P)                  :: n_threads
integer(I4P)                  :: mism_unsigned, mism_pseudo, mism_ray, mism_sa
character(len=512)            :: stl_path
logical                       :: are_tests_passed(4)
real(R8P), parameter          :: BBOX_INFLATE  = 1.5_R8P !< Probe grid spans 1.5x the surface bbox.
integer(I4P), parameter       :: GRID_EDGE     = 20_I4P  !< 20^3 = 8000 queries — enough contention, fast enough to repeat.
integer(I4P), parameter       :: N_SWEEPS      = 8_I4P   !< Parallel sweeps; races are nondeterministic, so repeat.
!$ integer, external          :: omp_get_max_threads

! ---- arguments ---------------------------------------------------------------
argc = command_argument_count()
if (argc >= 1) then
   call get_command_argument(1, stl_path)
else
   stl_path = 'src/tests/dragon.stl'
endif

n_threads = 1_I4P
!$ n_threads = omp_get_max_threads()

! ---- load + freeze the surface ----------------------------------------------
! From here on `surface` is READ-ONLY. Every `distance` call below — serial or
! parallel — only queries it; nothing mutates it. This is exactly the contract
! the ADAM IB marker pass must honour.
call surface%load_from_file(file_name=trim(stl_path), guess_format=.true.)
nf = surface%get_facets_number()

! ---- deterministic probe grid -----------------------------------------------
grid_n = GRID_EDGE
bmin = surface%get_bmin()
bmax = surface%get_bmax()
span = bmax - bmin
lo(1) = bmin%x - 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%x
lo(2) = bmin%y - 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%y
lo(3) = bmin%z - 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%z
hi(1) = bmax%x + 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%x
hi(2) = bmax%y + 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%y
hi(3) = bmax%z + 0.5_R8P * (BBOX_INFLATE - 1._R8P) * span%z

n_queries = grid_n * grid_n * grid_n
allocate(points(n_queries))
allocate(ref_unsigned(n_queries), ref_pseudo(n_queries), ref_ray(n_queries), ref_sa(n_queries))
allocate(par_unsigned(n_queries), par_pseudo(n_queries), par_ray(n_queries), par_sa(n_queries))

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

! ---- golden reference: serial sweep over every sign path --------------------
do q = 1_I4P, n_queries
   ref_unsigned(q) = surface%distance(point=points(q), is_square_root=.true.)
   ref_pseudo(q)   = surface%distance(point=points(q), is_signed=.true., &
                                      sign_algorithm=SIGN_PSEUDO_NORMAL, is_square_root=.true.)
   ref_ray(q)      = surface%distance(point=points(q), is_signed=.true., &
                                      sign_algorithm=SIGN_RAY_INTERSECTIONS, is_square_root=.true.)
   ref_sa(q)       = surface%distance(point=points(q), is_signed=.true., &
                                      sign_algorithm=SIGN_SOLID_ANGLE, is_square_root=.true.)
enddo

! ---- contended sweeps: same queries, same shared surface, many threads ------
mism_unsigned = 0_I4P
mism_pseudo   = 0_I4P
mism_ray      = 0_I4P
mism_sa       = 0_I4P

do sweep = 1_I4P, N_SWEEPS
   ! chunk-1 dynamic scheduling: threads grab single queries, maximising
   ! interleaving so any shared-state corruption has the widest window.
   !$omp parallel do default(none) schedule(dynamic, 1) &
   !$omp   shared(surface, points, n_queries, par_unsigned, par_pseudo, par_ray, par_sa) &
   !$omp   private(q)
   do q = 1_I4P, n_queries
      par_unsigned(q) = surface%distance(point=points(q), is_square_root=.true.)
      par_pseudo(q)   = surface%distance(point=points(q), is_signed=.true., &
                                         sign_algorithm=SIGN_PSEUDO_NORMAL, is_square_root=.true.)
      par_ray(q)      = surface%distance(point=points(q), is_signed=.true., &
                                         sign_algorithm=SIGN_RAY_INTERSECTIONS, is_square_root=.true.)
      par_sa(q)       = surface%distance(point=points(q), is_signed=.true., &
                                         sign_algorithm=SIGN_SOLID_ANGLE, is_square_root=.true.)
   enddo
   !$omp end parallel do

   ! Bit-exact check after every sweep — a race that fires on sweep 5 but not
   ! sweeps 1-4 still gets caught.
   do q = 1_I4P, n_queries
      if (par_unsigned(q) /= ref_unsigned(q)) mism_unsigned = mism_unsigned + 1_I4P
      if (par_pseudo(q)   /= ref_pseudo(q))   mism_pseudo   = mism_pseudo   + 1_I4P
      if (par_ray(q)      /= ref_ray(q))      mism_ray      = mism_ray      + 1_I4P
      if (par_sa(q)       /= ref_sa(q))       mism_sa       = mism_sa       + 1_I4P
   enddo
enddo

! ---- report ------------------------------------------------------------------
print '(A,A)',         '--- FOSSIL distance concurrency-safety test: ', trim(stl_path)
print '(A,I0)',        'facets:                       ', nf
print '(A,I0,A,I0)',   'probe grid:                   ', grid_n, ' ^3  ->  queries = ', n_queries
print '(A,I0)',        'parallel sweeps:              ', N_SWEEPS
print '(A,I0)',        'OpenMP threads:               ', n_threads
if (n_threads == 1_I4P) then
   print '(A)',        'NOTE: 1 thread — built without -fopenmp, or OMP_NUM_THREADS=1.'
   print '(A)',        '      This run only proves serial self-consistency. The real'
   print '(A)',        '      contention check needs the tests-gnu-openmp build with'
   print '(A)',        '      OMP_NUM_THREADS > 1.'
endif
print '(A)',           ''
print '(A,I0)',        'mismatches vs serial — unsigned:            ', mism_unsigned
print '(A,I0)',        'mismatches vs serial — signed pseudo-normal: ', mism_pseudo
print '(A,I0)',        'mismatches vs serial — signed ray-intersect: ', mism_ray
print '(A,I0)',        'mismatches vs serial — signed solid-angle:   ', mism_sa
print '(A)',           ''

are_tests_passed(1) = (mism_unsigned == 0_I4P)
are_tests_passed(2) = (mism_pseudo   == 0_I4P)
are_tests_passed(3) = (mism_ray      == 0_I4P)
are_tests_passed(4) = (mism_sa       == 0_I4P)

print '(A,4L2)',       'per-path results: ', are_tests_passed
print '(A,L1)',        'Are all tests passed? ', all(are_tests_passed)

endprogram fossil_test_distance_concurrency
