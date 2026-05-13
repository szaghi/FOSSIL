!< FOSSIL test: QEM mesh decimation (issue #18 §1.3).
!<
!< Three invariants on analytic and small inputs:
!<
!<   1. Sphere via §1.5 marching cubes → decimate to 50% / 25% / 10%:
!<      assert facet count ≤ target; volume preserved within ratio-dependent
!<      tolerance (10% / 20% / 30%); manifoldness preserved (zero
!<      non-manifold edges introduced). The tolerances are MVP-realistic for
!<      a QEM implementation without tie-breaking heuristics — production
!<      decimators add a length penalty for ties and an optional volume
!<      preservation term, both deferred for §1.3 follow-up.
!<      Note: `is_watertight` is NOT asserted because the upstream MC sphere
!<      itself is not watertight (vertex-pool dedup tolerance vs MC's
!<      per-edge interpolated vertex jitter — separate concern from §1.3).
!<
!<   2. Cube.stl → decimate to its current count (no-op): assert
!<      facet_count and volume unchanged. Sanity check on the trivial path.
!<
!<   3. Cube.stl → decimate to 4 facets (unreachable for the cube without
!<      severe geometric distortion): assert no crash. The QEM cost is zero
!<      everywhere on a perfectly-flat-faced input (planar quadrics evaluate
!<      to zero at all cube vertices), so the algorithm picks arbitrarily
!<      and produces a heavily distorted result. Documenting that the MVP
!<      handles the case without crashing is the assertion of value here.

program fossil_test_decimate

use fossil, only : surface_stl_object, &
                   extract_isosurface, &
                   DEC_STATUS_OK, DEC_STATUS_NO_PROGRESS
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: sphere, cube
type(facet_object), allocatable :: facets(:)
real(R8P), allocatable   :: values(:, :, :)
type(vector_R8P)         :: bmin, bmax
real(R8P)                :: vol, vol_target, vol_before
integer(I4P)             :: status, n0, n_target, n_before, i, j, k
real(R8P)                :: x, y, z, dx
logical                  :: are_tests_passed(3)

integer(I4P), parameter :: N_GRID = 32_I4P
real(R8P),    parameter :: PI     = 3.141592653589793_R8P

are_tests_passed = .false.

! ----- 1. Sphere SDF → MC → decimate, checking volume + manifoldness.
dx = 4._R8P / real(N_GRID - 1, R8P)
allocate(values(N_GRID, N_GRID, N_GRID))
do k = 1, N_GRID
   z = -2._R8P + (k - 1) * dx
   do j = 1, N_GRID
      y = -2._R8P + (j - 1) * dx
      do i = 1, N_GRID
         x = -2._R8P + (i - 1) * dx
         values(i, j, k) = sqrt(x**2 + y**2 + z**2) - 1._R8P
      enddo
   enddo
enddo
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=facets, status=status)
call sphere%adopt_facets(facets=facets)
deallocate(values)
n0         = sphere%get_facets_number()
vol_target = 4._R8P / 3._R8P * PI

! Decimate to 50%; assert facet count and volume.
n_target = n0 / 2
call sphere%decimate(target_facets=n_target, status=status)
vol = sphere%get_volume()
print '(A,I0,A,I0,A,I0,A,F8.4,A,F8.4,A,I0)', &
      'sphere 50%: n0=', n0, ' target=', n_target, ' got=', sphere%get_facets_number(), &
      '  vol=', vol, ' target=', vol_target, ' nm_edges=', sphere%get_non_manifold_edges_number()
if (.not. (status == DEC_STATUS_OK .and. &
           sphere%get_facets_number() <= n_target .and. &
           abs(vol - vol_target) <= 0.10_R8P * vol_target .and. &
           sphere%get_non_manifold_edges_number() == 0_I4P)) goto 100

n_target = n0 / 4
call sphere%decimate(target_facets=n_target, status=status)
vol = sphere%get_volume()
print '(A,I0,A,I0,A,F8.4,A,I0)', 'sphere 25%: target=', n_target, ' got=', sphere%get_facets_number(), &
      '  vol=', vol, ' nm_edges=', sphere%get_non_manifold_edges_number()
if (.not. (status == DEC_STATUS_OK .and. &
           sphere%get_facets_number() <= n_target .and. &
           abs(vol - vol_target) <= 0.20_R8P * vol_target .and. &
           sphere%get_non_manifold_edges_number() == 0_I4P)) goto 100

n_target = n0 / 10
call sphere%decimate(target_facets=n_target, status=status)
vol = sphere%get_volume()
print '(A,I0,A,I0,A,F8.4,A,I0)', 'sphere 10%: target=', n_target, ' got=', sphere%get_facets_number(), &
      '  vol=', vol, ' nm_edges=', sphere%get_non_manifold_edges_number()
if (.not. ((status == DEC_STATUS_OK .or. status == DEC_STATUS_NO_PROGRESS) .and. &
           sphere%get_facets_number() <= n_target .and. &
           abs(vol - vol_target) <= 0.30_R8P * vol_target .and. &
           sphere%get_non_manifold_edges_number() == 0_I4P)) goto 100

are_tests_passed(1) = .true.
100 continue

! ----- 2. Cube.stl → decimate to its current count (no-op).
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
n_before   = cube%get_facets_number()
vol_before = cube%get_volume()
call cube%decimate(target_facets=n_before, status=status)
vol = cube%get_volume()
are_tests_passed(2) = (status == DEC_STATUS_OK .and. &
                       cube%get_facets_number() == n_before .and. &
                       abs(vol - vol_before) <= 1.0e-12_R8P)
print '(A,I0,A,I0,A,F8.4)', 'cube no-op: target=', n_before, ' got=', cube%get_facets_number(), '  vol=', vol

! ----- 3. Cube.stl → decimate to 4 (unreachable cleanly): assert no crash.
call cube%destroy
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%decimate(target_facets=4_I4P, status=status)
vol = cube%get_volume()
! Acceptable: any status, as long as we got a non-empty mesh out and didn't crash.
are_tests_passed(3) = (status == DEC_STATUS_OK .or. status == DEC_STATUS_NO_PROGRESS) .and. &
                      cube%get_facets_number() >= 4_I4P
print '(A,I0,A,I0,A,F8.4)', 'cube to 4 (degenerate): status=', status, &
      ' got=', cube%get_facets_number(), '  vol=', vol

print '(A,3L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_decimate
