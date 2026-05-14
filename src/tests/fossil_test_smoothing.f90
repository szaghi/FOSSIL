!< FOSSIL test: explicit + Taubin mesh smoothing (issue #18 §2.3).
!<
!< Four invariants on the public contract of `surface%smooth`:
!<
!<   1. **Taubin volume preservation on bunny**: |Δvolume/volume| ≤
!<      5% after 5 Taubin pairs. The bunny has real high-frequency
!<      surface variation that Taubin is designed to smooth without
!<      shrinking. This is the headline guarantee.
!<   2. **Taubin reduces surface area on bunny**: Taubin smoothing
!<      attenuates high-frequency surface "wiggle", which lowers the
!<      total surface area. Assert post-smooth area < pre-smooth area
!<      by at least 0.1% (a small but nonzero drop confirms the
!<      smoother is doing visible work).
!<   3. **Explicit shrinks the bunny more than Taubin**: with the
!<      same lambda and iteration count, explicit Laplacian (no
!<      counter-shrinking) should lose more volume than Taubin. This
!<      is the "wrong tool for the job" comparison: catches a
!<      regression where explicit accidentally stops shrinking or
!<      Taubin starts shrinking.
!<   4. **Status OK on cube degenerate input**: the cube has only 8
!<      corner vertices, no flat-region interior vertices, so
!<      smoothing collapses it (corner pulls toward face centroid).
!<      We don't assert anything about the resulting geometry — only
!<      that the routine returns SMOOTH_STATUS_OK without crashing.
!<      This guards against null-vertex / division-by-zero issues
!<      on minimal inputs.

program fossil_test_smoothing

use fossil, only : surface_stl_object, &
                   SMOOTH_METHOD_EXPLICIT, SMOOTH_METHOD_TAUBIN, &
                   SMOOTH_STATUS_OK
use penf, only : I4P, R8P

implicit none

type(surface_stl_object) :: bunny_taubin, bunny_explicit, cube
real(R8P) :: vol_before, vol_taubin_after, vol_explicit_after
real(R8P) :: area_before, area_taubin_after
integer(I4P) :: status
real(R8P), parameter :: TAUBIN_VOL_TOL_REL  = 0.05_R8P  ! ≤ 5% volume drift
real(R8P), parameter :: AREA_DROP_MIN_REL   = 0.001_R8P ! ≥ 0.1% area drop
logical :: are_tests_passed(4)

are_tests_passed = .false.

! ---- Inv 1 & 2: Taubin on bunny.
call bunny_taubin%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
vol_before  = bunny_taubin%get_volume()
area_before = bunny_taubin%get_area()
print '(A,F12.6,A,F12.6)', 'bunny before: vol=', vol_before, '  area=', area_before

call bunny_taubin%smooth(method=SMOOTH_METHOD_TAUBIN, iterations=5_I4P, status=status)
vol_taubin_after  = bunny_taubin%get_volume()
area_taubin_after = bunny_taubin%get_area()
print '(A,F12.6,A,F12.6,A,I0)', 'bunny Taubin: vol=', vol_taubin_after, &
      '  area=', area_taubin_after, '  status=', status

print '(A,F8.4,A)', 'inv 1 (Taubin vol): rel change = ', &
      abs(vol_taubin_after - vol_before) / abs(vol_before), ' (target ≤ 0.05)'
are_tests_passed(1) = (status == SMOOTH_STATUS_OK) .and. &
                      (abs(vol_taubin_after - vol_before) <= TAUBIN_VOL_TOL_REL * abs(vol_before))

print '(A,F8.4,A,F8.4)', 'inv 2 (area drop): rel drop = ', &
      (area_before - area_taubin_after) / area_before, '  target ≥ ', AREA_DROP_MIN_REL
are_tests_passed(2) = ((area_before - area_taubin_after) / area_before >= AREA_DROP_MIN_REL)

! ---- Inv 3: Explicit shrinks more than Taubin.
call bunny_explicit%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call bunny_explicit%smooth(method=SMOOTH_METHOD_EXPLICIT, lambda=0.5_R8P, &
                            iterations=5_I4P, status=status)
vol_explicit_after = bunny_explicit%get_volume()
print '(A,F12.6,A,F12.6,A,F12.6)', 'inv 3 (explicit vs Taubin): explicit_vol=', vol_explicit_after, &
      '  Taubin_vol=', vol_taubin_after, '  explicit_loss=', (vol_before - vol_explicit_after)
are_tests_passed(3) = (status == SMOOTH_STATUS_OK) .and. &
                      ((vol_before - vol_explicit_after) > (vol_before - vol_taubin_after))

! ---- Inv 4: Cube — status OK only, no geometry assertion.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%smooth(method=SMOOTH_METHOD_TAUBIN, iterations=3_I4P, status=status)
print '(A,I0)', 'inv 4 (cube status): status=', status
are_tests_passed(4) = (status == SMOOTH_STATUS_OK)

print '(A,4L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
endprogram fossil_test_smoothing
