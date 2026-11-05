!< FOSSIL test: explicit + Taubin mesh smoothing (issue #18 §2.3).
!<
!< Regression variant (issue #45): invariants 1-3 on dragon.stl (6588 facets) instead of bunny, same thresholds
!< (measured on dragon: volume drift 0.22 %, area drop 7.2 %, explicit loss 78.0 vs Taubin 1.04). The full bunny
!< version is src/tests/validation/fossil_val_smoothing.f90.
!<
!< Four invariants on the public contract of `surface%smooth`:
!<
!<   1. **Taubin volume preservation on dragon**: |Δvolume/volume| ≤
!<      5% after 5 Taubin pairs. The dragon has real high-frequency
!<      surface variation that Taubin is designed to smooth without
!<      shrinking. This is the headline guarantee.
!<   2. **Taubin reduces surface area on dragon**: Taubin smoothing
!<      attenuates high-frequency surface "wiggle", which lowers the
!<      total surface area. Assert post-smooth area < pre-smooth area
!<      by at least 0.1% (a small but nonzero drop confirms the
!<      smoother is doing visible work).
!<   3. **Explicit shrinks the dragon more than Taubin**: with the
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
                   SMOOTH_STATUS_OK, DEV_STATUS_OK
use penf, only : I4P, R8P

implicit none

type(surface_stl_object) :: dragon_taubin, dragon_explicit, cube, dragon0
real(R8P) :: vol_before, vol_taubin_after, vol_explicit_after
real(R8P) :: area_before, area_taubin_after
integer(I4P) :: status
real(R8P), parameter :: TAUBIN_VOL_TOL_REL  = 0.05_R8P  ! ≤ 5% volume drift
real(R8P), parameter :: AREA_DROP_MIN_REL   = 0.001_R8P ! ≥ 0.1% area drop
logical :: are_tests_passed(5)
real(R8P) :: lt, ut, le, ue                    !< Hausdorff bounds, Taubin and explicit (issue #60).
integer(I4P) :: st, se
real(R8P), parameter :: H_TAUBIN_MAX = 0.45_R8P !< Bound of the Taubin dragon's Hausdorff distance (0.333 when written).

are_tests_passed = .false.

! ---- Inv 1 & 2: Taubin on dragon.
call dragon_taubin%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
dragon0 = dragon_taubin
vol_before  = dragon_taubin%get_volume()
area_before = dragon_taubin%get_area()
print '(A,F12.6,A,F12.6)', 'dragon before: vol=', vol_before, '  area=', area_before

call dragon_taubin%smooth(method=SMOOTH_METHOD_TAUBIN, iterations=5_I4P, status=status)
vol_taubin_after  = dragon_taubin%get_volume()
area_taubin_after = dragon_taubin%get_area()
print '(A,F12.6,A,F12.6,A,I0)', 'dragon Taubin: vol=', vol_taubin_after, &
      '  area=', area_taubin_after, '  status=', status

print '(A,F8.4,A)', 'inv 1 (Taubin vol): rel change = ', &
      abs(vol_taubin_after - vol_before) / abs(vol_before), ' (target ≤ 0.05)'
are_tests_passed(1) = (status == SMOOTH_STATUS_OK) .and. &
                      (abs(vol_taubin_after - vol_before) <= TAUBIN_VOL_TOL_REL * abs(vol_before))

print '(A,F8.4,A,F8.4)', 'inv 2 (area drop): rel drop = ', &
      (area_before - area_taubin_after) / area_before, '  target ≥ ', AREA_DROP_MIN_REL
are_tests_passed(2) = ((area_before - area_taubin_after) / area_before >= AREA_DROP_MIN_REL)

! ---- Inv 3: Explicit shrinks more than Taubin.
call dragon_explicit%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call dragon_explicit%smooth(method=SMOOTH_METHOD_EXPLICIT, lambda=0.5_R8P, &
                            iterations=5_I4P, status=status)
vol_explicit_after = dragon_explicit%get_volume()
print '(A,F12.6,A,F12.6,A,F12.6)', 'inv 3 (explicit vs Taubin): explicit_vol=', vol_explicit_after, &
      '  Taubin_vol=', vol_taubin_after, '  explicit_loss=', (vol_before - vol_explicit_after)
are_tests_passed(3) = (status == SMOOTH_STATUS_OK) .and. &
                      ((vol_before - vol_explicit_after) > (vol_before - vol_taubin_after))

! ---- Inv 4: Cube — status OK only, no geometry assertion.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%smooth(method=SMOOTH_METHOD_TAUBIN, iterations=3_I4P, status=status)
print '(A,I0)', 'inv 4 (cube status): status=', status
are_tests_passed(4) = (status == SMOOTH_STATUS_OK)

! ---- Inv 5 (issue #60): how far smoothing moved the dragon. Taubin's lambda|mu steps keep the surface in place, the
!      explicit Laplacian shrinks it: the Hausdorff distance to the input, certified to the default tolerance, is
!      smaller for Taubin (0.307-0.333 against 0.565-0.589 of the 26.7 diagonal when written) and below 0.45.
call dragon_taubin%hausdorff(dragon0, lower=lt, upper=ut, status=st)
call dragon_explicit%hausdorff(dragon0, lower=le, upper=ue, status=se)
print '(A,2ES11.3,A,2ES11.3)', 'inv 5 (Hausdorff to the input): Taubin ', lt, ut, '  explicit ', le, ue
are_tests_passed(5) = st == DEV_STATUS_OK .and. se == DEV_STATUS_OK .and. ut < le .and. ut <= H_TAUBIN_MAX

print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. (all(are_tests_passed))) error stop 1
endprogram fossil_test_smoothing
