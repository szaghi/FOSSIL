!< FOSSIL test: alpha wrapping, `surface%alpha_wrap` (issue #18 §1.6, issue #50).
!<
!< REGRESSION variant (issue #45): the validation copy src/tests/validation/fossil_val_alpha_wrap.f90 runs the same
!< wraps at a finer alpha, on a marching-cubes sphere and on dragon.stl, and checks them for self-crossings
!< (find_self_intersections, 2-5 s each in the debug build) and for their distance from the input.
!<
!<   1. The cube: status OK, watertight, manifold, volume within 0.5% of the offset cube (the Minkowski sum of the
!<      cube and a ball of radius offset).
!<   2. The cube minus one facet, hole wider than the ball (alpha = 0.4): one closed shell lining both sides of the
!<      walls (the offset-surface barrier keeps the two layers apart), OK, manifold, volume of a hollow shell.
!<   3. The same, hole narrower than the ball (alpha = 0.6): the hole is capped, volume just below the offset cube's.
!<   4. One triangle (zero thickness): a closed shell around both sides, volume 2 offset x area plus the rim.
!<   5. Bad input: alpha larger than the bounding-box diagonal is DEGENERATE with an empty wrap; an empty surface and
!<      alpha <= 0 are BAD_INPUT.

program fossil_test_alpha_wrap

use fossil, only : surface_stl_object, AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT, AWRAP_STATUS_DEGENERATE
use fossil_facet_object, only : facet_object
use fossil_utils, only : PI
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube, cube_holed, wrapped_tbp, plate, empty
type(facet_object), allocatable :: plate_facets(:)
type(facet_object), allocatable, target :: holed_facets(:)
type(facet_object), pointer :: facets(:)
real(R8P) :: wrap_volume, expected_volume, offset
integer(I4P) :: status, i
logical :: are_tests_passed(5)

are_tests_passed = .false.

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
facets => cube%facets_ref()
allocate(holed_facets(size(facets, kind=I4P) - 1_I4P))
do i = 1_I4P, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1_I4P)
enddo
call cube_holed%adopt_facets(facets=holed_facets)

! ---- Invariant 1: capstone TBP on the cube (issue #50): status OK, closed, manifold, and its
!      volume is the cube grown by the offset (Minkowski sum with a ball: 1 + 6e + 3 pi e^2 + 4/3 pi e^3).
offset = 0.033_R8P
call cube%alpha_wrap(alpha=0.3_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
expected_volume = 1._R8P + 6._R8P * offset + 3._R8P * PI * offset**2 + 4._R8P / 3._R8P * PI * offset**3
print '(A,I0,A,I0,2(A,L1),A,F8.5,A,F8.5)', &
      'inv 1 (cube TBP): status=', status, ' n_facets=', wrapped_tbp%get_facets_number(), &
      ' watertight=', wrapped_tbp%is_watertight(), ' manifold=', wrapped_tbp%is_manifold(), &
      ' volume=', wrap_volume, ' offset cube=', expected_volume
are_tests_passed(1) = (status == AWRAP_STATUS_OK) .and. wrapped_tbp%is_watertight() .and. &
                       wrapped_tbp%is_manifold() .and. &
                       (abs(wrap_volume - expected_volume) <= 0.005_R8P * expected_volume)

! ---- Invariant 2: the holed cube (a hole of inradius 0.29) at alpha = 0.4: the ball (diameter alpha, kept offset
!      from the input) gets through the hole, so the wrap lines both sides of every wall: one closed shell, two offset
!      layers that the barrier keeps apart (they collapsed onto one wall before issue #50).
call cube_holed%alpha_wrap(alpha=0.4_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
print '(A,I0,A,I0,2(A,L1),A,F8.5)', &
      'inv 2 (holed cube, hole open): status=', status, ' n_facets=', wrapped_tbp%get_facets_number(), &
      ' watertight=', wrapped_tbp%is_watertight(), ' manifold=', wrapped_tbp%is_manifold(), &
      ' volume=', wrap_volume
! hollow: at least the 2-offset layer over the 5.5 units of wall area, and far below the solid's 1.2 (the ball, radius
! 0.2, rounds the inner corners, so the inner layer leaves the walls there)
are_tests_passed(2) = (status == AWRAP_STATUS_OK) .and. wrapped_tbp%is_watertight() .and. &
                       wrapped_tbp%is_manifold() .and. &
                       (wrap_volume >= 11._R8P * offset) .and. (wrap_volume <= 0.7_R8P)

! ---- Invariant 3: the same holed cube at alpha = 0.6: the ball does not fit through the hole, the wrap closes it
!      with a cap that dips into it, so the volume is just below the offset cube's.
offset = 0.02_R8P
call cube_holed%alpha_wrap(alpha=0.6_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
expected_volume = 1._R8P + 6._R8P * offset + 3._R8P * PI * offset**2 + 4._R8P / 3._R8P * PI * offset**3
print '(A,I0,A,I0,2(A,L1),A,F8.5,A,F8.5)', &
      'inv 3 (holed cube, hole capped): status=', status, ' n_facets=', wrapped_tbp%get_facets_number(), &
      ' watertight=', wrapped_tbp%is_watertight(), ' manifold=', wrapped_tbp%is_manifold(), &
      ' volume=', wrap_volume, ' offset cube=', expected_volume
are_tests_passed(3) = (status == AWRAP_STATUS_OK) .and. wrapped_tbp%is_watertight() .and. &
                       wrapped_tbp%is_manifold() .and. &
                       (wrap_volume < expected_volume) .and. (wrap_volume > 0.98_R8P * expected_volume)

! ---- Invariant 4: one triangle (a sheet of zero thickness, area 1/2): a closed shell around both of its sides,
!      volume 2 offset x area plus the rim.
allocate(plate_facets(1))
plate_facets(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
plate_facets(1)%vertex(2) = vector_R8P(1._R8P, 0._R8P, 0._R8P)
plate_facets(1)%vertex(3) = vector_R8P(0._R8P, 1._R8P, 0._R8P)
call plate_facets(1)%compute_metrix
call plate%adopt_facets(facets=plate_facets)
call plate%alpha_wrap(alpha=0.3_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
print '(A,I0,A,I0,2(A,L1),A,F8.5)', &
      'inv 4 (one triangle): status=', status, ' n_facets=', wrapped_tbp%get_facets_number(), &
      ' watertight=', wrapped_tbp%is_watertight(), ' manifold=', wrapped_tbp%is_manifold(), &
      ' volume=', wrap_volume
! 2 e A = 0.02; the rim (perimeter 3.41, half a cylinder of radius e) adds at most 0.0022
are_tests_passed(4) = (status == AWRAP_STATUS_OK) .and. wrapped_tbp%is_watertight() .and. &
                       wrapped_tbp%is_manifold() .and. &
                       (wrap_volume >= 0.02_R8P) .and. (wrap_volume <= 0.0225_R8P)

! ---- Invariant 5: bad input.
call cube%alpha_wrap(alpha=100._R8P, wrapped=wrapped_tbp, status=status)
are_tests_passed(5) = (status == AWRAP_STATUS_DEGENERATE) .and. (wrapped_tbp%get_facets_number() == 0_I4P)
print '(A,I0,A,I0)', 'inv 5 (alpha > bbox diagonal): status=', status, ' n_facets=', wrapped_tbp%get_facets_number()
call empty%alpha_wrap(wrapped=wrapped_tbp, status=status)
are_tests_passed(5) = are_tests_passed(5) .and. (status == AWRAP_STATUS_BAD_INPUT)
print '(A,I0)', 'inv 5 (empty surface): status=', status
call cube%alpha_wrap(alpha=-1._R8P, wrapped=wrapped_tbp, status=status)
are_tests_passed(5) = are_tests_passed(5) .and. (status == AWRAP_STATUS_BAD_INPUT)
print '(A,I0)', 'inv 5 (alpha <= 0): status=', status

print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. (all(are_tests_passed))) error stop 1
endprogram fossil_test_alpha_wrap
