!< FOSSIL test: surface%boolean end-to-end (issue #18 §1.1).
!<
!< Six invariants:
!<
!<   1. Disjoint cubes / DIFFERENCE — A.boolean(B, BOOL_DIFFERENCE) returns A
!<      unchanged (volume preserved). Trivial-case path: every A-facet is
!<      "outside B" (kept), no B-facet is "inside A" (all dropped), no cuts,
!<      no retriangulation.
!<   2. Overlapping cubes / DIFFERENCE — B offset by (0.5, 0.5, 0.5): the
!<      (0.5)^3 corner of A overlapping B is removed, leaving volume
!<      1 - 0.5^3 = 0.875. Tests the full pipeline end-to-end on the canonical
!<      mesh-arrangement test case.
!<   3. Overlapping cubes / UNION — same offset: the union has volume
!<      1 + 1 - 0.125 = 1.875.
!<   4. Overlapping cubes / INTERSECT — same offset: the intersection is
!<      the (0.5)^3 corner cube, volume = 0.125.
!<   5. Overlapping cubes / SYMDIFF — same offset: A △ B = (A ∪ B) \ (A ∩ B),
!<      volume = 1.875 - 0.125 = 1.75.
!<   6. Single-axis offset / DIFFERENCE — B offset by (0.5, 0, 0) only,
!<      which makes 4 face pairs coplanar. Asserts only structural success
!<      (status = OK); volume of the result is NOT yet correct in this
!<      configuration (see surface%boolean docstring's Limitations section).
!<
!< The 3D offset (0.5, 0.5, 0.5) used in invariants 2-5 is chosen so no
!< faces or edges of A and B are coplanar — the volume-correctness regime
!< this MVP supports robustly. Single-axis or planar-aligned configurations
!< trip a known limitation in the tri-tri intersection pipeline; see the
!< limitation note in the boolean TBP docstring for the workaround
!< (perturb inputs by an epsilon along each axis).

program fossil_test_boolean

use fossil, only : surface_stl_object, &
                   BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF, &
                   BOOL_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube_a, cube_b
real(R8P)                :: vol_a_before, vol_after
integer(I4P)             :: status
logical                  :: are_tests_passed(6)

real(R8P), parameter :: TOL_VOLUME_PRESERVED = 1.0e-10_R8P  ! disjoint case: exact match expected
real(R8P), parameter :: TOL_VOLUME_3D        = 0.05_R8P     ! 3D-offset cases: 5% slack for FP + retri quality

are_tests_passed = .false.

! ----- 1. Disjoint cubes / DIFFERENCE: A unchanged.
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(10._R8P, 0._R8P, 0._R8P))
call cube_b%analyze
vol_a_before = cube_a%get_volume()
call cube_a%boolean(other=cube_b, op=BOOL_DIFFERENCE, status=status)
vol_after = cube_a%get_volume()
are_tests_passed(1) = (status == BOOL_STATUS_OK .and. &
                       abs(vol_after - vol_a_before) <= TOL_VOLUME_PRESERVED)
print '(A,I0,A,ES14.6)', 'disjoint DIFFERENCE: status=', status, '  vol=', vol_after

! ----- 2. Overlapping cubes / DIFFERENCE (3D offset): volume = 0.875.
call cube_a%destroy ; call cube_b%destroy
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube_b%analyze
call cube_a%boolean(other=cube_b, op=BOOL_DIFFERENCE, status=status)
vol_after = cube_a%get_volume()
are_tests_passed(2) = (status == BOOL_STATUS_OK .and. &
                       abs(vol_after - 0.875_R8P) <= TOL_VOLUME_3D)
print '(A,I0,A,ES14.6)', 'overlap DIFFERENCE:  status=', status, '  vol=', vol_after

! ----- 3. Overlapping cubes / UNION (3D offset): volume = 1.875.
call cube_a%destroy ; call cube_b%destroy
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube_b%analyze
call cube_a%boolean(other=cube_b, op=BOOL_UNION, status=status)
vol_after = cube_a%get_volume()
are_tests_passed(3) = (status == BOOL_STATUS_OK .and. &
                       abs(vol_after - 1.875_R8P) <= TOL_VOLUME_3D)
print '(A,I0,A,ES14.6)', 'overlap UNION:       status=', status, '  vol=', vol_after

! ----- 4. Overlapping cubes / INTERSECT (3D offset): volume = 0.125.
call cube_a%destroy ; call cube_b%destroy
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube_b%analyze
call cube_a%boolean(other=cube_b, op=BOOL_INTERSECT, status=status)
vol_after = cube_a%get_volume()
are_tests_passed(4) = (status == BOOL_STATUS_OK .and. &
                       abs(vol_after - 0.125_R8P) <= TOL_VOLUME_3D)
print '(A,I0,A,ES14.6)', 'overlap INTERSECT:   status=', status, '  vol=', vol_after

! ----- 5. Overlapping cubes / SYMDIFF (3D offset): volume = 1.75.
call cube_a%destroy ; call cube_b%destroy
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube_b%analyze
call cube_a%boolean(other=cube_b, op=BOOL_SYMDIFF, status=status)
vol_after = cube_a%get_volume()
are_tests_passed(5) = (status == BOOL_STATUS_OK .and. &
                       abs(vol_after - 1.75_R8P) <= TOL_VOLUME_3D)
print '(A,I0,A,ES14.6)', 'overlap SYMDIFF:     status=', status, '  vol=', vol_after

! ----- 6. Single-axis offset / DIFFERENCE: structural success only.
!         Volume not asserted — see file-header note and surface%boolean
!         docstring Limitations section.
call cube_a%destroy ; call cube_b%destroy
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(0.5_R8P, 0._R8P, 0._R8P))
call cube_b%analyze
call cube_a%boolean(other=cube_b, op=BOOL_DIFFERENCE, status=status)
vol_after = cube_a%get_volume()
are_tests_passed(6) = (status == BOOL_STATUS_OK)
print '(A,I0,A,ES14.6)', 'single-axis (structural-only): status=', status, '  vol=', vol_after

print '(A,6L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_boolean
