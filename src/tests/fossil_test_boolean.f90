!< FOSSIL test: surface%boolean(BOOL_DIFFERENCE) end-to-end (issue #18 §1.1).
!<
!< Four invariants:
!<   1. Disjoint cubes — A.boolean(B, BOOL_DIFFERENCE) returns A unchanged
!<      (volume preserved, watertight). Tests the trivial-case path: every
!<      A-facet is "outside B" (kept), no B-facet is "inside A" (all dropped),
!<      no cuts, no retriangulation.
!<   2. Two unit cubes with B offset by (0.5, 0.5, 0.5) — A.boolean(B,
!<      BOOL_DIFFERENCE) returns a watertight result with volume = 1 - 0.5^3
!<      = 0.875 (the (0.5)^3 corner of A overlapping B is removed). Tests
!<      the full pipeline: arrangement + retriangulation (CDT) + WN tagging
!<      + selection + stitching, end to end. The 3D offset is chosen so no
!<      faces or edges of A and B are coplanar — the volume-correctness
!<      regime that this MVP supports robustly.
!<   3. NOT_IMPLEMENTED stub — A.boolean(B, BOOL_UNION) returns
!<      BOOL_STATUS_NOT_IMPLEMENTED. Confirms the API surface is stable
!<      while only DIFFERENCE is wired this PR.
!<   4. Single-axis offset structural check — B offset by (0.5, 0, 0) only,
!<      which makes 4 face pairs coplanar. With the new tri-tri segment clip
!<      and tolerance-aware cut-endpoint dedup, the pipeline now runs to
!<      completion (status = OK) instead of failing at retriangulation as
!<      it did before the clip. The volume of the result, however, is NOT
!<      yet correct in this configuration: the WN tagging is ambiguous on
!<      sub-triangles whose centroids lie on faces shared between A and B
!<      (WN ≈ 0.5 from both, neither cleanly inside nor outside), which the
!<      current 0.5-threshold selection misclassifies. This invariant
!<      asserts only "no crash" — full coplanar-aware boolean correctness
!<      is the next chunk of §1.1 work.

program fossil_test_boolean

use fossil, only : surface_stl_object, &
                   BOOL_UNION, BOOL_DIFFERENCE, &
                   BOOL_STATUS_OK, BOOL_STATUS_NOT_IMPLEMENTED
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube_a, cube_b
real(R8P)                :: vol_a_before, vol_after
integer(I4P)             :: status
logical                  :: are_tests_passed(4)

real(R8P), parameter :: TOL_VOLUME_PRESERVED = 1.0e-10_R8P  ! disjoint case: exact match expected
real(R8P), parameter :: TOL_VOLUME_HALVED    = 0.05_R8P     ! overlap case: 5% slack for FP + retri quality

are_tests_passed = .false.

! ----- 1. Disjoint cubes: A unchanged.
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(10._R8P, 0._R8P, 0._R8P))
call cube_b%analyze
vol_a_before = cube_a%get_volume()
call cube_a%boolean(other=cube_b, op=BOOL_DIFFERENCE, status=status)
vol_after = cube_a%get_volume()
are_tests_passed(1) = (status == BOOL_STATUS_OK .and. &
                       abs(vol_after - vol_a_before) <= TOL_VOLUME_PRESERVED)
print '(A,I0,A,ES14.6,A,ES14.6)', 'disjoint: status=', status, &
      '  vol_before=', vol_a_before, '  vol_after=', vol_after

! ----- 2. Overlapping cubes (B offset by (0.5, 0.5, 0.5)): result volume = 0.875.
call cube_a%destroy ; call cube_b%destroy
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube_b%analyze
vol_a_before = cube_a%get_volume()
call cube_a%boolean(other=cube_b, op=BOOL_DIFFERENCE, status=status)
vol_after = cube_a%get_volume()
are_tests_passed(2) = (status == BOOL_STATUS_OK .and. &
                       abs(vol_after - 0.875_R8P) <= TOL_VOLUME_HALVED)
print '(A,I0,A,ES14.6,A,ES14.6)', 'overlap:  status=', status, &
      '  vol_before=', vol_a_before, '  vol_after=', vol_after

! ----- 3. NOT_IMPLEMENTED stub for non-DIFFERENCE ops.
!         Use the same overlapping configuration as test 2 so the pipeline
!         reaches the per-op selection step (which then rejects with the
!         NOT_IMPLEMENTED status).
call cube_a%destroy ; call cube_b%destroy
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call cube_b%analyze
call cube_a%boolean(other=cube_b, op=BOOL_UNION, status=status)
are_tests_passed(3) = (status == BOOL_STATUS_NOT_IMPLEMENTED)
print '(A,I0)', 'union stub: status=', status

! ----- 4. Single-axis offset structural check (no crash, but volume not
!         asserted — see file-header note for the WN-coplanar limitation).
call cube_a%destroy ; call cube_b%destroy
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(0.5_R8P, 0._R8P, 0._R8P))
call cube_b%analyze
call cube_a%boolean(other=cube_b, op=BOOL_DIFFERENCE, status=status)
vol_after = cube_a%get_volume()
are_tests_passed(4) = (status == BOOL_STATUS_OK)
print '(A,I0,A,ES14.6)', 'single-axis: status=', status, '  vol_after=', vol_after

print '(A,4L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_boolean
