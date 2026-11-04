!< FOSSIL test: optional self-intersection checks of the curve validity and of the swept-geometry inputs (issue #44).
program fossil_test_self_intersection_checks
!< FOSSIL test: optional self-intersection checks of the curve validity and of the swept-geometry inputs (issue #44).
!<
!< 1. A counter-clockwise square with a knot on its top edge (the segment (2, 5)-(2, 3) crosses the edge y = 4 at
!<    (2, 4)) is valid by default, as before, and CURVE_SELF_INTERSECTING in loop 1 with
!<    `check_self_intersections`; the plain square stays CURVE_VALID with the check.
!< 2. swept_cells_2d(check_input=.true.) flags the knotted curve, not the square; without check_input nothing is
!<    flagged (default unchanged).
!< 3. swept_cells_3d(check_input=.true.) flags two overlapping cubes merged into one surface, not a single cube.
!<    The offset is generic: at (0.5, 0.5, 0.5) the cubes share vertex coordinates and find_self_intersections,
!<    which skips facet pairs sharing a vertex, finds no pair although the faces cross (a known limitation).

use fossil,             only : surface_stl_object, curve_object, CURVE_VALID, CURVE_SELF_INTERSECTING, rigid_pose_t, &
                               rigid_pose_2d_t, pose_identity, pose_2d, swept_cells_2d_t, swept_cells_2d, &
                               swept_cells_3d_t, swept_cells_3d
use penf,               only : I4P, R8P
use vecfor,             only : ex_R8P, ey_R8P, ez_R8P

implicit none
type(curve_object)       :: knot, square
type(surface_stl_object) :: cube, other, cubes
type(swept_cells_2d_t)   :: sw2
type(swept_cells_3d_t)   :: sw3
real(R8P)                :: e(9)
integer(I4P)             :: loop, i
logical                  :: ok(3)

! 1. validity
call knot%adopt_loops(u=[0._R8P, 4._R8P, 4._R8P, 1._R8P, 2._R8P, 2._R8P, 0._R8P], &
                      v=[0._R8P, 0._R8P, 4._R8P, 4._R8P, 5._R8P, 3._R8P, 4._R8P], loop_sizes=[7_I4P])
call square%adopt_loops(u=[0._R8P, 4._R8P, 4._R8P, 0._R8P], v=[0._R8P, 0._R8P, 4._R8P, 4._R8P], loop_sizes=[4_I4P])
ok(1) = knot%check_validity() == CURVE_VALID .and. &
        knot%check_validity(loop=loop, check_self_intersections=.true.) == CURVE_SELF_INTERSECTING .and. loop == 1 .and. &
        square%check_validity(check_self_intersections=.true.) == CURVE_VALID
print '(A,L1)', 'knot valid by default, self-intersecting when checked; square valid: ', ok(1)

! 2. swept_cells_2d input check
e = [(-0.5_R8P + 0.625_R8P * i, i = 0, 8)]
call swept_cells_2d(knot, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P, translation=[0.1_R8P, 0._R8P]), e, e, sw2, &
                    check_input=.true.)
ok(2) = sw2%self_intersecting
call swept_cells_2d(square, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P, translation=[0.1_R8P, 0._R8P]), e, e, sw2, &
                    check_input=.true.)
ok(2) = ok(2) .and. .not. sw2%self_intersecting
call swept_cells_2d(knot, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P, translation=[0.1_R8P, 0._R8P]), e, e, sw2)
ok(2) = ok(2) .and. .not. sw2%self_intersecting
print '(A,L1)', 'swept_cells_2d: knot flagged, square not, no flag without check_input: ', ok(2)

! 3. swept_cells_3d input check
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
other = cube
call other%translate(delta=0.37_R8P * ex_R8P + 0.41_R8P * ey_R8P + 0.29_R8P * ez_R8P) ! generic: no shared vertices
cubes = cube
call cubes%merge_solids(other)
call cubes%analyze ! merge_solids does not rebuild the tree the broad phase uses
e(1:5) = [(-0.5_R8P + 0.5_R8P * i, i = 0, 4)]
call swept_cells_3d(cubes, pose_identity(), pose_identity(), e(1:5), e(1:5), e(1:5), sw3, check_input=.true.)
ok(3) = sw3%self_intersecting
call swept_cells_3d(cube, pose_identity(), pose_identity(), e(1:5), e(1:5), e(1:5), sw3, check_input=.true.)
ok(3) = ok(3) .and. .not. sw3%self_intersecting
print '(A,L1)', 'swept_cells_3d: overlapping cubes flagged, one cube not: ', ok(3)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1
endprogram fossil_test_self_intersection_checks
