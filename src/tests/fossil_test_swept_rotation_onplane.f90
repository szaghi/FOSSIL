!< FOSSIL test: exact swept geometry under rotation when the wall starts on grid lines / planes (issues #84, #86).
!<
!< A wall side that starts on a grid line (2-D) or a body face that starts on a grid plane (3-D) is sliced nearly
!< along itself for a while after s = 0: the crossings of the moved copy carry round-off amplified by 1 / (angle
!< between the wall and the grid line or plane), so on those pieces two successive Gauss values stop agreeing at a
!< floor above 64 eps of the face measure times the piece length, while the result is exact to round-off (the
!< rigid-field residual stays ~1e-16). Before the fix `sw%exact` was false in both cases below. Each case checks
!< `sw%exact`, the residual (<= 1e-14 in every closed cell) and the projected walls summing to zero:
!< - 2-D: the square [-0.5, 0.5]^2 (its sides on the grid lines x, y = +-0.5 of a 4 x 4 grid of spacing 0.5) turning
!<   0.1 rad about the origin while moving by (0.05, 0.02): the corners (0.5, 0.5) and (-0.5, -0.5) also leave
!<   their lines tangentially (zero zones of width ~1e-6 at s = 0);
!< - 3-D: the unit cube centred at the origin turning 0.35 rad about (0.3, 0.2, 1) while moving by (0.15, 0.05, 0)
!<   (issue #84's motion), grid x, y = -1, -0.5, 0, 0.5, 1 (the faces y = +-0.5 start on grid planes), z = +-0.75.
program fossil_test_swept_rotation_onplane

use fossil,                   only : surface_stl_object, curve_object, rigid_pose_t, pose_identity, pose_from_axis_angle,&
                                     pose_2d, swept_cells_2d_t, swept_cells_2d, swept_cells_3d_t, swept_cells_3d, &
                                     STATUS_OK
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(curve_object)        :: sq
type(swept_cells_2d_t)    :: s2
type(swept_cells_3d_t)    :: s3
type(rigid_pose_t)        :: p1
real(R8P)                 :: e(5), r, w(3)
integer(I4P)              :: i, st
logical                   :: ok(2)

! ---- 2-D: a square on the grid lines
call sq%adopt_loops(u=[-0.5_R8P, 0.5_R8P, 0.5_R8P, -0.5_R8P], v=[-0.5_R8P, -0.5_R8P, 0.5_R8P, 0.5_R8P], &
                    loop_sizes=[4_I4P])
e = [(-1._R8P + 0.5_R8P * i, i = 0, 4)]
call swept_cells_2d(sq, pose_2d(angle=0._R8P), pose_2d(angle=0.1_R8P, translation=[0.05_R8P, 0.02_R8P]), e, e, s2, &
                    status=st)
r = maxval(abs(s2%residual))
w(1:2) = [sum(s2%wall(1, :, :)), sum(s2%wall(2, :, :))]
ok(1) = st == STATUS_OK .and. s2%exact .and. r <= 1e-14_R8P .and. all(abs(w(1:2)) <= 1e-14_R8P)
print '(A,L1,A,ES9.2,A,2ES9.1,A,I0,A,L1)', '2-D square on grid lines: exact ', s2%exact, ', residual ', r, &
      ', walls ', w(1:2), ', evaluations/face ', s2%evaluations, ', ok ', ok(1)

! ---- 3-D: the cube with faces on grid planes, tilted axis (issue #84)
call bench_bodies(b)
call to_surface(b(1)%soup, s) ! unit cube [0, 1]^3
call s%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
call s%analyze
p1 = pose_from_axis_angle(axis=[0.3_R8P, 0.2_R8P, 1._R8P], angle=0.35_R8P, translation=[0.15_R8P, 0.05_R8P, 0._R8P])
call swept_cells_3d(s, pose_identity(), p1, e, e, [-0.75_R8P, 0.75_R8P], s3, status=st)
r = maxval(abs(s3%residual), mask=s3%alpha0 >= 0._R8P .and. s3%alpha1 >= 0._R8P)
w = [sum(s3%wall(1, :, :, :)), sum(s3%wall(2, :, :, :)), sum(s3%wall(3, :, :, :))]
ok(2) = st == STATUS_OK .and. s3%exact .and. r <= 1e-14_R8P .and. all(abs(w) <= 1e-14_R8P)
print '(A,L1,A,ES9.2,A,3ES9.1,A,I0,A,L1)', '3-D cube on grid planes, tilted axis: exact ', s3%exact, ', residual ', r, &
      ', walls ', w, ', evaluations/plane ', s3%evaluations, ', ok ', ok(2)
print '(A,L1)', 'Are all tests passed? ', all(ok)
endprogram fossil_test_swept_rotation_onplane
