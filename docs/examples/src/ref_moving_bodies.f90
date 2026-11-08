!run ref_moving_bodies ref_moving_bodies
!< Reference page "Moving bodies": pose algebra and the screw motion, body-frame queries against a moved copy and a
!< mesh moved and rebuilt, and the swept geometry of one step of the unit cube: exact translation, the Gauss path as
!< a negative control, and an exact rotation.
program ref_moving_bodies
!region poses
use fossil, only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, pose_compose, pose_inverse, &
                   pose_interpolate, to_world, to_body, swept_cells_3d_t, swept_cells_3d
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(rigid_pose_t)       :: p0, p1, q, half
type(surface_stl_object) :: cube, moved, rebuilt
type(swept_cells_3d_t)   :: sw
type(vector_R8P)         :: points(27)
real(R8P)                :: d_posed(27), d_moved(27), d_rebuilt(27), x(7), b(3)
integer(I4P)             :: i, j, k, l

p0 = pose_identity()
p1 = pose_from_axis_angle(axis=[1._R8P, 1._R8P, 1._R8P], angle=PI / 6._R8P, translation=[0.1_R8P, -0.05_R8P, 0.02_R8P])
q = pose_compose(pose_inverse(p1), p1)                                          ! p1^-1 o p1 = identity
print '(A,ES9.2)', 'p1^-1 o p1 - identity: ', max(maxval(abs(q%r - p0%r)), maxval(abs(q%t)))
b = to_body(p1, to_world(p1, [0.3_R8P, 0.2_R8P, 0.1_R8P]))
print '(A,ES9.2)', 'to_body(to_world(x)) - x: ', maxval(abs(b - [0.3_R8P, 0.2_R8P, 0.1_R8P]))
half = pose_interpolate(p0, p1, 0.5_R8P)                                       ! the screw motion at mid-step
print '(A,ES9.2,A,ES9.2)', 'R(1/2)^2 - R1: ', maxval(abs(matmul(half%r, half%r) - p1%r)), &  ! angle linear in s
                           ', t(1/2) - t1/2: ', maxval(abs(half%t - 0.5_R8P * p1%t))           ! translation linear in s
q = pose_interpolate(p0, p1, 1._R8P)
print '(A,L1)', 'interpolate(p0, p1, 1) = p1 bit for bit: ', all(q%r == p1%r) .and. all(q%t == p1%t)
!endregion poses
!region queries
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
i = 0
do k = 0, 2 ; do j = 0, 2 ; do l = 0, 2                                     ! 27 world points across the moved cube
   i = i + 1
   points(i) = vector_R8P(-0.1_R8P + 0.6_R8P * l, -0.1_R8P + 0.6_R8P * j, -0.1_R8P + 0.6_R8P * k)
enddo ; enddo ; enddo
call cube%distance_many(points=points, distances=d_posed, is_signed=.true., is_square_root=.true., pose=p1)
moved = cube
call moved%move(p1)                                                              ! refit, not rebuilt
call moved%distance_many(points=points, distances=d_moved, is_signed=.true., is_square_root=.true.)
rebuilt = moved
call rebuilt%analyze                                                             ! the tree rebuilt from scratch
call rebuilt%distance_many(points=points, distances=d_rebuilt, is_signed=.true., is_square_root=.true.)
print '(A,ES9.2,A,I0)', 'posed - moved: ', maxval(abs(d_posed - d_moved)), ', inside ', count(d_posed < 0._R8P)
print '(A,L1)', 'moved = rebuilt, bit for bit: ', all(d_moved == d_rebuilt)
!endregion queries
!region swept
x = [(-0.22_R8P + 0.25_R8P * i, i = 0, 6)]                                     ! 6^3 cells of side 1/4
call swept_cells_3d(cube, p0, pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, &
                    translation=[0.1_R8P, 0.05_R8P, 0._R8P]), x, x, x, sw)
print '(A,L1,A,I0,A,ES9.2,A,ES9.2)', 'translation: exact ', sw%exact, ', evaluations ', sw%evaluations, &
                                     ', residual ', maxval(abs(sw%residual)), ', sum of swept ', sum(sw%swept)
call swept_cells_3d(cube, p0, pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, &
                    translation=[0.1_R8P, 0.05_R8P, 0._R8P]), x, x, x, sw, exact=.false.)
print '(A,L1,A,I0,A,ES9.2,A,ES9.2)', '  Gauss:     exact ', sw%exact, ', nodes ', sw%evaluations, &
                                     ', residual ', maxval(abs(sw%residual)), ', estimate ', sw%error_estimate
call swept_cells_3d(cube, p0, pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0.2_R8P, &
                    translation=[0.1_R8P, 0.05_R8P, 0._R8P]), x, x, x, sw)
print '(A,L1,A,I0,A,ES9.2)', 'rotation:    exact ', sw%exact, ', evaluations ', sw%evaluations, &
                             ', residual ', maxval(abs(sw%residual))
print '(A,I0,A,I0,A,ES10.3)', '  cleared ', count(sw%cleared), ', covered ', count(sw%covered), &
                              ', wall sum ', maxval(abs(sum(sum(sum(sw%wall, 4), 3), 2)))
!endregion swept
endprogram ref_moving_bodies
