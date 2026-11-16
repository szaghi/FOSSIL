!run ck_swept_bodies ck_swept_bodies
program ck_swept_bodies
use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, swept_cells_3d_t, swept_cells_3d, &
                   swept_cells_3d_bodies, curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, &
                   swept_cells_2d_bodies, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: cube(2)
type(rigid_pose_t)       :: pose0(2), pose1(2)
type(swept_cells_3d_t)   :: sw, alone
type(curve_object)       :: disc
type(swept_cells_2d_t)   :: sw2
real(R8P)                :: x(17), u(32), v(32)
integer(I4P)             :: i, status

call cube(1)%load_from_file(file_name='cube-01.stl', guess_format=.true.)
cube(2) = cube(1)
x = [(-2._R8P + 0.25_R8P * i, i = 0, 16)]                      ! 16^3 cells over [-2, 2]^3
pose0(1) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, translation=[-1.5_R8P, -0.4_R8P, 0._R8P])
pose1(1) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0.2_R8P, translation=[-1.3_R8P, -0.4_R8P, 0._R8P])
pose0(2) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=0._R8P, translation=[0.3_R8P, 0.4_R8P, 0._R8P])
pose1(2) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=-0.1_R8P, translation=[0.1_R8P, 0.4_R8P, 0._R8P])
call swept_cells_3d_bodies(cube, pose0, pose1, x, x, x, sw, exact=.false., status=status) ! each cube on its own motion
if (status /= STATUS_OK) error stop 'swept_cells_3d_bodies failed'
print '(A,I0,A,ES10.2,A,ES10.2)', 'cells with a swept volume ', count(sw%swept /= 0._R8P), ', sum of swept ', &
                                  sum(sw%swept), ', largest residual ', maxval(abs(sw%residual))
call swept_cells_3d(cube(1), pose0(1), pose1(1), x, x, x, alone, exact=.false.)
print '(A,ES10.2)', 'cube 1 alone, its half of the block: largest face difference ', &
                    max(maxval(abs(sw%face_x(1:8, :, :) - alone%face_x(1:8, :, :))), &
                        maxval(abs(sw%face_y(1:8, :, :) - alone%face_y(1:8, :, :))))
call swept_cells_3d_bodies(cube, pose0, pose1, x, x, x, sw, status=status)
print '(A,I0)', 'exact path asked: status ', status
do i = 1, 32
   u(i) = 0.2_R8P * cos(0.19634954084936207_R8P * (i - 1)) ; v(i) = 0.2_R8P * sin(0.19634954084936207_R8P * (i - 1))
enddo
call disc%adopt_loops(u=u, v=v, loop_sizes=[32_I4P])            ! 2-D: two discs passing each other
call swept_cells_2d_bodies([disc, disc], [pose_2d(0._R8P, [-0.3_R8P, 0.22_R8P]), pose_2d(0._R8P, [0.3_R8P, -0.22_R8P])], &
                           [pose_2d(0.3_R8P, [-0.1_R8P, 0.22_R8P]), pose_2d(-0.3_R8P, [0.1_R8P, -0.22_R8P])], &
                           x(5:13), x(5:13), sw2, exact=.false., status=status)
print '(A,I0,A,ES10.2,A,ES10.2)', '2-D: status ', status, ', sum of swept ', sum(sw2%swept), ', error estimate ', &
                                  sw2%error_estimate
endprogram ck_swept_bodies
