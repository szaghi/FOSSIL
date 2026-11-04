program ck_pose
use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle
use penf,   only : R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: cube
type(rigid_pose_t)       :: pose
type(vector_R8P)         :: p(1)
real(R8P)                :: d(1)

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
! 45 degrees about the z axis through the cube centre: x_world = R (x_body - c) + c
pose = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI / 4)
pose%t = [0.5_R8P, 0.5_R8P, 0._R8P] - matmul(pose%r, [0.5_R8P, 0.5_R8P, 0._R8P])
p = vector_R8P(1.2_R8P, 0.5_R8P, 0.5_R8P) ! world frame
call cube%distance_many(points=p, distances=d, is_signed=.true., is_square_root=.true., pose=pose)
print '(A,F8.4)', 'signed distance, rotated cube ', d(1)
print '(A,L1)',   'inside, rotated cube          ', cube%is_point_inside(point=p(1), pose=pose)
print '(A,F8.4)', 'signed distance, cube at rest ', cube%distance(point=p(1), is_signed=.true., is_square_root=.true.)
endprogram ck_pose
