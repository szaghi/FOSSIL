!run ck_body_distance ck_body_distance
program ck_body_distance
use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: a, b
type(rigid_pose_t)       :: pose
real(R8P)                :: d, pa(3), pb(3)
integer(I4P)             :: fa, fb
logical                  :: hit

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
! b turned 45 degrees about z, its centre at (2, 0.5, 0.5): a corner points at a's face x = 1
pose = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI / 4)
pose%t = [2._R8P, 0.5_R8P, 0.5_R8P] - matmul(pose%r, [0.5_R8P, 0.5_R8P, 0.5_R8P])
call a%distance_to(b, distance=d, point_self=pa, point_other=pb, facet_self=fa, facet_other=fb, &
                   intersecting=hit, pose_other=pose)
print '(A,F8.4,A,F8.4)', 'distance ', d, ', expected 1 - sqrt(2)/2 = ', 1._R8P - sqrt(2._R8P) / 2
print '(A,3F8.4,A,I0)', 'on a ', pa, '  facet ', fa
print '(A,3F8.4,A,I0)', 'on b ', pb, '  facet ', fb
print '(A,L1)',         'touching or crossing ', hit
pose%t(1) = pose%t(1) - 0.4_R8P ! 0.4 closer: the corner enters a
call a%distance_to(b, distance=d, intersecting=hit, pose_other=pose)
print '(A,F8.4,A,L1)', '0.4 closer: distance ', d, ', touching or crossing ', hit
endprogram ck_body_distance
