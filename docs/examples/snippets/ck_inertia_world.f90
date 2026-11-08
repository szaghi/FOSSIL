program ck_inertia_world
use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, inertia_to_world, principal_frame
use penf,   only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: block
type(rigid_pose_t)       :: pose, frame
real(R8P)                :: center(3), inertia(3,3), world(3,3), direct(3,3), moments(3), axes(3,3), diag(3,3)
integer(I4P)             :: i

call block%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call block%resize(factor=0.1_R8P * ex_R8P + 0.2_R8P * ey_R8P + 0.3_R8P * ez_R8P) ; call block%analyze
call block%mass_properties(inertia=inertia, density=7850._R8P)            ! body frame, about the centre
pose = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI / 6._R8P)  ! the body turns 30 degrees about z
world = inertia_to_world(inertia, pose)                                    ! R I R^T, no new integration
call block%rotate(axis=ez_R8P, angle=PI / 6._R8P) ; call block%analyze     ! check: integrate the turned body
call block%mass_properties(center=center, inertia=direct, principal_moments=moments, principal_axes=axes, &
                           density=7850._R8P)
print '(A,3F10.6)', 'body frame:  Ixx, Iyy, Izz ', (inertia(i,i), i = 1, 3)
print '(A,3F10.6)', 'world frame: Ixx, Iyy, Ixy ', world(1,1), world(2,2), world(1,2)
print '(A,ES10.3)', 'largest difference from the turned body ', maxval(abs(world - direct))
frame = principal_frame(center, axes)                                      ! principal axes as a pose
diag = 0._R8P ; diag(1,1) = moments(1) ; diag(2,2) = moments(2) ; diag(3,3) = moments(3)
print '(A,F9.6,A,2F9.4)', 'largest moment ', moments(3), ' about the axis (x, y) ', frame%r(1:2,3)
print '(A,ES10.3)', 'R diag(moments) R^T vs world ', maxval(abs(inertia_to_world(diag, frame) - world))
endprogram ck_inertia_world
