!run figck-body-distance fig_ck_body_distance
!fig ck-body-distance
!< Figure of the recipe "How close are my two bodies?": the unit cube a at rest and the cube b turned 45 degrees about
!< the vertical axis with its centre at (2, 0.5, 0.5), the pose of the recipe (b moved by move only to draw it), and the
!< closest pair distance_to returns: the point on a (yellow), the point on b (teal) and the segment joining them.
program fig_ck_body_distance
use fossil,     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: a, b
type(rigid_pose_t)       :: pose
type(vector_R8P)         :: qa(1), qb(1)
real(R8P)                :: d, pa(3), pb(3)
character(80)            :: caption
integer(I4P)             :: u

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
pose = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI / 4)                ! the recipe's pose
pose%t = [2._R8P, 0.5_R8P, 0.5_R8P] - matmul(pose%r, [0.5_R8P, 0.5_R8P, 0.5_R8P])
call a%distance_to(b, distance=d, point_self=pa, point_other=pb, pose_other=pose)
call b%move(pose)                                                                        ! only to draw it
qa(1) = vector_R8P(pa(1), pa(2), pa(3)) ; qb(1) = vector_R8P(pb(1), pb(2), pb(3))
write(caption, '(A,F7.4)') 'distance_to: ', d
call fig_open(u, 'ck-body-distance', trim(caption), -90._R8P, 60._R8P)
call fig_surface(u, 'body', a)
call fig_surface(u, 'body', b)
call fig_segment(u, 'accent', qa(1), qb(1))
call fig_points(u, 'foot', qa)
call fig_points(u, 'probe', qb)
call fig_close(u)
endprogram fig_ck_body_distance
