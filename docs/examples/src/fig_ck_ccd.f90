!run figck-ccd fig_ck_ccd
!fig ck-ccd
!< Figure of the recipe "When do my moving bodies first touch?": the cube b at rest, the cube a at the start of its
!< step (ghost; at s = 1 it would sit exactly on b) and at the time of first contact (solid), the witness points on a
!< (yellow) and b (teal).
program fig_ck_ccd
use fossil,     only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, pose_interpolate
use fossil_fig, only : fig_open, fig_surface, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: a, b, am
type(rigid_pose_t)       :: a0, a1, rest
type(vector_R8P)         :: qa(1), qb(1)
real(R8P)                :: toi, pa(3), pb(3)
character(80)            :: caption
integer(I4P)             :: status, u

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(x=3._R8P) ; call b%analyze
rest = pose_identity() ; a0 = pose_identity()
a1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI / 2, translation=[4._R8P, 0._R8P, 0._R8P])
call a%time_of_impact(b, a0, a1, rest, rest, toi, status, point_self=pa, point_other=pb)
qa(1) = vector_R8P(pa(1), pa(2), pa(3)) ; qb(1) = vector_R8P(pb(1), pb(2), pb(3))
write(caption, '(A,F7.4)') 'time_of_impact: s = ', toi
call fig_open(u, 'ck-ccd', trim(caption), 100._R8P, 20._R8P)
am = a ; call am%move(a0) ; call fig_surface(u, 'ghost', am)
am = a ; call am%move(pose_interpolate(a0, a1, toi)) ; call fig_surface(u, 'body', am)
call fig_surface(u, 'body', b)
call fig_points(u, 'foot', qa)
call fig_points(u, 'probe', qb)
call fig_close(u)
endprogram fig_ck_ccd
