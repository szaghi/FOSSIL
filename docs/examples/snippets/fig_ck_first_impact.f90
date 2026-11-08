!< Figure of the cookbook recipe "Which of my moving bodies touch first?": the four cubes of ck_first_impact at the
!< first contact s = toi, the pair in contact (3 and 4) in ember, the others in stone; the two moving cubes at s = 0
!< translucent; seen from above.
program fig_ck_first_impact
use fossil,     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, pose_interpolate, first_impact
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: cubes(4), moved
type(rigid_pose_t)       :: p0(4), p1(4)
real(R8P)                :: toi
integer(I4P)             :: i, pair(2), status, u
character(80)            :: caption

do i = 1, 4
   call cubes(i)%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   p0(i) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [2._R8P * (i - 1), 0._R8P, 0._R8P]) ; p1(i) = p0(i)
enddo
p1(2)%t = [2._R8P, 2._R8P, 0._R8P]
p1(4) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], PI / 4, [4.5_R8P, 0._R8P, 0._R8P])
call first_impact(cubes, p0, p1, toi, pair, status)
write(caption, '(A,F8.6,A,I0,A,I0,A)') 'first contact at s = ', toi, ': cubes ', pair(1), ' and ', pair(2), &
                                       '; moving cubes at s = 0 faint'
call fig_open(u, 'ck-first-impact', trim(caption), -90._R8P, 90._R8P, height=300_I4P)
do i = 2, 4, 2
   moved = cubes(i) ; call moved%move(p0(i)) ; call fig_surface(u, 'ghost', moved)
enddo
do i = 1, 4
   moved = cubes(i) ; call moved%move(pose_interpolate(p0(i), p1(i), toi))
   call fig_surface(u, merge('patch', 'body ', any(pair == i)), moved)
enddo
call fig_close(u)
endprogram fig_ck_first_impact
