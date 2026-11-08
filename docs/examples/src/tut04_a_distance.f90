!run tut04 tut04_a_distance
!fig tut04-axis
!< Tutorial 4, first run: signed distance of five points on the y axis, by two sign algorithms, and the inside test.
!< The figure shows the points, teal inside and ember outside, each joined to its closest point on the dragon.
program tut04_a_distance
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_points, fig_label, fig_close
!region distance
use fossil,     only : surface_stl_object, SIGN_FAST_WINDING_NUMBER
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: p(0:4), foot(0:4)
real(R8P)                :: d_pn, d_wn
logical                  :: inside(0:4)
integer(I4P)             :: k, u
character(16)            :: label

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
print '(A)', '    y      d (pseudo-n.)  d (winding n.) inside'
do k = 0, 4
   p(k) = 0._R8P * ex_R8P + (-0.6_R8P + 0.3_R8P * k) * ey_R8P + 0._R8P * ez_R8P
   d_pn = surface%distance(point=p(k), is_signed=.true., is_square_root=.true.)
   d_wn = surface%distance(point=p(k), is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   inside(k) = surface%is_point_inside(point=p(k), sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   print '(F7.3,2F15.6,L5)', p(k)%y, d_pn, d_wn, inside(k)
enddo
!endregion distance
call fig_open(u, 'tut04-axis', 'the y axis: points inside (teal), outside (ember), closest points (gold)', &
              -20._R8P, 15._R8P)
call fig_surface(u, 'ghost', surface)
call fig_segment(u, 'box', -0.6_R8P * ey_R8P, 0.6_R8P * ey_R8P)
do k = 0, 4
   call surface%compute_distance(point=p(k), distance=d_pn, closest_point=foot(k))
   call fig_segment(u, 'accent', p(k), foot(k))
   write(label, '(A,F4.1)') 'y ', real(nint(10._R8P * p(k)%y), R8P) / 10._R8P
   call fig_label(u, 'box', p(k), trim(label))
enddo
call fig_points(u, 'inside', pack(p, inside))
call fig_points(u, 'outside', pack(p, .not. inside))
call fig_points(u, 'foot', foot)
call fig_close(u)
endprogram tut04_a_distance
