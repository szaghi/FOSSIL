!run figck-rotate fig_ck_rotate
!fig ck-rotate
!< Figure of the cookbook recipe "Rotate a body": the NACA 0012 wing before (ghost) and after (body) the 10-degree
!< nose-up turn about the span axis through the quarter chord, seen along the span.
program fig_ck_rotate
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ez_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: wing, turned
integer(I4P)             :: u

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
turned = wing
call turned%rotate(axis=ez_R8P, angle=-10._R8P * PI / 180._R8P, center=vector_R8P(0._R8P, 0._R8P, 0._R8P))
call turned%analyze
call fig_open(u, 'ck-rotate', 'rotate(axis=z, angle=-10 deg, center=quarter chord): before (ghost) and after', &
              -90._R8P, 90._R8P, height=300_I4P)
call fig_surface(u, 'ghost', wing)
call fig_surface(u, 'body', turned)
call fig_points(u, 'accent', [vector_R8P(0._R8P, 0._R8P, 1.01_R8P)])
call fig_label(u, 'accent', vector_R8P(0._R8P, 0.03_R8P, 1.01_R8P), 'centre')
call fig_close(u)
endprogram fig_ck_rotate
