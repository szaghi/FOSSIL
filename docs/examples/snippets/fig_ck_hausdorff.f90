!< Figure of the cookbook recipe "Compare two versions of a part": the NACA 0012 wing (ghost) and the same wing at 2
!< degrees more incidence (body), seen along the span and zoomed on the trailing edge, where the Hausdorff distance
!< is reached.
program fig_ck_hausdorff
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ez_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: wing, turned
type(vector_R8P)         :: te
real(R8P)                :: lower, upper
integer(I4P)             :: u
character(90)            :: caption

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
turned = wing
call turned%rotate(axis=ez_R8P, angle=-2._R8P * PI / 180._R8P, center=vector_R8P(0._R8P, 0._R8P, 0._R8P))
call turned%analyze
call wing%hausdorff(other=turned, lower=lower, upper=upper, tolerance=1e-6_R8P)
te = vector_R8P(0.75_R8P * cos(2._R8P * PI / 180._R8P), -0.75_R8P * sin(2._R8P * PI / 180._R8P), 1._R8P)
write(caption, '(A,F9.7,A)') 'the wing (ghost) and the wing at +2 deg: Hausdorff distance ', upper
call fig_open(u, 'ck-hausdorff', trim(caption), -90._R8P, 90._R8P, height=320_I4P)
call fig_zoom(u, vector_R8P(0.69_R8P, -0.015_R8P, 1._R8P), 0.07_R8P)
call fig_surface(u, 'ghost', wing)
call fig_surface(u, 'body', turned)
call fig_points(u, 'accent', [te])
call fig_label(u, 'accent', te + vector_R8P(-0.045_R8P, -0.012_R8P, 0._R8P), 'trailing edge at +2 deg')
call fig_close(u)
endprogram fig_ck_hausdorff
