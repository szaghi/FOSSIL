!< Figure of the cookbook recipe "Snap points onto a body with sharp edges": the unit cube (translucent), the five points of
!< ck_project_points (teal) joined to their feet (yellow), each foot labelled with the constraint classify_points gave
!< it: the face, the edge x = 1, z = 1, or the corner (1, 1, 1).
program fig_ck_project_points
use fossil,     only : surface_stl_object, projection_target_t, classify_points, project_points
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_points, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
character(6), parameter   :: KIND_NAME(0:3) = [character(6) :: 'free', 'face', 'edge', 'corner']
type(surface_stl_object)  :: cube
type(projection_target_t) :: target
real(R8P)                 :: x(5), y(5), z(5), fx(5), fy(5), fz(5)
integer(I4P)              :: kinds(5), ids(5), i, u
type(vector_R8P)          :: p(5), f(5)

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%projection_target(target)
x = [0.40_R8P, 0.92_R8P, 0.92_R8P, 1.08_R8P, 0.93_R8P]
y = [0.50_R8P, 0.50_R8P, 0.20_R8P, 0.65_R8P, 1.07_R8P]
z = [1.08_R8P, 1.07_R8P, 0.92_R8P, 1.04_R8P, 0.95_R8P]
call classify_points(target, x, y, z, tolerance=0.12_R8P, kinds=kinds, ids=ids)
call project_points(target, x, y, z, kinds, ids, fx, fy, fz)
do i = 1, 5
   p(i) = vector_R8P(x(i), y(i), z(i)) ; f(i) = vector_R8P(fx(i), fy(i), fz(i))
enddo
call fig_open(u, 'ck-project-points', 'classify_points + project_points: face, edge or corner, tolerance 0.12', &
              -30._R8P, 35._R8P)
call fig_zoom(u, vector_R8P(0.75_R8P, 0.6_R8P, 0.95_R8P), 0.42_R8P)
call fig_surface(u, 'ghost', cube)
call fig_segments(u, 'probe', p, f)
call fig_points(u, 'probe', p)
call fig_points(u, 'foot', f)
do i = 1, 5
   call fig_label(u, 'foot', f(i) + vector_R8P(0.015_R8P, -0.03_R8P, 0.02_R8P), trim(KIND_NAME(kinds(i))))
enddo
call fig_close(u)
endprogram fig_ck_project_points
