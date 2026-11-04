!run fig07 fig_tut07
!fig tut07-section
!fig tut07-field
!< Figures of tutorial 7: the section of the dragon by the plane z = 0 and the 2-D signed distance of that curve.
!< The dragon is placed as in tutorials 2 and 3 (sanitize with fill_holes, centre, quarter turn, longest side 1,
!< saved and read back as dragon-placed.stl is), sliced by slice_surface, and the field is the 81^2 one of tutorial 7.
program fig_tut07
use fossil,     only : surface_stl_object, curve_object, distance_field_object, slice_surface
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_segment, fig_segments, fig_label, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter          :: pi = 4._R8P * atan(1._R8P)
real(R8P), parameter          :: s = 0.05_R8P                     ! colour scale of the field: tanh(phi / s)
type(surface_stl_object)      :: surface, lower, upper
type(curve_object)            :: section
type(distance_field_object)   :: field
type(vector_R8P)              :: bmin, bmax, o
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: x(81), uv(2), scale, h
integer(I4P)                  :: u, i, j, n, status

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)        ! dragon-placed.stl, as tutorials 2-3
call surface%sanitize(fill_holes=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))
call surface%rotate(axis=ez_R8P, angle=pi / 2._R8P)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)
call surface%analyze
call surface%save_into_file(file_name='fig07-dragon-placed.stl')
call surface%load_from_file(file_name='fig07-dragon-placed.stl', guess_format=.true.)

call slice_surface(surface, 3_I4P, 0._R8P, section, merge_collinear=.true., status=status)
n = section%get_vertices_number()
allocate(a(n), b(n))
do i = 1, n                                                                       ! segment i: vertex i to its next
   uv = section%get_vertex(i) ; a(i) = uv(1) * ex_R8P + uv(2) * ey_R8P
   uv = section%get_vertex(section%get_next(i)) ; b(i) = uv(1) * ex_R8P + uv(2) * ey_R8P
enddo

bmin = surface%get_bmin() ; bmax = surface%get_bmax()
lower = surface
call lower%clip(bmin=bmin, bmax=bmax%x * ex_R8P + bmax%y * ey_R8P + 0._R8P * ez_R8P, remainder=upper)
call fig_open(u, 'tut07-section', 'slice_surface(z = 0): the section, two closed loops', -60._R8P, 50._R8P)
call fig_surface(u, 'body', lower)
call fig_surface(u, 'ghost', upper)
call fig_box(u, 'box', bmin%x * ex_R8P + bmin%y * ey_R8P, bmax%x * ex_R8P + bmax%y * ey_R8P)  ! the plane z = 0
call fig_segments(u, 'hole', a, b)
call fig_close(u)

x = [(-0.6_R8P + 0.015_R8P * i, i = 0, 80)]
call field%initialize(x=x, y=x)
call field%compute(curve=section, status=status)
h = 0.5_R8P * 0.015_R8P                                                           ! nodes around the curve only
call fig_open(u, 'tut07-field', 'compute(curve=): the 2-D signed distance, colour tanh(phi/0.05)', 0._R8P, 90._R8P)
do j = 15, 76
   do i = 21, 55                                                                  ! the node's cell, coloured by phi
      call fig_quad(u, 'field', (x(i) - h) * ex_R8P + (x(j) - h) * ey_R8P, (x(i) + h) * ex_R8P + (x(j) - h) * ey_R8P, &
                    (x(i) + h) * ex_R8P + (x(j) + h) * ey_R8P, (x(i) - h) * ex_R8P + (x(j) + h) * ey_R8P, &
                    value=tanh(field%phi(i, j, 1) / s))
   enddo
enddo
o = (x(21) + 0.02_R8P) * ex_R8P + (x(15) + 0.02_R8P) * ey_R8P + 0.001_R8P * ez_R8P  ! axes glyph at a corner
call fig_segment(u, 'box', o, o + 0.08_R8P * ex_R8P) ; call fig_label(u, 'box', o + 0.08_R8P * ex_R8P, 'x')
call fig_segment(u, 'box', o, o + 0.08_R8P * ey_R8P) ; call fig_label(u, 'box', o + 0.08_R8P * ey_R8P, 'y')
call fig_segments(u, 'accent', a + 0.001_R8P * ez_R8P, b + 0.001_R8P * ez_R8P)
call fig_close(u)
endprogram fig_tut07
