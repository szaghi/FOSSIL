!run figck-curve-corners fig_ck_curve_corners
!fig ck-curve-corners
!< Figure of the cookbook recipe "Find the corners of a section": the trailing edge of the NACA 0012 section, zoomed,
!< its vertices faint and the two corners feature_corners finds at the blunt base in ember.
program fig_ck_curve_corners
use fossil,     only : surface_stl_object, curve_object, slice_surface
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)      :: wing
type(curve_object)            :: section
integer(I4P), allocatable     :: corners(:)
type(vector_R8P), allocatable :: a(:), b(:), c(:)
real(R8P)                     :: uv(2)
integer(I4P)                  :: i, n, u
character(80)                 :: caption

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, section)
call section%feature_corners(corners=corners, angle=30._R8P)
n = section%get_vertices_number() ; allocate(a(n), b(n))
do i = 1, n
   uv = section%get_vertex(i) ; a(i) = vector_R8P(uv(1), uv(2), 0._R8P)
   uv = section%get_vertex(section%get_next(i)) ; b(i) = vector_R8P(uv(1), uv(2), 0._R8P)
enddo
allocate(c(size(corners)))
do i = 1, size(corners)
   uv = section%get_vertex(corners(i)) ; c(i) = vector_R8P(uv(1), uv(2), 0._R8P)
enddo
write(caption, '(A,I0,A,I0,A)') 'feature_corners(angle=30): ', size(corners), ' of ', n, &
                                ' vertices, at the blunt trailing edge'
call fig_open(u, 'ck-curve-corners', trim(caption), -90._R8P, 90._R8P, height=320_I4P)
call fig_zoom(u, vector_R8P(0.745_R8P, 0._R8P, 0._R8P), 0.007_R8P)
call fig_segments(u, 'body', a, b)
call fig_points(u, 'ghost', a)
call fig_points(u, 'hole', c)
call fig_close(u)
endprogram fig_ck_curve_corners
