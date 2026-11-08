!run figck-curve-markers fig_ck_curve_markers
!fig ck-curve-markers
!< Figure of the cookbook recipe "Place markers on a curve": the NACA 0012 section and the 41 markers curve_markers
!< places for a grid of spacing 0.05, each with its outward unit normal (scaled).
program fig_ck_curve_markers
use fossil,     only : surface_stl_object, curve_object, slice_surface
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)      :: wing
type(curve_object)            :: section
real(R8P), allocatable        :: centroids(:,:), lengths(:), normals(:,:)
type(vector_R8P), allocatable :: a(:), b(:), m(:), t(:)
real(R8P)                     :: uv(2)
integer(I4P)                  :: i, n, u
character(80)                 :: caption

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, section)
call section%curve_markers(dx=0.05_R8P, centroids=centroids, lengths=lengths, normals=normals)
n = section%get_vertices_number() ; allocate(a(n), b(n))
do i = 1, n
   uv = section%get_vertex(i) ; a(i) = vector_R8P(uv(1), uv(2), 0._R8P)
   uv = section%get_vertex(section%get_next(i)) ; b(i) = vector_R8P(uv(1), uv(2), 0._R8P)
enddo
n = size(lengths) ; allocate(m(n), t(n))
do i = 1, n
   m(i) = vector_R8P(centroids(1, i), centroids(2, i), 0._R8P)
   t(i) = m(i) + vector_R8P(normals(1, i), normals(2, i), 0._R8P) * 0.04_R8P
enddo
write(caption, '(A,I0,A,F8.6,A)') 'curve_markers(dx=0.05): ', n, ' markers of arc length ', lengths(1), ', normals'
call fig_open(u, 'ck-curve-markers', trim(caption), -90._R8P, 90._R8P, height=300_I4P)
call fig_segments(u, 'body', a, b)
call fig_segments(u, 'probe', m, t)
call fig_points(u, 'foot', m)
call fig_close(u)
endprogram fig_ck_curve_markers
