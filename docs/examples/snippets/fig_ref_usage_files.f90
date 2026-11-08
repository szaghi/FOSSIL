!< Figures of the usage reference, "Printing statistics" and "Saving an STL file", on the cube fixture after sanitize:
!< the quantities statistics() prints (bounding box, centroid, volume, facets) with the outward facet normals that make
!< the volume positive; then what save_into_file writes: every facet on its own, with its normal and three vertices,
!< drawn exploded (each facet shrunk towards its centroid), the file loaded back with the same facets and volume.
program fig_ref_usage_files
use fossil,     only : surface_stl_object, facet_object
use fossil_fig, only : fig_open, fig_triangle, fig_surface, fig_segments, fig_box, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P), parameter          :: DEG = acos(-1._R8P) / 180._R8P
type(surface_stl_object)      :: cube, back
type(facet_object), pointer   :: f
type(vector_R8P), allocatable :: a(:), b(:)
type(vector_R8P)              :: c, v(3), view
character(120)                :: caption
integer(I4P)                  :: u, i, j, n, distinct, front

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%sanitize
n = cube%get_facets_number()
allocate(a(n), b(n))
do i = 1, n
   f => cube%facet_at(i)
   a(i) = f%centroid ; b(i) = f%centroid + 0.18_R8P * f%normal
enddo
c = cube%get_centroid()
write(caption, '(A,I0,A,F6.4,A)') 'statistics() after sanitize: ', n, ' facets, volume ', cube%get_volume(), &
                                  ', every normal outward'
print '(A)', trim(caption)
call fig_open(u, 'ref-usage-statistics', trim(caption), -55._R8P, 24._R8P)
call fig_surface(u, 'ghost', cube)
call fig_segments(u, 'probe', a, b)
call fig_box(u, 'box', cube%get_bmin() - 0.04_R8P, cube%get_bmax() + 0.04_R8P)
call fig_points(u, 'accent', [c])
call fig_label(u, 'accent', c, 'centroid (0.5, 0.5, 0.5)')
call fig_label(u, 'box', cube%get_bmin() - 0.04_R8P, 'X, Y, Z extents [0, 1]')
call fig_close(u)

call cube%save_into_file(file_name='cube-saved.stl')
call back%load_from_file(file_name='cube-saved.stl', guess_format=.true.)
distinct = 0
do i = 1, n
   f => cube%facet_at(i)
   distinct = max(distinct, maxval(f%vertex_id))
enddo
write(caption, '(A,I0,A,I0,A,I0,A)') 'save_into_file: ', n, ' facets, ', 3 * n, ' vertex records for ', distinct, &
                                     ' distinct vertices'
print '(A)', trim(caption)
print '(A,I0,A,F6.4)', 'loaded back: ', back%get_facets_number(), ' facets, volume ', back%get_volume()
call fig_open(u, 'ref-usage-saving', trim(caption), -55._R8P, 24._R8P)
front = 0
view = vector_R8P(cos(24._R8P * DEG) * cos(-55._R8P * DEG), cos(24._R8P * DEG) * sin(-55._R8P * DEG), &
                  sin(24._R8P * DEG))
do i = 1, n
   f => cube%facet_at(i)
   do j = 1, 3
      v(j) = f%centroid + 0.78_R8P * (f%vertex(j) - f%centroid) + 0.06_R8P * f%normal
   enddo
   if ((f%normal .dot. view) <= 0._R8P) cycle                      ! the facets facing the camera
   call fig_triangle(u, 'body', v(1), v(2), v(3))
   call fig_points(u, 'foot', v)
   if (front == 0 .and. f%normal%z > 0.5_R8P .and. f%centroid%x < 0.5_R8P .and. f%centroid%y < 0.5_R8P) front = i
enddo
f => cube%facet_at(front)
call fig_segments(u, 'probe', [f%centroid + 0.06_R8P * f%normal], [f%centroid + 0.26_R8P * f%normal])
call fig_label(u, 'probe', f%centroid + 0.26_R8P * f%normal, 'facet normal')
call fig_label(u, 'foot', f%centroid + 0.78_R8P * (f%vertex(1) - f%centroid) + 0.06_R8P * f%normal, 'vertex')
call fig_close(u)
endprogram fig_ref_usage_files
