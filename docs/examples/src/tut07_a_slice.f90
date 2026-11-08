!run tut07 tut07_a_slice
!fig tut07-crossing
!< Tutorial 7, first run: the section of the placed dragon by the plane z = 0 as a curve_object, its metrics and
!< validity, saved for chapters 10 and 11. The figure zooms on the one pair of segments that cross, seen from above.
program tut07_a_slice
use fossil_fig, only : fig_open, fig_zoom, fig_segment, fig_points, fig_close
!region slice
use fossil,     only : surface_stl_object, curve_object, slice_surface, CURVE_VALID
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P

implicit none
type(surface_stl_object)  :: surface
type(curve_object)        :: section
integer(I4P)              :: open_loops, status, crossings, s, u
integer(I4P), allocatable :: pairs(:,:)
type(vector_R8P)          :: a, b, c
character(80)             :: caption

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call slice_surface(surface, 3_I4P, 0._R8P, section, merge_collinear=.true., open_loops=open_loops, status=status)
print '(A,I0,A,I0,A,I0,A)', 'section z=0: ', section%get_loops_number(), ' loops, ', section%get_vertices_number(), &
                            ' vertices, ', open_loops, ' open'
print '(A,F9.6,A,F9.6,A,L1)', 'area ', section%get_area(), '  perimeter ', section%get_perimeter(), &
                              '  valid ', section%check_validity() == CURVE_VALID
call section%self_intersections(pairs, crossings)
print '(A,I0)', 'self-intersecting segment pairs ', crossings
call section%save_into_file(file_name='dragon-section.dat')
!endregion slice
print '(A,I0,A,I0)', 'crossing segments ', pairs(1, 1), ' and ', pairs(2, 1)
c = 0.5_R8P * (xy(pairs(1, 1)) + xy(section%get_next(pairs(1, 1))))
write(caption, '(A,I0,A,I0,A)') 'a fold of the section, seen from above: segments ', pairs(1, 1), ' and ', &
                                pairs(2, 1), ' (ember) cross'
call fig_open(u, 'tut07-crossing', trim(caption), 0._R8P, 90._R8P)
call fig_zoom(u, c, 0.006_R8P)
do s = 1, section%get_vertices_number()
   a = xy(s) ; b = xy(section%get_next(s))
   call fig_segment(u, merge('hole  ', 'accent', any(pairs(:, 1) == s)), a, b)
enddo
call fig_points(u, 'foot', [(xy(s), s = 1, section%get_vertices_number())])
call fig_close(u)
contains
   function xy(s) result(q)
   !< Vertex s of the section as a point of the plane z = 0.
   integer(I4P), intent(in) :: s    !< Vertex.
   type(vector_R8P)         :: q    !< Point.
   real(R8P)                :: w(2) !< (u, v) = (x, y).

   w = section%get_vertex(s) ; q = w(1) * ex_R8P + w(2) * ey_R8P
   endfunction xy
endprogram tut07_a_slice
