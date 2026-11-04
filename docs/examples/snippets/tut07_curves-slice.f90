use fossil, only : surface_stl_object, curve_object, distance_field_object, slice_surface, CURVE_VALID, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(surface_stl_object)    :: surface
type(curve_object)          :: section
type(distance_field_object) :: field
real(R8P)                   :: x(81), cp(2), d
integer(I4P)                :: i, open_loops, segment, vertex, status, crossings
integer(I4P), allocatable   :: pairs(:,:)

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call slice_surface(surface, 3_I4P, 0._R8P, section, merge_collinear=.true., open_loops=open_loops, status=status)
print '(A,I0,A,I0,A,I0,A)', 'section z=0: ', section%get_loops_number(), ' loops, ', section%get_vertices_number(), &
                            ' vertices, ', open_loops, ' open'
print '(A,F9.6,A,F9.6,A,L1)', 'area ', section%get_area(), '  perimeter ', section%get_perimeter(), &
                              '  valid ', section%check_validity() == CURVE_VALID
call section%self_intersections(pairs, crossings)
print '(A,I0)', 'self-intersecting segment pairs ', crossings
call section%save_into_file(file_name='dragon-section.dat')
