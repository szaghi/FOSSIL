program ck_remesh
use fossil, only : surface_stl_object, REM_STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: cube
integer(I4P)             :: status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
print '(A,I0,A,F7.4,A,F7.4)', 'before: ', cube%get_facets_number(), ' facets, edges ', cube%smallest_edge_len(), &
      ' to ', cube%largest_edge_len()
call cube%isotropic_remesh(target_length=0.1_R8P, iterations=10_I4P, preserve_features=.true., status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
print '(A,I0,A,F7.4,A,F7.4)', 'after:  ', cube%get_facets_number(), ' facets, edges ', cube%smallest_edge_len(), &
      ' to ', cube%largest_edge_len()
print '(A,F8.5)', 'volume  ', cube%get_volume()
endprogram ck_remesh
