!run ck_decimate ck_decimate
program ck_decimate
use fossil, only : surface_stl_object, DEC_STATUS_OK
use penf,   only : I4P
implicit none
type(surface_stl_object) :: dragon
integer(I4P)             :: status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
print '(A,I0,A,F9.3)', 'loaded:    ', dragon%get_facets_number(), ' facets, area ', dragon%get_area()
call dragon%decimate(target_facets=dragon%get_facets_number() / 4, status=status)
if (status /= DEC_STATUS_OK) error stop 'decimate failed'
print '(A,I0,A,F9.3)', 'decimated: ', dragon%get_facets_number(), ' facets, area ', dragon%get_area()
call dragon%save_into_file(file_name='dragon-coarse.stl', is_ascii=.false.)
endprogram ck_decimate
