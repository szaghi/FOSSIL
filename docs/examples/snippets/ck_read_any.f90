program ck_read_any
use fossil, only : surface_stl_object
implicit none
type(surface_stl_object) :: binary, ascii

call binary%load_from_file(file_name='dragon.stl', guess_format=.true.)
call ascii%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
print '(A,I0,A)', 'dragon.stl:         ', binary%get_facets_number(), ' facets'
print '(A,I0,A)', 'naca0012-ascii.stl: ', ascii%get_facets_number(), ' facets'
endprogram ck_read_any
