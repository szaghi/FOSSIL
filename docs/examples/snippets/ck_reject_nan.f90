program ck_reject_nan
use fossil, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,   only : I4P
implicit none
type(surface_stl_object) :: surface
integer(I4P)             :: u, status

open(newunit=u, file='broken.stl', action='write')
write(u, '(A)') 'solid broken', ' facet normal 0 0 1', '  outer loop', '   vertex 0 0 0', '   vertex 1 0 0', &
                '   vertex NaN 1 0', '  endloop', ' endfacet', 'endsolid broken'
close(u)
call surface%load_from_file(file_name='broken.stl', is_ascii=.true., status=status)
if (status == STATUS_INVALID_INPUT) print '(A)', 'broken.stl rejected: non-finite coordinates'
print '(A,I0,A,I0)', 'status ', status, ', facets ', surface%get_facets_number()
endprogram ck_reject_nan
