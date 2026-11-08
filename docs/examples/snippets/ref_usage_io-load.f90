use fossil, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
type(surface_stl_object) :: surface, wing, copy, head
type(vector_R8P)         :: bmin, bmax
integer(I4P)             :: stat, u

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)  ! ASCII or binary, from the file size
print '(A,I0)', 'dragon.stl             facets ', surface%get_facets_number()
