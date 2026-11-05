use fossil, only : surface_stl_object
use penf,   only : I4P, R8P

implicit none
type(surface_stl_object) :: surface, original
integer(I4P)             :: facets_before
real(R8P)                :: lower, upper

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
facets_before = surface%get_facets_number()
