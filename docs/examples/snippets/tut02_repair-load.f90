use fossil, only : surface_stl_object
use penf,   only : I4P

implicit none
type(surface_stl_object) :: surface
integer(I4P)             :: facets_before

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
facets_before = surface%get_facets_number()
