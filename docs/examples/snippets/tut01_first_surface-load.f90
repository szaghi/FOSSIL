use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: bmin, bmax

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.) ! analyze runs inside
