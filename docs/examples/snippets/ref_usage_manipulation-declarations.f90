use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P, rotation_matrix_R8P, mirror_matrix_R8P

implicit none
real(R8P), parameter     :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: dragon, surface, remainder, surface_1, surface_2
type(vector_R8P)         :: delta, axis, normal, factor, bmin, bmax
real(R8P)                :: angle, matrix(3,3)

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
