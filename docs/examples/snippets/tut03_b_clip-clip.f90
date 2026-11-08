use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object) :: surface, half, rest, twin
type(vector_R8P)         :: bmin, bmax

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
half = surface
call half%clip(bmin=bmin, bmax=bmax%x * ex_R8P + 0._R8P * ey_R8P + bmax%z * ez_R8P, remainder=rest)  ! the box y <= 0
print '(A,I0,A,I0,A,I0)', 'clip y<=0 ', half%get_facets_number(), ' facets, remainder ', rest%get_facets_number(), &
                          ', sum ', half%get_facets_number() + rest%get_facets_number()
