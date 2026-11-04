!run ck_load_clipped ck_load_clipped
program ck_load_clipped
use fossil, only : surface_stl_object
use penf,   only : R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: head
type(vector_R8P)         :: bmin, bmax

call head%load_from_file(file_name='dragon.stl', guess_format=.true.,               &
                         clip_min=vector_R8P(-11._R8P, -5._R8P, 15._R8P), &
                         clip_max=vector_R8P( 10._R8P,  5._R8P, 20._R8P))
bmin = head%get_bmin() ; bmax = head%get_bmax()
print '(A,I0)',      'facets kept ', head%get_facets_number()
print '(A,2F9.4)',   'z range     ', bmin%z, bmax%z
endprogram ck_load_clipped
