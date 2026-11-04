!run ck_closest_point ck_closest_point
program ck_closest_point
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: cube
type(vector_R8P)         :: foot, normal
real(R8P)                :: d
integer(I4P)             :: facet

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%compute_distance(point=vector_R8P(1.5_R8P, 0.5_R8P, 2._R8P), distance=d, is_signed=.true., &
                           is_square_root=.true., facet_index=facet, closest_point=foot, normal=normal)
print '(A,F8.5)',  'signed distance ', d
print '(A,3F8.4)', 'closest point   ', foot%x, foot%y, foot%z
print '(A,3F8.4)', 'normal          ', normal%x, normal%y, normal%z
print '(A,I0)',    'facet           ', facet
endprogram ck_closest_point
