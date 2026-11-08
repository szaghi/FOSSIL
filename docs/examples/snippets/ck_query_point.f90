program ck_query_point
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
character(6), parameter  :: KIND_NAME(3) = [character(6) :: 'face', 'edge', 'vertex'] ! FEATURE_FACE, _EDGE, _VERTEX
type(surface_stl_object) :: cube
type(vector_R8P)         :: p(3), foot, normal
real(R8P)                :: d
integer(I4P)             :: i, feature

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
p = [vector_R8P(0.3_R8P, 0.15_R8P, 1.3_R8P), &   ! above the top face
     vector_R8P(1.3_R8P, 0.3_R8P, 1.3_R8P), &   ! off the edge x = 1, z = 1
     vector_R8P(1.3_R8P, 1.3_R8P, 1.3_R8P)]     ! off the corner (1, 1, 1)
do i = 1, 3
   call cube%query_point(point=p(i), distance=d, closest_point=foot, feature=feature, normal=normal)
   print '(A,I0,A,F7.4,A,3F7.3,3A,3F7.3)', 'point ', i, ': d ', d, ', foot', foot%x, foot%y, foot%z, &
         ', ', KIND_NAME(feature), ', normal', normal%x, normal%y, normal%z
enddo
endprogram ck_query_point
