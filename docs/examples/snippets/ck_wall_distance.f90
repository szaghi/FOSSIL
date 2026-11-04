program ck_wall_distance
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: cube
type(vector_R8P)         :: points(4)
real(R8P)                :: d(4)
integer(I4P)             :: body(4), i

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
points = [vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(1.1_R8P, 0.5_R8P, 0.5_R8P), &
          vector_R8P(1.3_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(3.0_R8P, 3.0_R8P, 3.0_R8P)]
call cube%wall_distance_many(points=points, distances=d, bodies=body, cutoff=0.25_R8P)
do i = 1, 4
   print '(A,3F5.1,A,F7.4,A,I0)', 'point', points(i)%x, points(i)%y, points(i)%z, '  distance', d(i), '  body ', body(i)
enddo
endprogram ck_wall_distance
