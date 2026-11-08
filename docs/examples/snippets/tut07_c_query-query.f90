use fossil,     only : curve_object
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P

implicit none
real(R8P), parameter :: PI = 4._R8P * atan(1._R8P)
type(curve_object)   :: section
real(R8P)            :: cp(2), d, ring(2, 12), dr(12), uv(2), cpk(2)
integer(I4P)         :: segment, vertex, k, s, u
type(vector_R8P)     :: a, b

call section%load_from_file(file_name='dragon-section.dat')
call section%compute_distance(p=[0._R8P, 0._R8P], distance=d, closest_point=cp, closest_segment=segment, &
                              closest_vertex=vertex)
print '(A,F10.6,A,2F10.6,A,I0,A,I0)', 'd(0,0) ', d, '  closest ', cp, '  segment ', segment, '  vertex ', vertex
do k = 1, 12                                                    ! an ellipse about the origin
   ring(:, k) = [0.12_R8P * cos(2._R8P * PI * k / 12._R8P), 0.4_R8P * sin(2._R8P * PI * k / 12._R8P)]
enddo
call section%distance_many(points=ring, distances=dr)
print '(A,I0,A,F9.5,A,F9.5)', 'ellipse: inside ', count(dr < 0._R8P), ' of 12, d from ', minval(dr), ' to ', maxval(dr)
