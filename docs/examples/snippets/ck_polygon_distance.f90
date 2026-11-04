program ck_polygon_distance
use fossil, only : curve_object
use penf,   only : I4P, R8P
implicit none
type(curve_object) :: square
real(R8P)          :: d, foot(2)

! the unit square, counter-clockwise: the solid is on the left of every segment
call square%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
print '(A,F8.4)', 'd(0.5, 0.25) ', square%distance([0.5_R8P, 0.25_R8P])
print '(A,F8.4)', 'd(2.0, 2.00) ', square%distance([2._R8P, 2._R8P])
call square%compute_distance([1.5_R8P, 0.5_R8P], distance=d, closest_point=foot)
print '(A,F8.4,A,2F8.4)', 'd(1.5, 0.50) ', d, ', closest point', foot
endprogram ck_polygon_distance
