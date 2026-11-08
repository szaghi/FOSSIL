do i = 1, size(points)                                     ! 21 points on a line through the dragon
   points(i) = (-12._R8P + 1.1_R8P * (i - 1)) * ex_R8P + 0.5_R8P * ey_R8P + 12._R8P * ez_R8P
enddo
call dragon%distance_many(points=points, distances=d_bvh, is_signed=.true., is_square_root=.true., status=status)
if (status /= STATUS_OK) error stop 'distance_many: points and distances differ in size'
same = .true.
do i = 1, size(points)
   same = same .and. (d_bvh(i) == dragon%distance(point=points(i), is_signed=.true., is_square_root=.true.))
enddo
print '(A,I0,A,L1)', 'distance_many: ', size(points), ' points, equal to distance bit for bit ', same
print '(A,I0)', '  inside ', count(d_bvh < 0._R8P)
