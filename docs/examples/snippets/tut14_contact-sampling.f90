! one step: thrown along x through the dragon's body, turning a quarter turn about y
p0 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=0._R8P, translation=[-0.6_R8P, -0.3_R8P, -0.2_R8P])
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=PI / 2, translation=[0.6_R8P, -0.3_R8P, -0.2_R8P])
do k = 0, 4                                                                      ! the gap at five instants
   s = 0.25_R8P * k
   call dragon%distance_to(pebble, distance=d, intersecting=hit, pose_other=pose_interpolate(p0, p1, s))
   print '(A,F5.2,A,F8.5,A,L1)', 's = ', s, ': gap ', d, ', touching ', hit
enddo
