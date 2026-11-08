do k = 0, 4                                                                      ! the gap at five instants
   s = 0.25_R8P * k
   call dragon%distance_to(pebble, distance=d, intersecting=hit, pose_other=pose_interpolate(p0, p1, s))
   print '(A,F5.2,A,F8.5,A,L1)', 's = ', s, ': gap ', d, ', touching ', hit
enddo
