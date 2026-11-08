touching = 0 ; was = .false.                                                     ! 200 instants, 1/199 apart
do k = 0, 199
   s = k / 199._R8P
   call dragon%distance_to(pebble, distance=d, intersecting=hit, pose_other=pose_interpolate(p0, p1, s))
   if (hit) touching = touching + 1
   if (hit .and. .not. was) first = s                                            ! an interval of contact starts
   if (was .and. .not. hit) print '(A,F7.4,A,F7.4)', '  touching from s = ', first, ' to s = ', last
   if (hit) last = s
   was = hit
enddo
print '(A,I0,A)', '200 samples: touching at ', touching, ' of them'
