bodies(1) = dragon ; poses0(1) = pose_identity() ; poses1(1) = poses0(1)       ! the dragon at rest
bodies(2:4) = pebble                                                             ! three pebbles, three throws
poses0(2) = p0 ; poses1(2) = p1
poses0(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=0._R8P, translation=[0._R8P, 0.3_R8P, 0.6_R8P])
poses1(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=PI, translation=[0._R8P, 0.3_R8P, -0.6_R8P])
poses0(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, translation=[-0.15_R8P, 0.3_R8P, 0.45_R8P])
poses1(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI, translation=[1.05_R8P, 0.3_R8P, 0.45_R8P])
call first_impact(bodies, poses0, poses1, toi, pair, status, tolerance=1.e-6_R8P, tested=tested)
print '(A,L1,A,F9.6,A,2I2,A,I0,A)', 'first contact ', status == CCD_CONTACT, ' at s = ', toi, ', bodies', pair, &
                                    ' (', tested, ' of 6 pairs tested)'
toi_all = 1._R8P
do i = 1, 3
   do j = i + 1, 4
      call bodies(i)%time_of_impact(bodies(j), poses0(i), poses1(i), poses0(j), poses1(j), s, status, &
                                     tolerance=1.e-6_R8P)
      if (status == CCD_CONTACT) then
         toi_all = min(toi_all, s)
         print '(A,I0,A,I0,A,F9.6)', '  pair ', i, '-', j, ': contact at s = ', s
      else
         print '(A,I0,A,I0,A)', '  pair ', i, '-', j, ': separate'
      endif
   enddo
enddo
print '(A,L1)', 'equal to the loop over all pairs, bit for bit: ', toi_all == toi
