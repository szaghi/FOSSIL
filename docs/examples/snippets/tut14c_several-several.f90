bodies(1) = dragon ; poses0(1) = pose_identity() ; poses1(1) = poses0(1)       ! the dragon at rest
bodies(2:4) = pebble                                                             ! three pebbles, three throws
poses0(2) = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=0._R8P, translation=[-0.6_R8P, -0.3_R8P, -0.2_R8P])
poses1(2) = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=PI / 2, translation=[0.6_R8P, -0.3_R8P, -0.2_R8P])
poses0(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=0._R8P, translation=[0._R8P, 0.3_R8P, 0.6_R8P])
poses1(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=PI, translation=[0._R8P, 0.3_R8P, -0.6_R8P])
poses0(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, translation=[-0.15_R8P, 0.3_R8P, 0.45_R8P])
poses1(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI, translation=[1.05_R8P, 0.3_R8P, 0.45_R8P])
call first_impact(bodies, poses0, poses1, toi, pair, status, tolerance=1.e-6_R8P, tested=tested)
print '(A,L1,A,F9.6,A,2I2,A,I0,A)', 'first contact ', status == CCD_CONTACT, ' at s = ', toi, ', bodies', pair, &
                                    ' (', tested, ' of 6 pairs tested)'
