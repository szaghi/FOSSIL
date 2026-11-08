call wall%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[-0.1_R8P, -0.1_R8P, 0._R8P, 0._R8P], loop_sizes=[4_I4P])
do k = 1, 2                                                                      ! sliding 0.8 along the wall
   gap = merge(0.02_R8P, 4.e-6_R8P, k == 1)                                      ! at a gap of 0.02, then 4 tolerances
   q0 = pose_2d(angle=0._R8P, translation=[0.1_R8P, 0.03_R8P + gap])
   q1 = pose_2d(angle=0._R8P, translation=[0.9_R8P, 0.03_R8P + gap])
   call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                            tolerance=1.e-6_R8P, iterations=steps, directional=.false.)
   print '(A,ES8.1,A,A,A,F9.6,A,I0,A)', 'gap ', gap, ', plain:       ', status_name(status), ' at s = ', toi, &
                                      ' (', steps, ' steps)'
   call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                            tolerance=1.e-6_R8P, iterations=steps)
   print '(A,ES8.1,A,A,A,F9.6,A,I0,A)', 'gap ', gap, ', default:     ', status_name(status), ' at s = ', toi, &
                                      ' (', steps, ' steps)'
   call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                            tolerance=1.e-6_R8P, iterations=steps, directional=.true.)
   print '(A,ES8.1,A,A,A,F9.6,A,I0,A)', 'gap ', gap, ', directional: ', status_name(status), ' at s = ', toi, &
                                      ' (', steps, ' steps)'
enddo
