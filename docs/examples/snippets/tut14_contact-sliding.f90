! a particle sliding along a flat wall at a gap of 4 tolerances
call wall%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[-0.1_R8P, -0.1_R8P, 0._R8P, 0._R8P], loop_sizes=[4_I4P])
q0 = pose_2d(angle=0._R8P, translation=[0.1_R8P, 0.03_R8P + 4.e-6_R8P])
q1 = pose_2d(angle=0._R8P, translation=[0.9_R8P, 0.03_R8P + 4.e-6_R8P])
call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                         tolerance=1.e-6_R8P, iterations=steps)
print '(A,L1,A,F9.6,A,I0,A)', 'sliding, plain: not converged ', status == CCD_NOT_CONVERGED, ', s >= ', toi, ' (', &
                              steps, ' steps)'
call wall%time_of_impact(disc, pose_2d(angle=0._R8P), pose_2d(angle=0._R8P), q0, q1, toi, status, &
                         tolerance=1.e-6_R8P, iterations=steps, directional=.true.)
print '(A,L1,A,F9.6,A,I0,A)', 'sliding, directional: separate ', status == CCD_SEPARATE, ', s = ', toi, ' (', &
                              steps, ' steps)'
