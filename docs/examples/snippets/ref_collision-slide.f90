call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)        ! b rests on [0, 1]^3
rest = pose_identity()
a0 = pose_identity() ; a0%t = [-2._R8P, 0._R8P, 1._R8P + 4 * TOL]          ! a slides over b's top face, 4 tol above it
a1 = a0 ; a1%t(1) = 2._R8P
call a%time_of_impact(b, a0, a1, rest, rest, toi, status, tolerance=TOL, distance=d, iterations=steps)
print '(3A,F9.6,A,I0,A,ES10.3)', 'plain:       ', STATUS_NAME(status), ' toi ', toi, ', steps ', steps, &
                                 ', distance ', d
call a%time_of_impact(b, a0, a1, rest, rest, toi, status, tolerance=TOL, distance=d, iterations=steps, &
                      directional=.true., pairs=pairs)
print '(3A,F9.6,A,I0,A,ES10.3,A,I0)', 'directional: ', STATUS_NAME(status), ' toi ', toi, ', steps ', steps, &
                                      ', distance ', d, ', pairs ', pairs
