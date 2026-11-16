places(1) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [-0.3_R8P, 0._R8P, 0._R8P])
places(2) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.4_R8P, [0.3_R8P, 0._R8P, 0._R8P]) ! they interpenetrate
call balls(1)%enter_device ; call balls(2)%enter_device         ! map the array elements themselves
call fd_closure_3d_bodies(balls, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], table, poses=places, lengths=[0.05_R8P], &
                          ilw=.true., device=.true.)
print '(A,I0,A,I0,A,I0)', 'two balls: ghosts ', table%n, ', points inside both ', table%overlap_points, &
      ', overlapping ghosts ', count(table%overlap)
call fd_closure_3d_bodies(balls, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], twin, poses=places, lengths=[0.05_R8P], ilw=.true.)
print '(A,I0)', '  device path vs host path: mismatches ', ib_closure_mismatches(table, twin)
call balls(1)%exit_device ; call balls(2)%exit_device
