h = 0.02_R8P
xc = [(-0.6_R8P + (i - 0.5_R8P) * h, i = 1, 60)] ; yc = [(-0.3_R8P + (j - 0.5_R8P) * h, j = 1, 30)]
poses = [pose_2d(0._R8P, [-0.265_R8P, 0._R8P]), pose_2d(0._R8P, [0.265_R8P, 0._R8P])] ! a gap of 0.03 = 1.5 h
call fd_closure_2d_bodies(discs, xc, yc, [2_I4P, 2_I4P], table, poses=poses, &
                          lengths=[(k * sqrt(2._R8P) * h, k = 1, 3)], &     ! an image line of 3 probes
                          linear_velocities=reshape([0.1_R8P, 0._R8P, -0.1_R8P, 0._R8P], [2, 2]))
print '(A,I0,A,I0,A,I0)', 'two discs: ghosts ', table%n, ', of body 1 ', count(table%body == 1), ', of body 2 ', &
      count(table%body == 2)
print '(A,F6.3,A,I0,A,F7.4)', '  shared distance ', table%shared_distance, ', shared ghosts ', count(table%shared), &
      ', least clearance ', minval(table%clearance)
call fd_closure_2d(discs(1), xc, yc, [2_I4P, 2_I4P], alone, pose=poses(1), lengths=[(k * sqrt(2._R8P) * h, k = 1, 3)])
print '(A,I0,A,I0)', '  probes blocked: ', count(.not. table%visible), ', disc 1 alone: ', count(.not. alone%visible)
g = minloc(abs(table%points(1, :) + 0.03_R8P) + abs(table%points(2, :) - 0.01_R8P), 1) ! a ghost of the gap
print '(A,2F8.4,A,I0,A,F8.4,A,2F8.4)', 'ghost at', table%points(:, g), ': body ', table%body(g), ', distance', &
      table%distance(g), ', foot', table%foot(:, g)
print '(A,F7.4,A,I0,A,2F8.4,A,L1,A,3L2)', '  clearance', table%clearance(g), ' to body ', table%second_body(g), &
      ' at', table%second_foot(:, g), ', shared ', table%shared(g), ', visible', table%visible(:, g)
