call ellipse%build_curvature                                      ! the wall curvature, once per geometry
h = 0.02_R8P
xc = [(-0.5_R8P + (i - 0.5_R8P) * h, i = 1, 50)] ; yc = [(-0.3_R8P + (j - 0.5_R8P) * h, j = 1, 30)]
pose = pose_2d(0.3_R8P, [0.05_R8P, 0._R8P])                       ! the body at its pose: turned, moved
call fd_closure_2d(ellipse, xc, yc, [2_I4P, 2_I4P], table, pose=pose, &
                   lengths=[(k * sqrt(2._R8P) * h, k = 1, 3)], &  ! an image line: 3 samples along the normal
                   linear_velocity=[0.5_R8P, 0._R8P], angular_velocity=2._R8P, & ! v = V + omega x (foot - c)
                   ilw=.true.)                                    ! the grid-line crossings of inverse Lax-Wendroff
print '(A,I0,A,I0,A,I0,A,I0)', 'ellipse: ghosts ', table%n, ', probes ', table%n * table%m, ', visible ', &
      count(table%visible), ', ILW crossings ', count(table%ilw_element > 0)
g = 1
print '(A,2I3,A,F9.5,A,2F9.5,A,2F8.4)', 'ghost 1: ij', table%ijk(:, g), '  distance', table%distance(g), &
      '  foot', table%foot(:, g), '  normal', table%normal(:, g)
print '(A,F8.4,A,2F8.4,A,3(2I3,:,","))', '  curvature', table%kappa(1, g), '  velocity', table%velocity(:, g), &
      '  probe cells', table%probe_cell(:, :, g)
do d = 1, 4                                                       ! its ILW slots: -u, +u, -v, +v
   if (table%ilw_element(d, g) > 0) print '(A,I0,A,F9.5,A,I0,A,2F8.4)', '  ILW slot ', d, ': crossing at ', &
      table%ilw_t(d, g), ' on segment ', table%ilw_element(d, g), ', normal', table%ilw_normal(:, d, g)
enddo
call ellipse%enter_device
call fd_closure_2d(ellipse, xc, yc, [2_I4P, 2_I4P], twin, pose=pose, lengths=[(k * sqrt(2._R8P) * h, k = 1, 3)], &
                   linear_velocity=[0.5_R8P, 0._R8P], angular_velocity=2._R8P, ilw=.true., device=.true.)
print '(A,I0)', 'device path vs host path: mismatches ', ib_closure_mismatches(table, twin)
call ellipse%exit_device
