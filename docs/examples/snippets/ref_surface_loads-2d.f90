h = 0.02_R8P
xc = [(-0.5_R8P + (i - 0.5_R8P) * h, i = 1, 50)] ; yc = [(-0.3_R8P + (j - 0.5_R8P) * h, j = 1, 30)]
pose = pose_2d(0.3_R8P, [0.05_R8P, 0._R8P])                       ! the body at its pose: turned, moved
call load_samples_2d(ellipse, xc, yc, samples, order=3_I4P, &     ! three Gauss points per segment
                     lengths=[h, 2 * h], pose=pose)               ! two probes along the normal, + into the fluid
print '(A,I0,A,I0,A,I0,A,F10.7)', 'ellipse: segments ', samples%ne, ', samples ', samples%n, ', visible probes ', &
      count(samples%visible), ', sum of weights ', sum(samples%weight)
print '(A,2F9.5,A,F9.6,A,2F9.5,A,I0)', 'sample 1: point', samples%points(:, 1), '  weight', samples%weight(1), &
      '  normal', samples%normal(:, 1), '  segment ', samples%element(1)
print '(A,2F9.5,A,2F9.5,A,2(2I3,:,","))', '  probes', samples%probe(:, 1, 1), ',', samples%probe(:, 2, 1), &
      '  cells', samples%probe_cell(:, :, 1)
call ellipse%enter_device
call load_samples_2d(ellipse, xc, yc, twin, order=3_I4P, lengths=[h, 2 * h], pose=pose, device=.true.)
print '(A,I0)', 'device path vs host path: mismatches ', load_samples_mismatches(samples, twin)
call ellipse%exit_device
