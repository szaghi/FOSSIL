call load_samples_3d(sphere, samples=ls, order=2_I4P, pose=pose)  ! the samples of the body at its pose
p = 1._R8P - ls%points(3, :)                                      ! the solver's pressure at the samples
call integrate_loads(ls, force, moment, pressure=p, about=[0._R8P, 0._R8P, 0._R8P], area=area)
call sphere%mass_properties(volume=vol, center=xc) ; xc = to_world(pose, xc)
print '(A,3F11.7,A,F10.7,A)', 'force ', clean(force), '   buoyancy (0, 0, V), V =', vol, ' (faceted)'
print '(A,3F11.7,A,3F11.7)', 'moment', clean(moment), '   V x_c x (0, 0, 1) =', clean(vol * [xc(2), -xc(1), 0._R8P])
print '(A,F10.7,A,L1)', 'wetted area', area, ';  the identities to 1e-15: ', &
      max(maxval(abs(force - [0._R8P, 0._R8P, vol])), maxval(abs(moment - vol * [xc(2), -xc(1), 0._R8P]))) < 1.e-15_R8P
call sphere%enter_device
call load_samples_3d(sphere, samples=lr, order=2_I4P, pose=pose, resident=.true.)
call integrate_loads(lr, fd, md, pressure=p, about=[0._R8P, 0._R8P, 0._R8P], device=.true.)  ! read on the device
print '(A,L1)', 'device = host, bit for bit: ', all(transfer([fd, md], 0_I8P, 6) == transfer([force, moment], 0_I8P, 6))
call load_samples_reset(lr)
call sphere%exit_device
