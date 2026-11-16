call sphere%enter_device                                          ! the surface on the device, once
call load_samples_3d(sphere, xn, xn, xn, samples, order=3_I4P, lengths=[0.05_R8P, 0.1_R8P], &
                     resident=.true.)                             ! the samples stay on the device
print '(A,I0,A,I0,A,L1)', 'sphere: facets ', samples%ne, ', samples ', samples%n, ', resident ', samples%on_device
call samples%update_host                                          ! (a device solver would read them there)
call load_samples_3d(sphere, xn, xn, xn, twin, order=3_I4P, lengths=[0.05_R8P, 0.1_R8P])
print '(A,F9.6,A,F9.6,A,I0)', '  area ', sum(samples%weight), ' (4 pi r^2 =', PI, '), mismatches against the host path ', &
      load_samples_mismatches(samples, twin)
call load_samples_reset(samples)                                  ! releases the device copies
call sphere%exit_device
