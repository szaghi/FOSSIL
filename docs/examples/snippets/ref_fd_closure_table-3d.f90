call sphere%enter_device                                          ! the surface on the device, once
call fd_closure_3d(sphere, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], table, lengths=[0.05_R8P, 0.1_R8P], ilw=.true., &
                   resident=.true.)                               ! the table stays on the device
print '(A,I0,A,L1)', 'sphere: ghosts ', table%n, ', resident ', table%on_device
call table%update_host                                            ! (a device solver would read it there)
call fd_closure_3d(sphere, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], twin, lengths=[0.05_R8P, 0.1_R8P], ilw=.true.)
print '(A,I0,A,I0)', '  ILW crossings ', count(table%ilw_element > 0), ', mismatches against the host path ', &
      ib_closure_mismatches(table, twin)
call ib_closure_reset(table)                                      ! releases the device copies
call sphere%exit_device
