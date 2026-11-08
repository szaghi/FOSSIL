call block_report_2d(grid(1), wall, nu, estimate, u_inf=u_inf, x_le=-0.25_R8P, target_y1p=1._R8P)
print '(A,L1,A,3F8.4,A,F6.3)', 'flat-plate u_tau (estimated ', estimate%utau_estimated, '): y1+ min mean max', &
                               estimate%y1p_min, estimate%y1p_mean, estimate%y1p_max, ', above target ', estimate%above
m = maxloc(estimate%y1p, dim=1)
print '(A,F7.4,A,2F9.5,A)', '  max ', estimate%y1p(m), ' at (x, y) = (', grid(1)%u(estimate%point(1, m), 1), &
                            grid(1)%v(estimate%point(1, m), 1), ')'
