rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [grid(1)%ni, 1, 1] ; rec%name = 'wall'
call grid(1)%add_record(rec)                                            ! side j = 1 is the wall
call wall_cutter_2d(grid, wall)                                         ! the wall as a curve
call block_report_2d(grid(1), wall, nu, report, u_tau=u_tau, target_y1p=1._R8P)
print '(A,3F8.4,A,F6.3)', 'u_tau given:      y1+ min mean max', report%y1p_min, report%y1p_mean, report%y1p_max, &
                          ', above target ', report%above
m = maxloc(report%y1p, dim=1)
print '(A,F7.4,A,I0,A,2F9.5,A)', '  max ', report%y1p(m), ' at wall point ', report%point(1, m), ', (x, y) = (', &
                                 grid(1)%u(report%point(1, m), 1), grid(1)%v(report%point(1, m), 1), ')'
call block_report_2d(grid(1), wall, nu, report, u_inf=u_inf, x_le=-0.25_R8P, target_y1p=1._R8P)
print '(A,3F8.4,A,F6.3)', 'flat-plate u_tau: y1+ min mean max', report%y1p_min, report%y1p_mean, report%y1p_max, &
                          ', above target ', report%above
m = maxloc(report%y1p, dim=1)
print '(A,F7.4,A,I0,A,2F9.5,A)', '  max ', report%y1p(m), ' at wall point ', report%point(1, m), ', (x, y) = (', &
                                 grid(1)%u(report%point(1, m), 1), grid(1)%v(report%point(1, m), 1), ')'
print '(A,F7.4,A,F7.4,A,F7.4,A,I0)', 'growth max ', report%growth_max, ', skewness max ', report%skew_max, &
                                     ', orthogonality min ', report%orth_min, ', not valid ', report%not_valid
call save_block_report(report, 'tut15-report.txt')
open(newunit=unit, file='tut15-report.txt', action='read')
do i = 1, 7
   read(unit, '(A)') line
   if (i /= 5) print '(A)', trim(line)                                  ! line 5, the cell quality, is long
enddo
close(unit)
