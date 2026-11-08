rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [grid(1)%ni, 1, 1] ; rec%name = 'wall'
call grid(1)%add_record(rec)                                            ! side j = 1 is the wall
call wall_cutter_2d(grid, wall)                                         ! the wall as a curve
call block_report_2d(grid(1), wall, nu, report, u_tau=u_tau, target_y1p=1._R8P)
print '(A,I0,A,3F8.4,A,F6.3)', 'u_tau given, ', report%nw, ' wall points: y1+ min mean max', report%y1p_min, &
                               report%y1p_mean, report%y1p_max, ', above target ', report%above
m = maxloc(report%y1p, dim=1)
print '(A,F7.4,A,I0,A,2F9.5,A)', '  max ', report%y1p(m), ' at wall point ', report%point(1, m), ', (x, y) = (', &
                                 grid(1)%u(report%point(1, m), 1), grid(1)%v(report%point(1, m), 1), ')'
print '(A,ES10.3,A,ES10.3)', '  y1_line there ', report%y1_line(m), ', y1_wall ', report%y1_wall(m)
print '(A,F7.4,A,F7.4,A,F7.4,A,I0)', '  growth max ', report%growth_max, ', skewness max ', report%skew_max, &
                                     ', orthogonality min ', report%orth_min, ', not valid ', report%not_valid
call grid(1)%quality(orthogonality=orth)                                ! the cells behind the minimum (chapter 13)
print '(A,I0,A,I0,A,2I3)', '  cells with orthogonality < 0.5: ', count(orth < 0.5_R8P), ' of ', size(orth), &
                           ', the worst (i, j) =', minloc(orth)
