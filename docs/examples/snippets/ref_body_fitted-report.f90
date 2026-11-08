do i = 1, 129 ! a NACA 0012 section, chord 1, counter-clockwise from the trailing edge
   xc = 0.5_R8P * (1._R8P + cos(2 * PI * (i - 1) / 129._R8P))
   yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
   u(i) = xc ; v(i) = merge(yt, -yt, i <= 65)
enddo
call hyperbolic_grid_2d(u, v, geometric_spacing(2.e-5_R8P, 1.15_R8P, 50_I4P), grid(1), invalid=invalid, status=status)
print '(A,I0,A,I0,A,I0,A,I0)', 'O-grid ', grid(1)%ni, ' x ', grid(1)%nj, ', status ', status, ', not valid ', invalid
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [grid(1)%ni, 1, 1] ; rec%name = 'wall'
call grid(1)%add_record(rec)
call wall_cutter_2d(grid, wall)
print '(A,F8.5)', 'flat-plate u_tau at mid-chord ', utau_flat_plate(U_INF, NU, 0.5_R8P)
call block_report_2d(grid(1), wall, NU, report, u_inf=U_INF, x_le=0._R8P, target_y1p=1._R8P, cell_centred=.true.)
print '(A,I0,A,L1,A,L1)', 'wall cells ', report%nw, ', cell centred ', report%cell_centred, ', u_tau estimated ', &
                          report%utau_estimated
print '(A,3F8.4,A,F6.3)', 'y1+ min mean max ', report%y1p_min, report%y1p_mean, report%y1p_max, &
                          ', fraction above 1 ', report%above
i = maxloc(abs(report%y1_line - report%y1_wall), 1)
print '(A,ES10.3,A,ES10.3,A,I0)', 'y1 along the line ', report%y1_line(i), ', to the wall ', report%y1_wall(i), &
                                  ' at wall cell ', i
print '(A,F7.4,A,F7.4,A,I0)', 'growth max ', report%growth_max, ', skewness max ', report%skew_max, &
                              ', cells not valid ', report%not_valid
call save_block_report(report, 'naca_report.txt')
open(newunit=i, file='naca_report.txt', action='read')
read(i, '(A)') line ; print '(A)', trim(line)
read(i, '(A)') line ; print '(A)', trim(line)
read(i, '(A)') line ; print '(A)', trim(line)
close(i)
