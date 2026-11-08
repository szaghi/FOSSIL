x = [(-0.22_R8P + 0.25_R8P * i, i = 0, 6)]                                     ! 6^3 cells of side 1/4
call swept_cells_3d(cube, p0, pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, &
                    translation=[0.1_R8P, 0.05_R8P, 0._R8P]), x, x, x, sw)
print '(A,L1,A,I0,A,ES9.2,A,ES9.2)', 'translation: exact ', sw%exact, ', evaluations ', sw%evaluations, &
                                     ', residual ', maxval(abs(sw%residual)), ', sum of swept ', sum(sw%swept)
call swept_cells_3d(cube, p0, pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, &
                    translation=[0.1_R8P, 0.05_R8P, 0._R8P]), x, x, x, sw, exact=.false.)
print '(A,L1,A,I0,A,ES9.2,A,ES9.2)', '  Gauss:     exact ', sw%exact, ', nodes ', sw%evaluations, &
                                     ', residual ', maxval(abs(sw%residual)), ', estimate ', sw%error_estimate
call swept_cells_3d(cube, p0, pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0.2_R8P, &
                    translation=[0.1_R8P, 0.05_R8P, 0._R8P]), x, x, x, sw)
print '(A,L1,A,I0,A,ES9.2)', 'rotation:    exact ', sw%exact, ', evaluations ', sw%evaluations, &
                             ', residual ', maxval(abs(sw%residual))
print '(A,I0,A,I0,A,ES10.3)', '  cleared ', count(sw%cleared), ', covered ', count(sw%covered), &
                              ', wall sum ', maxval(abs(sum(sum(sum(sw%wall, 4), 3), 2)))
