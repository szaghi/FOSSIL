call section%load_from_file(file_name='dragon-section.dat')
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]
q0 = pose_2d(angle=0._R8P)
q1 = pose_2d(angle=0._R8P, translation=[0.02_R8P, 0._R8P])            ! a pure translation, 0.8 of a cell
call swept_cells_2d(section, q0, q1, x, x, sw)
print '(A,L1,A,I0,A,I0)', 'section swept: exact ', sw%exact, ', cells with swept area ', count(sw%swept /= 0._R8P), &
                         ', cleared ', count(sw%cleared)
print '(A,L1,A,I0)', '  sum of swept = 0 to 1e-14: ', abs(sum(sw%swept)) < 1.e-14_R8P, &
                    ', cells with residual > 1e-14: ', count(abs(sw%residual) > 1.e-14_R8P)
print '(A,ES9.2,A,2I3,A)', '  largest residual ', maxval(abs(sw%residual)), ' in cell', maxloc(abs(sw%residual)), &
                           ': the knot of the section (chapter 7)'
