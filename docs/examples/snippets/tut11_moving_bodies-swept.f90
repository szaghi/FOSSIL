call section%load_from_file(file_name='dragon-section.dat')
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]
q0 = pose_2d(angle=0._R8P)
q1 = pose_2d(angle=0._R8P, translation=[0.005_R8P, 0._R8P])           ! a pure translation, a fifth of a cell
call swept_cells_2d(section, q0, q1, x, x, sw)
print '(A,L1,A,I0,A,I0)', 'section swept: exact ', sw%exact, ', cells with swept area ', count(sw%swept /= 0._R8P), &
                         ', cleared ', count(sw%cleared)
print '(A,L1,A,L1)', '  sum of swept = 0 to 1e-14: ', abs(sum(sw%swept)) < 1.e-14_R8P, &
                     ', residual = 0 to 1e-14: ', maxval(abs(sw%residual)) < 1.e-14_R8P
