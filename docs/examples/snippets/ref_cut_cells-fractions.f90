c = 0.5_R8P * (g(1:6) + g(2:7))                                       ! the cell centres
call cube%cell_wall_fractions(x=c, y=c, z=c, frac=frac, inside=inside)
print '(A,I0,A,I0)', 'centres inside ', count(inside), ' of ', size(inside)
print '(A,L1,A,6F7.3)', '  centre (2, 3, 3): inside ', inside(2, 3, 3), ', fractions (-x .. +z)', frac(:, 2, 3, 3)
