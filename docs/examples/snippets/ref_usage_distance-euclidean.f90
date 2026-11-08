p = 2._R8P * ex_R8P + 2._R8P * ey_R8P + 0.5_R8P * ez_R8P
d = surface%distance(point=p, is_square_root=.true.)
print '(A,F10.6)', 'distance of (2, 2, 0.5)               ', d
