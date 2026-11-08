px = 0.75_R8P ; py = 0.00626_R8P ; pz = [0.1_R8P, 0.5_R8P, 0.9_R8P]  ! 0.005 above the upper trailing-edge line
call classify_points(target, px, py, pz, tolerance=0.02_R8P, kinds=pk, ids=pid)
print '(A,3I2)', 'classified kinds (2 = PROJ_CURVE): ', pk
call project_points(target, px, py, pz, [PROJ_FREE, PROJ_FREE, PROJ_FREE], [0_I4P, 0_I4P, 0_I4P], fx, fy, fz)
print '(A,3(" (",F7.4,",",F8.5,")"))', 'free feet x, y:        ', (fx(i), fy(i), i = 1, 3)
call project_points(target, px, py, pz, pk, pid, fx, fy, fz)
print '(A,3(" (",F7.4,",",F8.5,")"))', 'constrained feet x, y: ', (fx(i), fy(i), i = 1, 3)
