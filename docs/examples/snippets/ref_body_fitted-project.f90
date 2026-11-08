px = [0.745_R8P, 0.74_R8P, 0.57_R8P, 0.2_R8P]  ! near the trailing-edge corner at the root, near the trailing edge,
py = [0.004_R8P, 0.004_R8P, 0.03_R8P, 0.07_R8P] ! near the upper patch boundary, above the upper surface
pz = [0.01_R8P, 0.5_R8P, 0.4_R8P, 0.5_R8P]
call classify_points(target, px, py, pz, tolerance=0.02_R8P, kinds=kinds, ids=ids)
call project_points(target, px, py, pz, kinds, ids, fx, fy, fz, distances=d)
call project_points(target, px, py, pz, [(PROJ_FREE, i = 1, 4)], [(0_I4P, i = 1, 4)], gx, gy, gz, distances=e)
do i = 1, 4
   print '(A,I0,3A,I0,A,3F9.5,A,F8.5)', 'point ', i, ': ', trim(KIND_NAME(kinds(i) + 1)), ' ', ids(i), ', foot', fx(i), fy(i), &
                                       fz(i), ', moved ', d(i)
   print '(A,3F9.5,A,F8.5)', '         free foot          ', gx(i), gy(i), gz(i), ', moved ', e(i)
enddo
! a point above the main patch, constrained to the flap
call project_points(target, [0.45_R8P], [0.05_R8P], [0.5_R8P], [PROJ_PATCH], [2_I4P], fx(1:1), fy(1:1), fz(1:1), &
                    distances=d(1:1))
print '(A,3F9.5,A,F8.5)', '(0.45, 0.05, 0.5) on the flap ', fx(1), fy(1), fz(1), ', moved ', d(1)
