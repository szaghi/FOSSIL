px = [0.735_R8P, 0.75_R8P, 0.757_R8P, 0.758_R8P, 0.762_R8P, 0.756_R8P]                  ! around the trailing edge, mid-span
py = [0.0075_R8P, 0.0063_R8P, 0.0002_R8P, 0.0040_R8P, -0.0005_R8P, -0.005_R8P]
pz = 0.5_R8P
call classify_points(target, px, py, pz, tolerance=0.01_R8P, kinds=kinds, ids=ids)
free = PROJ_FREE
call project_points(target, px, py, pz, free, ids, ffx, ffy, ffz)   ! anywhere on the wing
call project_points(target, px, py, pz, kinds, ids, fx, fy, fz)     ! on the constraint found
do i = 1, N
   print '(A,I0,A,F6.3,A,F7.4,A,A12,I3,A,F7.4,A,F8.5,A,F7.4,A,F8.5,A)', 'point ', i, ' (', px(i), ',', py(i), '): ', &
         NAME(kinds(i)), ids(i), '  free (', ffx(i), ',', ffy(i), ')  constrained (', fx(i), ',', fy(i), ')'
enddo
