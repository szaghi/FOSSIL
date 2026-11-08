do j = 1, 21                                 ! a 41 x 21 grid 0.005 above the upper side, cosine chordwise
   do i = 1, 41
      xc = 0.5_R8P * (1._R8P - cos(PI * (i - 1) / 40._R8P))
      yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
      x(i, j, :) = xc - 0.25_R8P ; y(i, j, :) = [yt + 0.005_R8P, yt + 0.1_R8P] ; z(i, j, :) = (j - 1) / 20._R8P
   enddo
enddo
call block%initialize(x, y, z)
call project_block_face(block, face=5_I4P, target=target, tolerance=0.01_R8P, moved=moved, kinds=kinds, status=status)
print '(A,I0,A,I0,A,I0,A,F7.4)', 'face: corners ', count(kinds == PROJ_CORNER), ', on curves ', &
      count(kinds == PROJ_CURVE), ', on patches ', count(kinds == PROJ_PATCH), ', largest move ', maxval(moved)
print '(A,L1,A,L1)', 'trailing-edge row on x = 0.75: ', all(block%x(41, :, 1) == 0.75_R8P), &
                     ', root and tip rows on z = 0, 1: ', all(block%z(:, 1, 1) == 0._R8P) .and. all(block%z(:, 21, 1) == 1._R8P)
