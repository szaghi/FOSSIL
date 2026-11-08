do q = 2, 4, 2
   call grid(1)%stencils(cell, xi, q, nodes, w, sst)
   e(q / 2) = 0._R8P
   do n = 1, 400 ! interpolate f = sin(x) cos(y) cos(2 pi z), periodic in z as the grid, from the vertices
      f = 0._R8P
      do k = 1, q**3
         i = nodes(1, k, n) ; j = nodes(2, k, n) ; l = nodes(3, k, n)
         f = f + w(k, n) * sin(x(i, j, l)) * cos(y(i, j, l)) * cos(2 * PI * z(i, j, l))
      enddo
      e(q / 2) = max(e(q / 2), abs(f - sin(px(n)) * cos(py(n)) * cos(2 * PI * pz(n))))
   enddo
enddo
print '(A,2ES10.2)', 'largest interpolation error, widths 2 and 4 ', e
