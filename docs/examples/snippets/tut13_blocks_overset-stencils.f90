do q = 2, 4
   call blocks(1)%stencils(cell, xi, q, nodes, w, sst)       ! q x q Lagrange stencils, order-matched xi
   s = 0._R8P
   do i = 1, n
      if (st(i) /= LOC_INSIDE) cycle
      s = max(s, abs(sum([(w(r, i) * f(blocks(1)%u(nodes(1, r, i), nodes(2, r, i)), &
                                       blocks(1)%v(nodes(1, r, i), nodes(2, r, i))), r = 1, q**2)]) - f(pu(i), pv(i))))
   enddo
   print '(A,I0,A,I0,A,ES9.2)', 'q = ', q, ': centred stencils ', count(sst == STENCIL_CENTRED .and. st == LOC_INSIDE), &
         ', max error of sin(2x) cos(y) ', s
enddo
