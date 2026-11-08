real(R8P),    allocatable    :: w(:,:)
integer(I4P), allocatable    :: nodes(:,:,:), sst(:)
real(R8P)                    :: err
integer(I4P)                 :: q, k
do q = 2, 4
   call ogrid(1)%stencils(cell, xi, q, nodes, w, sst)        ! q x q Lagrange stencils on the vertices
   err = 0._R8P
   do i = 1, N
      if (st(i) /= LOC_INSIDE) cycle
      err = max(err, abs(sum([(w(k, i) * fun(u(nodes(1, k, i), nodes(2, k, i)), v(nodes(1, k, i), nodes(2, k, i))), &
                                k = 1, q**2)]) - fun(pu(i), pv(i))))
   enddo
   print '(A,I0,A,I0,A,ES9.2)', 'q = ', q, ': centred stencils ', count(sst == STENCIL_CENTRED .and. st == LOC_INSIDE), &
         ', max error of sin(2x) cos(y) ', err
enddo
