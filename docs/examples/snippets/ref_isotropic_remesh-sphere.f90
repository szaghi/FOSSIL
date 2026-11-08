do k = 1, 25 ; do j = 1, 25 ; do i = 1, 25                       ! a unit sphere by marching cubes
   phi(i, j, k) = norm2(-1.5_R8P + 0.125_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call input%adopt_facets(facets=facets)
call input%sanitize(fill_holes=.true.)
print '(A,I4,A,F6.4,A,F6.4,A,F7.5)', 'sphere input:  ', input%get_facets_number(), ' facets, edges ', &
      input%smallest_edge_len(), ' to ', input%largest_edge_len(), ', volume ', input%get_volume()
sphere = input
do k = 1, 3                                                       ! one iteration per call, to watch it
   call sphere%isotropic_remesh(target_length=0.1_R8P, iterations=1_I4P, status=status)
   print '(A,I1,A,I4,A,F6.4,A,F6.4,A,F7.5)', 'iteration ', k, ':   ', sphere%get_facets_number(), &
         ' facets, edges ', sphere%smallest_edge_len(), ' to ', sphere%largest_edge_len(), ', volume ', sphere%get_volume()
enddo
call sphere%hausdorff(other=input, lower=lower, upper=upper)     ! certified distance to the input
print '(A,F6.4,A,F6.4)', 'Hausdorff distance to the input ', lower, ' .. ', upper
