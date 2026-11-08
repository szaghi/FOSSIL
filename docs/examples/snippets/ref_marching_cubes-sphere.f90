previous = 0._R8P
do m = 0, 2
   n = 16 * 2**m + 1                                           ! 17, 33, 65 samples per side of [-1.5, 1.5]^3
   h = 3._R8P / (n - 1)
   allocate(phi(n, n, n))
   do k = 1, n ; do j = 1, n ; do i = 1, n
      phi(i, j, k) = norm2(-1.5_R8P + h * ([i, j, k] - 1)) - 1._R8P   ! signed distance to the unit sphere
   enddo ; enddo ; enddo
   call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                           bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
   if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
   call sphere%adopt_facets(facets=facets)                     ! analyze merges the vertices shared by cells
   error = EXACT - sphere%get_volume()
   print '(A,F7.5,A,I5,A,L1,A,ES8.2)', 'sphere, h ', h, ': ', sphere%get_facets_number(), ' facets, watertight ', &
         sphere%is_watertight(), ', volume error ', error
   if (m > 0) print '(A,F4.2)', '                    error ratio to the previous h ', previous / error
   previous = error ; deallocate(phi)
enddo
