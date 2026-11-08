do k = 1, 41 ; do j = 1, 41 ; do i = 1, 41                     ! a unit sphere by marching cubes
   phi(i, j, k) = norm2(-1.5_R8P + 0.075_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call fine%adopt_facets(facets=facets)
call fine%sanitize(fill_holes=.true.)
print '(A,I5,A,F8.5,A)', 'input:          ', fine%get_facets_number(), ' facets, volume ', fine%get_volume(), &
                         ' (sphere: 4.18879)'
