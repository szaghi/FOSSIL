do k = 1, 41 ; do j = 1, 41 ; do i = 1, 41                     ! a unit sphere by marching cubes
   phi(i, j, k) = norm2(-1.5_R8P + 0.075_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call noisy%adopt_facets(facets=facets)
call noisy%sanitize(fill_holes=.true.)
allocate(facets(noisy%get_facets_number()))                    ! every vertex moved by Gaussian noise, sigma 0.035
do i = 1, size(facets)
   f => noisy%facet_at(i)
   facets(i) = f
   do j = 1, 3
      facets(i)%vertex(j) = facets(i)%vertex(j) + 0.035_R8P * noise(facets(i)%vertex(j))
   enddo
   call facets(i)%compute_normal
enddo
call noisy%adopt_facets(facets=facets)
call report('noisy           ', noisy)
