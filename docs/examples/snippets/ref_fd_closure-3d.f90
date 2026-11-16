do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33                     ! a sphere of radius 0.5 by marching cubes
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call sphere%adopt_facets(facets=facets)
xn = [(-0.6_R8P + 0.05_R8P * i, i = 0, 24)]                  ! node lines, the sphere's poles on them
call ghost_points_3d(sphere, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], gp)
call ghost_points_3d(sphere, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], gd, device=.true.)
print '(A,I0,A,I0)', 'sphere: points ', size(xn)**3, ', ghost points ', gp%n
print '(A,I0,A,I0)', '  first layer ', count([(any(gp%steps(:, g) == 1), g = 1, gp%n)]), &
                     ', second only ', count([(all(gp%steps(:, g) /= 1), g = 1, gp%n)])
print '(A,L1)', '  device path = host path: ', gd%n == gp%n .and. all(gd%ijk == gp%ijk) .and. &
                all(gd%mask == gp%mask) .and. all(gd%steps == gp%steps) .and. all(gd%points == gp%points)
