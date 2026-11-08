n = 3 * pebble%get_facets_number()
allocate(x(n), y(n), z(n))
do i = 1, pebble%get_facets_number()                                              ! the pebble's vertices (body frame)
   f => pebble%facet_at(i)
   x(3 * i - 2:3 * i) = f%vertex%x ; y(3 * i - 2:3 * i) = f%vertex%y ; z(3 * i - 2:3 * i) = f%vertex%z
enddo
mu = motion_bound_3d(p0, p1, x, y, z)                                             ! |dt| + |theta| max rho
da = dragon%get_bmax() - dragon%get_bmin() ; db = pebble%get_bmax() - pebble%get_bmin()
tol = 1.e-6_R8P * max(da%normL2(), db%normL2())                                   ! the default tolerance
print '(A,F7.4,A,F7.4,A,ES10.3)', 'pebble speed bound mu = ', mu, ' (translation 1.2, spin ', mu - 1.2_R8P, &
                                  '), tolerance ', tol
