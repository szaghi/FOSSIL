call slice_surface(wing, normal_axis=3_I4P, level=0.5_R8P, curve=section)   ! the mid-span section, (u, v) = (x, y)
call section%get_loop_range(1_I4P, first, last)
allocate(u(2 * (last - first + 1)), v(2 * (last - first + 1)))
do i = first, last                                                      ! each segment split in 2: the same polygon
   a = section%get_vertex(i) ; b = section%get_vertex(merge(first, i + 1_I4P, i == last))
   do j = 0, 1
      u(2 * (i - first) + j + 1) = a(1) + (b(1) - a(1)) * j / 2._R8P
      v(2 * (i - first) + j + 1) = a(2) + (b(2) - a(2)) * j / 2._R8P
   enddo
enddo
print '(A,I0,A,I0,A)', 'section: ', last - first + 1, ' vertices, ', size(u), ' wall points'
nu = 1.e-6_R8P ; u_inf = 1._R8P                                         ! chord 1, Re = u_inf c / nu = 1e6
u_tau = utau_flat_plate(u_inf, nu, 0.5_R8P)                             ! flat-plate estimate at mid-chord
h1 = 1._R8P * nu / u_tau                                                ! y1+ = 1 by design: y1 = y1+ nu / u_tau
spacing = geometric_spacing(h1, 1.2_R8P, 40_I4P)                        ! 40 layers, ratio 1.2
call hyperbolic_grid_2d(u, v, spacing, grid(1), invalid=invalid, status=status)
print '(A,ES10.3,A,ES10.3,A,F7.3)', 'u_tau ', u_tau, ', h1 ', h1, ', outer height ', sum(spacing)
print '(A,I0,A,I0,A,I0,A,I0)', 'O-grid ', grid(1)%ni, ' x ', grid(1)%nj, ', status ', status, &
                               ', cells not certified valid ', invalid
