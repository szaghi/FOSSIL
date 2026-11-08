allocate(points(n**3), d_pn(n**3), d_fwn(n**3), d_host(n**3))
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n                     ! x fastest, z slowest: slice k is contiguous
   p = p + 1
   points(p) = (-0.6_R8P + 1.2_R8P * (i - 1) / (n - 1)) * ex_R8P + (-0.6_R8P + 1.2_R8P * (j - 1) / (n - 1)) * ey_R8P + &
               (-0.6_R8P + 1.2_R8P * (k - 1) / (n - 1)) * ez_R8P
enddo ; enddo ; enddo
call surface%distance_many_device(points=points, distances=d_pn, is_signed=.true., is_square_root=.true., &
                                  sign_algorithm=SIGN_PSEUDO_NORMAL, status=status)
call surface%distance_many(points=points, distances=d_host, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_PSEUDO_NORMAL)
print '(A,I0,A,I0)', 'pseudo-normal : device vs host mismatches ', count(d_pn /= d_host), ' of ', size(points)
call surface%distance_many_device(points=points, distances=d_fwn, is_signed=.true., is_square_root=.true., &
                                  sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
call surface%distance_many(points=points, distances=d_host, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,I0,A,I0)', 'winding number: device vs host mismatches ', count(d_fwn /= d_host), ' of ', size(points)
print '(A,I0,A,I0)', 'the two signs disagree at ', count(sign(1._R8P, d_pn) /= sign(1._R8P, d_fwn)), &
                     ' points; |d| differs at ', count(abs(d_pn) /= abs(d_fwn))
