allocate(err(n, n, n))
err_band = 0._R8P ; err_far = 0._R8P
do k = 1, n ; do j = 1, n ; do i = 1, n
   d = surface%distance(point=x(i) * ex_R8P + x(j) * ey_R8P + x(k) * ez_R8P, is_signed=.true., is_square_root=.true., &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   err(i, j, k) = abs(field%phi(i, j, k) - d)
   if (field%in_band(i, j, k)) then
      err_band = max(err_band, err(i, j, k))
   else
      err_far = max(err_far, err(i, j, k))
   endif
enddo ; enddo ; enddo
print '(A,ES10.3)', 'max error in the band ', err_band
print '(A,ES10.3)', 'max error elsewhere   ', err_far
