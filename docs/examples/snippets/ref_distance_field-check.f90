allocate(err(n, n, n))
band_err = 0._R8P ; far_err = 0._R8P ; foot_err = 0._R8P
do k = 1, n ; do j = 1, n ; do i = 1, n
   e = abs([x(i), x(j), x(k)] - 0.5_R8P) - 0.5_R8P                        ! exact distance of the box
   exact = norm2(max(e, 0._R8P)) + min(maxval(e), 0._R8P)
   err(i, j, k) = abs(field%phi(i, j, k) - exact)
   if (field%in_band(i, j, k)) then
      band_err = max(band_err, err(i, j, k))
      foot_err = max(foot_err, abs(norm2([x(i), x(j), x(k)] - field%closest(:, i, j, k)) - abs(field%phi(i, j, k))))
   else
      far_err = max(far_err, err(i, j, k))
   endif
enddo ; enddo ; enddo
print '(A,ES9.2,A,F6.3,A,F6.3)', 'max |phi - exact|: band ', band_err, ', far field ', far_err, &
                                 ', largest spacing ', maxval(x(2:) - x(:n - 1))
print '(A,ES9.2)', 'band: max ||x - closest| - |phi|| ', foot_err
