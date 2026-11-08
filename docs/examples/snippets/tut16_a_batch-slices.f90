do k = 1, n                                                 ! the same points, one host call per slice z = const
   s = (k - 1) * n**2
   call surface%distance_many(points=points(s + 1:s + n**2), distances=d_host(s + 1:s + n**2), is_signed=.true., &
                              is_square_root=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER)
enddo
print '(A,I0,A,I0)', '24 host calls of one slice vs one device batch: mismatches ', count(d_host /= d_fwn), &
                     ' of ', size(points)
