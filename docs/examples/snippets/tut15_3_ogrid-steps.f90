do k = 10, 40, 10                                                      ! layer k: its height, the distance marched
   print '(A,I2,A,ES10.3,A,ES10.3)', 'layer ', k, ': height ', spacing(k), ', from the wall ', sum(spacing(1:k))
enddo
