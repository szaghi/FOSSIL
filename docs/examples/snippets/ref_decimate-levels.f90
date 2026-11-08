coarse = fine
target = fine%get_facets_number()
do k = 1, 3
   target = target / 4
   call coarse%decimate(target_facets=target, status=status)    ! in place: each level from the previous one
   if (status /= DEC_STATUS_OK) error stop 'decimate did not reach the target'
   call coarse%hausdorff(other=fine, lower=lower, upper=upper)    ! certified distance to the input
   print '(A,I5,A,I5,A,F8.5,A,L1,A,F6.4,A,F6.4)', 'target ', target, ': ', coarse%get_facets_number(), &
         ' facets, volume ', coarse%get_volume(), ', closed manifold ', coarse%is_watertight() .and. &
         coarse%is_manifold(), ', Hausdorff ', lower, ' .. ', upper
enddo
