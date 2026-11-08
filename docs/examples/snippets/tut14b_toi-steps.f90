call dragon%distance_to(pebble, distance=d, pose_other=p0)                        ! s = 0
print '(A,F9.6,A,ES10.3,A,ES10.3)', 'step 0: s = ', 0._R8P, ', gap ', d, ', next step (d - tol/2) / mu = ', &
                                    (d - 0.5_R8P * tol) / mu
do k = 1, 5                                                                       ! stop after k steps: s_k, d(s_k)
   call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, s, status, distance=d, &
                              max_iterations=k)
   if (d > tol) then
      print '(A,I0,A,F9.6,A,ES10.3,A,ES10.3)', 'step ', k, ': s = ', s, ', gap ', d, ', next step (d - tol/2) / mu = ', &
                                               (d - 0.5_R8P * tol) / mu
   else
      print '(A,I0,A,F9.6,A,ES10.3,A)', 'step ', k, ': s = ', s, ', gap ', d, ' <= tolerance: contact'
   endif
enddo
