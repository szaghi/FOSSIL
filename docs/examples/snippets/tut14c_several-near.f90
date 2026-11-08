do i = 1, 4                                                                      ! the bodies at the first contact
   at_toi(i) = bodies(i)
   call at_toi(i)%move(pose_interpolate(poses0(i), poses1(i), toi))
enddo
call proximity_pairs_of(at_toi, cutoff=0.3_R8P, pairs=near, distances=gaps)
do i = 1, size(gaps)
   print '(A,I0,A,I0,A,ES10.3)', 'within 0.3 at s = toi: bodies ', near(1, i), ' and ', near(2, i), ', gap ', gaps(i)
enddo
