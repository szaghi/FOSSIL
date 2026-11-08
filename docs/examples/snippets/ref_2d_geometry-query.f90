p(:, 1) = [1.5_R8P, 0.5_R8P] ; p(:, 2) = [0.6_R8P, 0.5_R8P] ; p(:, 3) = [2.3_R8P, 0.2_R8P]   ! solid, hole, outside
do k = 1, 3
   call c%compute_distance(p(:, k), distance=d, closest_point=foot, closest_segment=seg, closest_vertex=vtx)
   print '(A,2F5.2,A,F8.5,A,2F8.5,A,I0,A,I0)', 'p = ', p(:, k), ': d ', d, ', foot ', foot, ', segment ', seg, &
                                              ', vertex ', vtx
   print '(A,F8.5,A,F8.5,A,I0,A,L1,A,L1)', '  w fast ', c%winding_number(p(:, k)), ', w sum ', &
         c%winding_number(p(:, k), use_brute_force=.true.), ', w exact ', c%winding_number_exact(p(:, k)), &
         ', inside ', c%is_point_inside(p(:, k)), ', inside (fwn) ', &
         c%is_point_inside(p(:, k), sign_algorithm=SIGN_FAST_WINDING_NUMBER)
enddo
