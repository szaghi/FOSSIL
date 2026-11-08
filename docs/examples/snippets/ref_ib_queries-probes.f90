! three ghost points inside the first cube: near the face x = 1, deeper, near the top face
ghost = [vector_R8P(0.95_R8P, 0.3_R8P, 0.4_R8P), vector_R8P(0.75_R8P, 0.3_R8P, 0.4_R8P), &
         vector_R8P(0.9_R8P, 0.3_R8P, 0.97_R8P)]
call a%query_many(points=ghost, distances=d, bodies=bod)
call a%probe_many(points=ghost, lengths=[0.3_R8P, 0.3_R8P, 0.3_R8P], probe_points=probe, visible=vis)
do i = 1, 3
   print '(A,3F5.2,A,F6.3,A,I0,A,3F6.3,A,L1)', 'ghost', ghost(i)%x, ghost(i)%y, ghost(i)%z, '  d ', d(i), &
         '  body ', bod(i), '  probe', probe(:, i), '  visible ', vis(i)
enddo
call a%image_points(points=ghost, images=img, visible=vis)
do i = 1, 3
   print '(A,I0,A,3F6.3,A,L1)', 'image ', i, ':', img(:, i), '  visible ', vis(i)
enddo
