! unsigned wall distance and body, capped at 0.3
wall = [vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(1.3_R8P, 0.5_R8P, 0.5_R8P), &
        vector_R8P(1.7_R8P, 0.5_R8P, 1.6_R8P)]
call a%wall_distance_many(points=wall, distances=dw, bodies=bw, cutoff=0.3_R8P)
do i = 1, 3
   print '(A,3F5.2,A,F6.3,A,I0)', 'wall distance at', wall(i)%x, wall(i)%y, wall(i)%z, '  ', dw(i), '  body ', bw(i)
enddo
