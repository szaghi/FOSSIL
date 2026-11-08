! a point at mid-gap has two foot points, one on each cube; a point nearer the first cube has one
gap = [vector_R8P(1.1_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(1.05_R8P, 0.5_R8P, 0.5_R8P)]
call a%query_many(points=gap, distances=dg, ambiguous=amb, second_foot_points=second)
do i = 1, 2
   print '(A,F5.2,A,F6.3,A,L1,A,F5.2)', 'gap point x =', gap(i)%x, '  d ', dg(i), '  ambiguous ', amb(i), &
         '  second foot x =', second(1, i)
enddo
