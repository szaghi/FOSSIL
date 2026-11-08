k = cut%quad_start(p + 1) - cut%quad_start(p)
print '(A,I0,A,I0,A,I0,A)', 'area rule: ', k, ' points (', &
      count(cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1) < 0._R8P), ' negative), wall rule: ', &
      cut%wall_quad_start(p + 1) - cut%wall_quad_start(p), ' points'
qx = 0._R8P
do k = cut%quad_start(p), cut%quad_start(p + 1) - 1  ! integrate X^2 with the rule
   qx = qx + cut%quad_weight(k) * ((cut%quad_point(1, k) - xc(1)) / h)**2
enddo
print '(A,ES9.2,A,ES9.2)', '  sum of weights - area ', &
      sum(cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1)) - cut%piece_area(p), &
      ', rule on X^2 - moment ', qx - cut%piece_moments(4, p)
