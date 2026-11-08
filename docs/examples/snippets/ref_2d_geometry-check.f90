defect = c%check_validity(loop=l)
print '(A,I0,A,I0)', 'before repair: defect ', defect, ' at loop ', l
call c%repair(tolerance=1.e-6_R8P, gaps_closed=gaps, loops_reversed=reversed, loops_open=still_open)
print '(A,I0,A,I0,A,I0)', 'repair: gaps closed ', gaps, ', loops reversed ', reversed, ', still open ', still_open
defect = c%check_validity(loop=l, check_self_intersections=.true.)
print '(A,L1,A,I0,A,I0,A,F8.5)', 'after repair: valid ', defect == CURVE_VALID, ', loops ', c%get_loops_number(), &
                                 ', bodies ', c%get_bodies_number(), ', area ', c%get_area()
do l = 1, c%get_loops_number()
   print '(A,I0,A,I0,A,F9.5)', '  loop ', l, ': nesting depth ', c%nesting_depth(l), ', signed area ', c%get_loop_area(l)
enddo
