b = findloc(codes, BOX_CUT, dim=1)                                             ! the first cut block
print '(A,I0,A,L1)', 'block ', b, ' meets the dragon: ', surface%box_intersects_surface(lo(:, b), hi(:, b))
print '(A,L1)', 'its lower half meets it:  ', &
      surface%box_intersects_surface(lo(:, b), [hi(1, b), hi(2, b), lo(3, b) + 0.5_R8P * hb])
