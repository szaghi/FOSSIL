do v = 1, 3
   call make_wedge(f, v)
   call wedge%adopt_facets(facets=f)
   call wedge%analyze
   do a = 1, 4
      wrong(a) = 0
      do i = 1, size(p)
         if (wedge%is_point_inside(point=p(i), sign_algorithm=ALGO(a)) .neqv. truth(i)) wrong(a) = wrong(a) + 1
      enddo
   enddo
   wrong(5) = 0
   do i = 1, size(p)
      if (naive_inside(p(i)) .neqv. truth(i)) wrong(5) = wrong(5) + 1
   enddo
   print '(A,I13,I5,I13,I9,I14)', VARIANT(v), wrong
   if (v == 1) then                                ! the naive errors by decade, delta = 1e-1 first
      do i = 1, ND
         per_decade(i) = count([(naive_inside(p(a)) .neqv. truth(a), a = 17 * (i - 1) + 1, 17 * i)])
      enddo
   endif
enddo
print '(A,12I3)', 'naive errors on the clean wedge by decade (1e-1 to 1e-12):', per_decade
