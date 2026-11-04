half = surface
call half%clip(bmin=bmin, bmax=bmax%x * ex_R8P + 0._R8P * ey_R8P + bmax%z * ez_R8P, remainder=rest)
print '(A,I0,A,I0)', 'clip y<=0 ', half%get_facets_number(), ' facets, remainder ', rest%get_facets_number()
