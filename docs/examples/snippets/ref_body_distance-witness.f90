call cube(1)%distance_to(cube(3), distance=d, point_self=pa, point_other=pb, facet_self=fa, facet_other=fb)
print '(A,F7.4,A,I0,A,I0)', 'cubes 1 and 3: ', d, ', facets ', fa, ' and ', fb
print '(A,3F6.2,A,3F6.2)', '  from', pa, ' to', pb
