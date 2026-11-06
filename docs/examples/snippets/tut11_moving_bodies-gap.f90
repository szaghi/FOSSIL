call box%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! a unit cube beside the dragon
call box%translate(x=0.7_R8P, y=-0.5_R8P, z=-0.5_R8P)
call box%analyze
call surface%distance_to(box, distance=gap, facet_self=fa, facet_other=fb, intersecting=hit, pose_self=p1)
print '(A,F8.5,A,I0,A,I0,A,L1)', 'gap to the cube ', gap, ' (facets ', fa, ', ', fb, '), touching ', hit
call moved%distance_to(box, distance=gap1)                              ! the moved copy of the queries
print '(A,L1)', 'posed = moved, bit for bit: ', gap == gap1
