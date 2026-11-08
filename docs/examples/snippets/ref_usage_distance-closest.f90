call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
p = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 25._R8P * ez_R8P
call dragon%compute_distance(point=p, distance=dist, is_square_root=.true., facet_index=facet, edge_index=edge, &
                             vertex_index=vertex, closest_point=closest, normal=normal)
print '(A,F9.4)',   'dragon, (0, 0, 25): distance ', dist
print '(A,3F9.4)',  '  closest point ', closest%x, closest%y, closest%z
print '(A,3F9.4)',  '  normal        ', normal%x, normal%y, normal%z
print '(A,I0,A,I0,A,I0)', '  facet ', facet, ', local edge ', edge, ', local vertex ', vertex
