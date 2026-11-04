p = 0._R8P * ex_R8P + 0.6_R8P * ey_R8P + 0._R8P * ez_R8P
call surface%compute_distance(point=p, distance=d, is_signed=.true., is_square_root=.true., &
                              facet_index=facet, edge_index=edge, vertex_index=vertex, closest_point=closest, normal=normal)
print '(A,F10.6,A,I0,A,I0,A,I0)', 'd ', d, '  facet ', facet, '  edge ', edge, '  vertex ', vertex
print '(A,3F10.6)', 'closest point ', closest%x, closest%y, closest%z
print '(A,3F10.6)', 'normal        ', normal%x, normal%y, normal%z
