call section%compute_distance(p=[0._R8P, 0._R8P], distance=d, closest_point=cp, closest_segment=segment, &
                              closest_vertex=vertex)
print '(A,F10.6,A,2F10.6,A,I0)', 'd(0,0) ', d, '  closest ', cp, '  segment ', segment
