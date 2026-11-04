call bvh%distance_many(points=points, distances=d_wn, is_signed=.true., is_square_root=.true., &
                       sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,I5,A,F10.6)', 'pseudo-normal:  inside ', count(d_bvh < 0._R8P), '  min d ', minval(d_bvh)
print '(A,I5,A,F10.6)', 'winding number: inside ', count(d_wn  < 0._R8P), '  min d ', minval(d_wn)
print '(A,I0,A,ES10.3)', 'signs differ at ', count(d_bvh * d_wn < 0._R8P), ' points, max ||d_pn| - |d_wn|| ', &
                        maxval(abs(abs(d_bvh) - abs(d_wn)))
