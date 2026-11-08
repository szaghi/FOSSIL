call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, toi, status, distance=d, point_self=pa, &
                           point_other=pb, facet_self=fa, facet_other=fb, iterations=steps)
print '(A,L1,A,F9.6,A,I0,A)', 'contact ', status == CCD_CONTACT, ' at s = ', toi, ' (', steps, ' steps)'
print '(A,ES10.3,A,I0,A,I0)', '  gap there ', d, ', dragon facet ', fa, ', pebble facet ', fb
print '(A,3F8.4,A,3F8.4)', '  on the dragon', pa, ', on the pebble', pb
call dragon%time_of_impact(pebble, pose_identity(), pose_identity(), p0, p1, toi_dir, status, iterations=steps, &
                           directional=.true., pairs=pairs)
print '(A,F9.6,A,I0,A,I0,A,L1)', 'directional: s = ', toi_dir, ' (', steps, ' steps, ', pairs, &
                                 ' exact pair evaluations), same contact to 1e-6: ', abs(toi_dir - toi) < 1.e-6_R8P
