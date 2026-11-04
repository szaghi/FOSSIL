allocate(foot(3, m), normal(3, m), image(3, m), feature(m), body(m), ambiguous(m), visible(m))
deallocate(d) ; allocate(d(m))
call surface%query_many(points=ghosts, distances=d, closest_points=foot, features=feature, normals=normal, bodies=body, &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER, ambiguous=ambiguous)
print '(A,I0,A,I0,A,I0)', 'ghost points ', m, ', shells ', surface%get_shells_number(), ', bodies ', surface%get_bodies_number()
print '(A,I0,A,I0,A,I0)', 'feature: face ', count(feature == FEATURE_FACE), ', edge ', count(feature == FEATURE_EDGE), &
                          ', vertex ', count(feature == FEATURE_VERTEX)
print '(A,I0)', 'ambiguous foot points ', count(ambiguous)
print '(A,3F9.5)', 'ghost 1      ', ghosts(1)%x, ghosts(1)%y, ghosts(1)%z
print '(A,3F9.5)', '  wall point ', foot(:, 1)
print '(A,3F9.5)', '  normal     ', normal(:, 1)
