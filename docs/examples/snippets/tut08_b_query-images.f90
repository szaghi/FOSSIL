allocate(image(3, m), visible(m))
call surface%image_points(points=ghosts, images=image, visible=visible)
print '(A,3F9.5)', '  image      ', image(:, 1)
print '(A,I0,A,I0)', 'images visible from their wall point ', count(visible), ' of ', m
allocate(probe(3, m), length(m))
length = 2._R8P * h                                               ! two spacings out along the normal
call surface%probe_many(points=ghosts, lengths=length, probe_points=probe, visible=visible)
print '(A,3F9.5)', '  probe      ', probe(:, 1)
print '(A,I0,A,I0)', 'probes at 2h visible from their wall point ', count(visible), ' of ', m
