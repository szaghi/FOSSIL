call surface%image_points(points=ghosts, images=image, visible=visible)
print '(A,3F9.5)', '  image      ', image(:, 1)
print '(A,I0,A,I0)', 'images visible from their wall point ', count(visible), ' of ', m
