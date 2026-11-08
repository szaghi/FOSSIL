call surface%wall_distance_many_device(points=points, distances=d(:, 1), bodies=body(:, 1), cutoff=0.05_R8P)
call surface%wall_distance_many(points=points, distances=d(:, 2), bodies=body(:, 2), cutoff=0.05_R8P)
print '(A,I0,A,I0,A,I0)', 'wall distance, cutoff 0.05: mismatches ', count(d(:, 1) /= d(:, 2)), &
      ', probes nearer than the cutoff ', count(body(:, 1) /= 0), ' of ', size(points)
