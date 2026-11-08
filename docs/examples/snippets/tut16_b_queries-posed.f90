pose = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
call surface%query_many_device(points=points, distances=d(:, 1), facets=facet(:, 1), pose=pose)  ! chapter 11's step
call surface%query_many(points=points, distances=d(:, 2), facets=facet(:, 2), pose=pose)
print '(A,I0,A,I0)', 'posed query: distance mismatches ', count(d(:, 1) /= d(:, 2)), ', facet mismatches ', &
      count(facet(:, 1) /= facet(:, 2))
call surface%exit_device
