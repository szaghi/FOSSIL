call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)        ! [0, 1]^3: the frame origin at a corner
allocate(x(0), y(0), z(0))
do i = 1, a%get_facets_number()
   f => a%facet_at(i)
   x = [x, (f%vertex(k)%x, k = 1, 3)] ; y = [y, (f%vertex(k)%y, k = 1, 3)] ; z = [z, (f%vertex(k)%z, k = 1, 3)]
enddo
a0 = pose_identity()
a1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI / 2, translation=[1._R8P, 0._R8P, 0._R8P])
print '(A,F8.5)', 'mu, frame origin at a corner:   ', motion_bound_3d(a0, a1, x, y, z)
print '(A,F8.5)', 'mu, frame origin at the centre: ', motion_bound_3d(a0, a1, x - 0.5_R8P, y - 0.5_R8P, z - 0.5_R8P)
