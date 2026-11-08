p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
call wrap%move(p1)
call wrap%mass_properties(inertia=inertia1, density=1000._R8P)
print '(A,ES9.2)', 'moved wrap against inertia_to_world, relative ', &
                   maxval(abs(inertia1 - inertia_to_world(inertia, p1))) / maxval(abs(inertia))
