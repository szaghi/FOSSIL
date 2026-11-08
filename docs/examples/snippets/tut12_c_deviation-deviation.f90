call mesh%deviation(part, lower=lower, upper=upper, mean=mean, per_vertex=d_vertex, status=status)
if (status /= DEV_STATUS_OK) error stop 'deviation failed'
print '(A,ES9.2,A,ES9.2)', 'h(mesh, part): ', lower, ' <= h <= ', upper     ! how far the mesh strays
print '(A,ES9.2,A,ES9.2)', '  at the vertices ', maxval(d_vertex), ', mean over the mesh ', mean
call part%deviation(mesh, lower=lower, upper=upper, tolerance=5.e-3_R8P, status=status)
if (status /= DEV_STATUS_OK) error stop 'deviation failed'
print '(A,ES9.2,A,ES9.2)', 'h(part, mesh): ', lower, ' <= h <= ', upper     ! what of the part the mesh misses
print '(A,ES9.2)', 'volume change, relative: ', (mesh%get_volume() - part%get_volume()) / part%get_volume()
