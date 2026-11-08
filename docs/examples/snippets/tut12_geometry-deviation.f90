call mesh%deviation(part, lower=lower, upper=upper, per_vertex=d_vertex, status=status)  ! how far the mesh strays
if (status /= DEV_STATUS_OK) error stop 'deviation failed'
print '(A,ES9.2,A,ES9.2,A,ES9.2)', 'h(mesh, part): ', lower, ' <= h <= ', upper, ', at the vertices ', maxval(d_vertex)
call part%deviation(mesh, lower=lower, upper=upper, tolerance=5.e-3_R8P, status=status) ! what of the part it misses
if (status /= DEV_STATUS_OK) error stop 'deviation failed'
print '(A,ES9.2,A,ES9.2)', 'h(part, mesh): ', lower, ' <= h <= ', upper
print '(A,ES9.2)', 'volume change, relative: ', (mesh%get_volume() - v_part) / v_part
