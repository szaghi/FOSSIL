inside = surface%is_point_inside(point=0.5_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0.5_R8P * ez_R8P)
if (inside) print '(A)', 'point (0.5, 0.5, 0.5) is inside the surface'
