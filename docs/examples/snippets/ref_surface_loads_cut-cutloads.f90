call cut_cells_3d(sphere, x + 0.1_R8P, x, x - 0.1_R8P, cut, walls=walls)   ! the cut and its wall terms
p = 1._R8P - walls%centroid(3, :)                                 ! the pressure at the walls' centroids
call cut_loads_3d(walls, force, m0, pressure=p, about=[0._R8P, 0._R8P, 0._R8P], area=area)
call cut_walls_gyration(walls, sphere, x + 0.1_R8P, x, x - 0.1_R8P, gyration)   ! for the exact moment
g = spread([0._R8P, 0._R8P, -1._R8P], 2, walls%n)                 ! the pressure gradient on every wall
call cut_loads_3d(walls, force, m1, pressure=p, gradient=g, gyration=gyration, about=[0._R8P, 0._R8P, 0._R8P])
