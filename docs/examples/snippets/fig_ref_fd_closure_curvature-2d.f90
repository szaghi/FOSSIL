call ellipse%build_curvature                                      ! the vertex curvatures, once per geometry
call ghost_points_2d(ellipse, xc, yc, [2_I4P, 2_I4P], gp)         ! the ghost points of a stencil of half-width 2
allocate(d(gp%n), foot(2, gp%n), nrm(2, gp%n), kap(1, gp%n), dirs(2, 1, gp%n))
call ellipse%query_many(points=gp%points, distances=d, closest_points=foot, normals=nrm, principal_curvatures=kap, &
                        principal_directions=dirs)                ! kap(1, g): the wall curvature at the foot of ghost g
