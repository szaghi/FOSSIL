p = 0.5_R8P * ex_R8P + 0.5_R8P * ey_R8P + 0.25_R8P * ez_R8P
sd = surface%distance(point=p, is_signed=.true., is_square_root=.true.)
! Override the sign algorithm for benchmarking or for meshes whose orientation
! cannot be sanitised (e.g. open shells where outward is ill-defined).
sd_ray = surface%distance(point=p, is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_RAY_INTERSECTIONS)
print '(A,F10.6,A,F10.6)', 'signed distance of (0.5, 0.5, 0.25)   ', sd, ', by ray intersections ', sd_ray
