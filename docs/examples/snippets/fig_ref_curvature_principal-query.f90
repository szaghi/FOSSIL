call torus%build_curvature                    ! the per-vertex shape operator (cubic jet over the 2-ring)
call torus%query_many(points=q, distances=d, closest_points=cp, normals=nm, principal_curvatures=k, &
                      principal_directions=dirs)
