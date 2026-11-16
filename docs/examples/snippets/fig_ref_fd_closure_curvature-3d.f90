call spheroid%build_curvature                                     ! before enter_device: mapped with the surface
call spheroid%query_many(points=q, distances=d, closest_points=foot, principal_curvatures=kap, principal_directions=dirs)
do i = 1, gp%n
   vt = [-foot(2, i), foot(1, i), 0._R8P]                          ! the wall velocity of a spin about z, Omega = 1
   if (norm2(vt) > 0._R8P) kn = normal_curvature(kap(:, i), dirs(:, :, i), vt) ! then d_n v_t = kn v_t along v_t
enddo
