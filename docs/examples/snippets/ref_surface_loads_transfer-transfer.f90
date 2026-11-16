call fd_closure_3d(sphere, x + 0.3_R8P, x + 0.2_R8P, x + 0.1_R8P, [2_I4P, 2_I4P, 2_I4P], cl, pose=pose, ilw=.true.)
allocate(tf(1, cl%n), ti(1, 6, cl%n))
call transfer_to_closure(sphere, cl, tf, vertex_field=tw, ilw_values=ti, pose=pose)   ! T at the feet, the crossings
