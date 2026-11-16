call load_samples_3d(sphere, samples=ls, order=2_I4P, pose=pose)   ! a wall stress known at the samples (world frame)
allocate(tau(2, ls%n), tv(2, pool%vertex_count()))
tau(1, :) = 0.1_R8P * (1._R8P + ls%points(1, :) - ls%points(3, :))  ! linear
tau(2, :) = 0.1_R8P * cos(4 * ls%points(3, :))                     ! curved
call samples_to_vertices(sphere, ls, tau, tv)                       ! per vertex, area-weighted
