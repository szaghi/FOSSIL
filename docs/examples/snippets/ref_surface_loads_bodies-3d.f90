places(1) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [-0.55_R8P, 0._R8P, 0._R8P])
places(2) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.4_R8P, [0.55_R8P, 0._R8P, 0._R8P]) ! 0.1 apart
call balls(1)%enter_device ; call balls(2)%enter_device           ! map the array elements themselves
call load_samples_3d_bodies(balls, samples=ls, poses=places, lengths=[0.05_R8P, 0.15_R8P], device=.true.)
print '(A,I0,A,I0,A,F7.4)', 'two balls: samples ', ls%n, ', probes blocked ', count(.not. ls%visible), &
      ', least clearance ', minval(ls%clearance)
p = -ls%points(3, :)                                                ! gravity -z
call integrate_loads(ls, force, moment, pressure=p, per_body_force=fb, device=.true.)
print '(A,F9.5,A,F9.5)', '  buoyancy of ball 2 (z)', fb(3, 2), ', its volume', balls(2)%get_volume()
call load_samples_3d_bodies(balls, samples=twin, poses=places, lengths=[0.05_R8P, 0.15_R8P])
print '(A,I0)', '  device path vs host path: mismatches ', load_samples_mismatches(ls, twin)
call balls(1)%exit_device ; call balls(2)%exit_device
