h = 0.02_R8P
poses = [pose_2d(0._R8P, [-0.265_R8P, 0._R8P]), pose_2d(0._R8P, [0.265_R8P, 0._R8P])] ! a gap of 0.03
call load_samples_2d_bodies(discs, samples=ls, poses=poses, lengths=[h, 2 * h])   ! order 2, two probes
print '(A,I0,A,I0,A,I0)', 'two discs: samples ', ls%n, ', of body 1 ', count(ls%body == 1), ', of body 2 ', &
      count(ls%body == 2)
call load_samples_2d(discs(1), samples=alone, pose=poses(1), lengths=[h, 2 * h])
print '(A,I0,A,I0)', '  probes blocked: ', count(.not. ls%visible), ', disc 1 alone: ', count(.not. alone%visible)
print '(A,F7.4,A,I0,A,I0)', '  least clearance ', minval(ls%clearance), ', samples within 2 h of the other body ', &
      count(ls%clearance < 2 * h), ', second body of sample 1: ', ls%second_body(1)
s = minloc(ls%clearance, 1)                                         ! a sample of the gap
print '(A,2F8.4,A,I0,A,F7.4,A,2F8.4,A,2L2)', 'sample at', ls%points(:, s), ': body ', ls%body(s), ', clearance', &
      ls%clearance(s), ' to', ls%second_foot(:, s), ', visible', ls%visible(:, s)
p = -ls%points(2, :)                                                ! hydrostatic p = -v (rho g = 1, gravity -v)
call integrate_loads(ls, force2, moment2, pressure=p, per_body_force=fb, per_body_moment=mb)
print '(A,2F9.5,A,F9.5)', '  buoyancy per body (v)', fb(2, :), ', area of a disc', discs(1)%get_area()
