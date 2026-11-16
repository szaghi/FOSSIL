!run ref_surface_loads_bodies ref_surface_loads_bodies
!< Reference page "Surface loads", several bodies (#109 WP3 step 02): two circles of radius 0.25 (200 vertices) 0.03
!< apart: the samples of both in body order, the probes the other body blocks, the clearance across the gap, one sample
!< of the gap in full, the buoyancy of each body; then two marching-cubes spheres of radius 0.5, 0.1 apart, on the
!< device: the least clearance, the buoyancy per body, the host path equal.
program ref_surface_loads_bodies
use fossil, only : surface_stl_object, curve_object, facet_object, extract_isosurface, load_samples_t, load_samples_2d, &
                   load_samples_2d_bodies, load_samples_3d_bodies, load_samples_mismatches, integrate_loads, &
                   rigid_pose_2d_t, pose_2d, rigid_pose_t, pose_from_axis_angle
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P)
type(curve_object)              :: discs(2)
type(surface_stl_object)        :: balls(2)
type(load_samples_t)            :: ls, alone, twin
type(rigid_pose_2d_t)           :: poses(2)
type(rigid_pose_t)              :: places(2)
type(facet_object), allocatable :: facets(:)
real(R8P), allocatable          :: p(:), fb(:,:), mb(:)
real(R8P)                       :: h, u(200), v(200), f(33, 33, 33), force(3), moment(3), force2(2), moment2
integer(I4P)                    :: i, j, k, s

u = [(0.25_R8P * cos(2 * PI * (i - 1) / 200), i = 1, 200)] ; v = [(0.25_R8P * sin(2 * PI * (i - 1) / 200), i = 1, 200)]
call discs(1)%adopt_loops(u=u, v=v, loop_sizes=[200_I4P]) ; call discs(2)%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
!region 2d
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
!endregion 2d
do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33                          ! a sphere of radius 0.5 by marching cubes
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call balls(1)%adopt_facets(facets=facets)
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call balls(2)%adopt_facets(facets=facets)
!region 3d
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
!endregion 3d
endprogram ref_surface_loads_bodies
