!< Reference page "Finite-difference wall closure", several bodies (#109 WP3 step 01): two circles of radius 0.25 (200
!< vertices) 0.03 apart, on cell-centred lines of spacing 0.02 with reach 2: the table of both, ghosts per body, the
!< shared ghosts of the narrow gap, the probes the other body blocks, one ghost of the gap in full; then two
!< marching-cubes spheres of radius 0.5 that interpenetrate, on the device: the overlap reported, the host path equal.
program ref_fd_closure_bodies
use fossil, only : surface_stl_object, curve_object, facet_object, extract_isosurface, ib_closure_t, fd_closure_2d, &
                   fd_closure_2d_bodies, fd_closure_3d_bodies, ib_closure_mismatches, rigid_pose_2d_t, pose_2d, &
                   rigid_pose_t, pose_from_axis_angle
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P)
type(curve_object)              :: discs(2)
type(surface_stl_object)        :: balls(2)
type(ib_closure_t)              :: table, alone, twin
type(rigid_pose_2d_t)           :: poses(2)
type(rigid_pose_t)              :: places(2)
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: h, u(200), v(200), xc(60), yc(30), xn(41), f(33, 33, 33)
integer(I4P)                    :: i, j, k, g

u = [(0.25_R8P * cos(2 * PI * (i - 1) / 200), i = 1, 200)] ; v = [(0.25_R8P * sin(2 * PI * (i - 1) / 200), i = 1, 200)]
call discs(1)%adopt_loops(u=u, v=v, loop_sizes=[200_I4P]) ; call discs(2)%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
h = 0.02_R8P
xc = [(-0.6_R8P + (i - 0.5_R8P) * h, i = 1, 60)] ; yc = [(-0.3_R8P + (j - 0.5_R8P) * h, j = 1, 30)]
poses = [pose_2d(0._R8P, [-0.265_R8P, 0._R8P]), pose_2d(0._R8P, [0.265_R8P, 0._R8P])] ! a gap of 0.03 = 1.5 h
call fd_closure_2d_bodies(discs, xc, yc, [2_I4P, 2_I4P], table, poses=poses, &
                          lengths=[(k * sqrt(2._R8P) * h, k = 1, 3)], &     ! an image line of 3 probes
                          linear_velocities=reshape([0.1_R8P, 0._R8P, -0.1_R8P, 0._R8P], [2, 2]))
print '(A,I0,A,I0,A,I0)', 'two discs: ghosts ', table%n, ', of body 1 ', count(table%body == 1), ', of body 2 ', &
      count(table%body == 2)
print '(A,F6.3,A,I0,A,F7.4)', '  shared distance ', table%shared_distance, ', shared ghosts ', count(table%shared), &
      ', least clearance ', minval(table%clearance)
call fd_closure_2d(discs(1), xc, yc, [2_I4P, 2_I4P], alone, pose=poses(1), lengths=[(k * sqrt(2._R8P) * h, k = 1, 3)])
print '(A,I0,A,I0)', '  probes blocked: ', count(.not. table%visible), ', disc 1 alone: ', count(.not. alone%visible)
g = minloc(abs(table%points(1, :) + 0.03_R8P) + abs(table%points(2, :) - 0.01_R8P), 1) ! a ghost of the gap
print '(A,2F8.4,A,I0,A,F8.4,A,2F8.4)', 'ghost at', table%points(:, g), ': body ', table%body(g), ', distance', &
      table%distance(g), ', foot', table%foot(:, g)
print '(A,F7.4,A,I0,A,2F8.4,A,L1,A,3L2)', '  clearance', table%clearance(g), ' to body ', table%second_body(g), &
      ' at', table%second_foot(:, g), ', shared ', table%shared(g), ', visible', table%visible(:, g)
do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33                         ! a sphere of radius 0.5 by marching cubes
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call balls(1)%adopt_facets(facets=facets)
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call balls(2)%adopt_facets(facets=facets)
xn = [(-1._R8P + 0.05_R8P * i, i = 0, 40)]
places(1) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [-0.3_R8P, 0._R8P, 0._R8P])
places(2) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.4_R8P, [0.3_R8P, 0._R8P, 0._R8P]) ! they interpenetrate
call balls(1)%enter_device ; call balls(2)%enter_device         ! map the array elements themselves
call fd_closure_3d_bodies(balls, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], table, poses=places, lengths=[0.05_R8P], &
                          ilw=.true., device=.true.)
print '(A,I0,A,I0,A,I0)', 'two balls: ghosts ', table%n, ', points inside both ', table%overlap_points, &
      ', overlapping ghosts ', count(table%overlap)
call fd_closure_3d_bodies(balls, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], twin, poses=places, lengths=[0.05_R8P], ilw=.true.)
print '(A,I0)', '  device path vs host path: mismatches ', ib_closure_mismatches(table, twin)
call balls(1)%exit_device ; call balls(2)%exit_device
endprogram ref_fd_closure_bodies
