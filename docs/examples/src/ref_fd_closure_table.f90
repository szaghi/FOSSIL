!run ref_fd_closure_table ref_fd_closure_table
!< Reference page "Finite-difference wall closure", the closure table: an ellipse (semi-axes 0.4, 0.2, 400 vertices,
!< curvature built) at a pose (turned by 0.3 rad, moved by 0.05 along u) and spinning, on cell-centred lines of spacing
!< 0.02 with reach 2: the table with three probes along the normal (an image line of spacing sqrt(2) h), the wall
!< velocity and the ILW crossings, its first ghost, the device path against the host; then a marching-cubes sphere on
!< node lines, the table left resident on the device and brought back.
program ref_fd_closure_table
use fossil, only : surface_stl_object, curve_object, facet_object, extract_isosurface, ib_closure_t, fd_closure_2d, &
                   fd_closure_3d, ib_closure_mismatches, ib_closure_reset, rigid_pose_2d_t, pose_2d
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P)
type(curve_object)              :: ellipse
type(surface_stl_object)        :: sphere
type(ib_closure_t)              :: table, twin
type(rigid_pose_2d_t)           :: pose
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: h, u(400), v(400), xc(50), yc(30), xn(25), f(33, 33, 33)
integer(I4P)                    :: i, j, k, d, g

u = [(0.4_R8P * cos(2 * PI * (i - 1) / 400), i = 1, 400)] ; v = [(0.2_R8P * sin(2 * PI * (i - 1) / 400), i = 1, 400)]
call ellipse%adopt_loops(u=u, v=v, loop_sizes=[400_I4P])
!region 2d
call ellipse%build_curvature                                      ! the wall curvature, once per geometry
h = 0.02_R8P
xc = [(-0.5_R8P + (i - 0.5_R8P) * h, i = 1, 50)] ; yc = [(-0.3_R8P + (j - 0.5_R8P) * h, j = 1, 30)]
pose = pose_2d(0.3_R8P, [0.05_R8P, 0._R8P])                       ! the body at its pose: turned, moved
call fd_closure_2d(ellipse, xc, yc, [2_I4P, 2_I4P], table, pose=pose, &
                   lengths=[(k * sqrt(2._R8P) * h, k = 1, 3)], &  ! an image line: 3 samples along the normal
                   linear_velocity=[0.5_R8P, 0._R8P], angular_velocity=2._R8P, & ! v = V + omega x (foot - c)
                   ilw=.true.)                                    ! the grid-line crossings of inverse Lax-Wendroff
print '(A,I0,A,I0,A,I0,A,I0)', 'ellipse: ghosts ', table%n, ', probes ', table%n * table%m, ', visible ', &
      count(table%visible), ', ILW crossings ', count(table%ilw_element > 0)
g = 1
print '(A,2I3,A,F9.5,A,2F9.5,A,2F8.4)', 'ghost 1: ij', table%ijk(:, g), '  distance', table%distance(g), &
      '  foot', table%foot(:, g), '  normal', table%normal(:, g)
print '(A,F8.4,A,2F8.4,A,3(2I3,:,","))', '  curvature', table%kappa(1, g), '  velocity', table%velocity(:, g), &
      '  probe cells', table%probe_cell(:, :, g)
do d = 1, 4                                                       ! its ILW slots: -u, +u, -v, +v
   if (table%ilw_element(d, g) > 0) print '(A,I0,A,F9.5,A,I0,A,2F8.4)', '  ILW slot ', d, ': crossing at ', &
      table%ilw_t(d, g), ' on segment ', table%ilw_element(d, g), ', normal', table%ilw_normal(:, d, g)
enddo
call ellipse%enter_device
call fd_closure_2d(ellipse, xc, yc, [2_I4P, 2_I4P], twin, pose=pose, lengths=[(k * sqrt(2._R8P) * h, k = 1, 3)], &
                   linear_velocity=[0.5_R8P, 0._R8P], angular_velocity=2._R8P, ilw=.true., device=.true.)
print '(A,I0)', 'device path vs host path: mismatches ', ib_closure_mismatches(table, twin)
call ellipse%exit_device
!endregion 2d
do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33                         ! a sphere of radius 0.5 by marching cubes
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call sphere%adopt_facets(facets=facets)
xn = [(-0.6_R8P + 0.05_R8P * i, i = 0, 24)]
!region 3d
call sphere%enter_device                                          ! the surface on the device, once
call fd_closure_3d(sphere, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], table, lengths=[0.05_R8P, 0.1_R8P], ilw=.true., &
                   resident=.true.)                               ! the table stays on the device
print '(A,I0,A,L1)', 'sphere: ghosts ', table%n, ', resident ', table%on_device
call table%update_host                                            ! (a device solver would read it there)
call fd_closure_3d(sphere, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], twin, lengths=[0.05_R8P, 0.1_R8P], ilw=.true.)
print '(A,I0,A,I0)', '  ILW crossings ', count(table%ilw_element > 0), ', mismatches against the host path ', &
      ib_closure_mismatches(table, twin)
call ib_closure_reset(table)                                      ! releases the device copies
call sphere%exit_device
!endregion 3d
endprogram ref_fd_closure_table
