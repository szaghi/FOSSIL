!< Reference page "Surface loads", the load samples: an ellipse (semi-axes 0.4, 0.2, 200 vertices) at a pose (turned by
!< 0.3 rad, moved by 0.05 along u) on cell-centred lines of spacing 0.02: three Gauss points per segment with two probes
!< along the normal, their visibility and cells, the first sample, the device path against the host; then a
!< marching-cubes sphere, the samples of order 3 left resident on the device and brought back.
program ref_surface_loads
use fossil, only : surface_stl_object, curve_object, facet_object, extract_isosurface, load_samples_t, load_samples_2d, &
                   load_samples_3d, load_samples_mismatches, load_samples_reset, rigid_pose_2d_t, pose_2d
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P)
type(curve_object)              :: ellipse
type(surface_stl_object)        :: sphere
type(load_samples_t)            :: samples, twin
type(rigid_pose_2d_t)           :: pose
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: h, u(200), v(200), xc(50), yc(30), xn(25), f(33, 33, 33)
integer(I4P)                    :: i, j, k

u = [(0.4_R8P * cos(2 * PI * (i - 1) / 200), i = 1, 200)] ; v = [(0.2_R8P * sin(2 * PI * (i - 1) / 200), i = 1, 200)]
call ellipse%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
h = 0.02_R8P
xc = [(-0.5_R8P + (i - 0.5_R8P) * h, i = 1, 50)] ; yc = [(-0.3_R8P + (j - 0.5_R8P) * h, j = 1, 30)]
pose = pose_2d(0.3_R8P, [0.05_R8P, 0._R8P])                       ! the body at its pose: turned, moved
call load_samples_2d(ellipse, xc, yc, samples, order=3_I4P, &     ! three Gauss points per segment
                     lengths=[h, 2 * h], pose=pose)               ! two probes along the normal, + into the fluid
print '(A,I0,A,I0,A,I0,A,F10.7)', 'ellipse: segments ', samples%ne, ', samples ', samples%n, ', visible probes ', &
      count(samples%visible), ', sum of weights ', sum(samples%weight)
print '(A,2F9.5,A,F9.6,A,2F9.5,A,I0)', 'sample 1: point', samples%points(:, 1), '  weight', samples%weight(1), &
      '  normal', samples%normal(:, 1), '  segment ', samples%element(1)
print '(A,2F9.5,A,2F9.5,A,2(2I3,:,","))', '  probes', samples%probe(:, 1, 1), ',', samples%probe(:, 2, 1), &
      '  cells', samples%probe_cell(:, :, 1)
call ellipse%enter_device
call load_samples_2d(ellipse, xc, yc, twin, order=3_I4P, lengths=[h, 2 * h], pose=pose, device=.true.)
print '(A,I0)', 'device path vs host path: mismatches ', load_samples_mismatches(samples, twin)
call ellipse%exit_device
do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33                         ! a sphere of radius 0.5 by marching cubes
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call sphere%adopt_facets(facets=facets)
xn = [(-0.6_R8P + 0.05_R8P * i, i = 0, 24)]
call sphere%enter_device                                          ! the surface on the device, once
call load_samples_3d(sphere, xn, xn, xn, samples, order=3_I4P, lengths=[0.05_R8P, 0.1_R8P], &
                     resident=.true.)                             ! the samples stay on the device
print '(A,I0,A,I0,A,L1)', 'sphere: facets ', samples%ne, ', samples ', samples%n, ', resident ', samples%on_device
call samples%update_host                                          ! (a device solver would read them there)
call load_samples_3d(sphere, xn, xn, xn, twin, order=3_I4P, lengths=[0.05_R8P, 0.1_R8P])
print '(A,F9.6,A,F9.6,A,I0)', '  area ', sum(samples%weight), ' (4 pi r^2 =', PI, '), mismatches against the host path ', &
      load_samples_mismatches(samples, twin)
call load_samples_reset(samples)                                  ! releases the device copies
call sphere%exit_device
endprogram ref_surface_loads
