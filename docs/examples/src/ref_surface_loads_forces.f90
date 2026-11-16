!run ref_surface_loads_forces ref_surface_loads_forces
!< Reference page "Surface loads", forces and moments (`integrate_loads`): a marching-cubes sphere of radius 0.5 at a pose
!< (turned, moved to (0.3, 0.2, 0.1)) in the pressure p = 1 - z of a fluid at rest under gravity -z, samples of order 2:
!< the force against the buoyancy (0, 0, V), the moment about the origin against V x_c x (0, 0, 1), the wetted area; then
!< the same sums with the samples resident on the device, compared by bits.
program ref_surface_loads_forces
use fossil, only : surface_stl_object, facet_object, extract_isosurface, load_samples_t, load_samples_3d, &
                   load_samples_reset, integrate_loads, rigid_pose_t, pose_from_axis_angle, to_world
use penf,   only : I4P, I8P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: sphere
type(load_samples_t)            :: ls, lr
type(rigid_pose_t)              :: pose
type(facet_object), allocatable :: facets(:)
real(R8P),          allocatable :: p(:)
real(R8P)                       :: f(33, 33, 33), force(3), moment(3), fd(3), md(3), area, vol, xc(3)
integer(I4P)                    :: i, j, k

do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33                         ! a sphere of radius 0.5 by marching cubes
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call sphere%adopt_facets(facets=facets)
pose = pose_from_axis_angle([1._R8P, 1._R8P, 0._R8P], 0.5_R8P, [0.3_R8P, 0.2_R8P, 0.1_R8P])
!region loads
call load_samples_3d(sphere, samples=ls, order=2_I4P, pose=pose)  ! the samples of the body at its pose
p = 1._R8P - ls%points(3, :)                                      ! the solver's pressure at the samples
call integrate_loads(ls, force, moment, pressure=p, about=[0._R8P, 0._R8P, 0._R8P], area=area)
call sphere%mass_properties(volume=vol, center=xc) ; xc = to_world(pose, xc)
print '(A,3F11.7,A,F10.7,A)', 'force ', clean(force), '   buoyancy (0, 0, V), V =', vol, ' (faceted)'
print '(A,3F11.7,A,3F11.7)', 'moment', clean(moment), '   V x_c x (0, 0, 1) =', clean(vol * [xc(2), -xc(1), 0._R8P])
print '(A,F10.7,A,L1)', 'wetted area', area, ';  the identities to 1e-15: ', &
      max(maxval(abs(force - [0._R8P, 0._R8P, vol])), maxval(abs(moment - vol * [xc(2), -xc(1), 0._R8P]))) < 1.e-15_R8P
call sphere%enter_device
call load_samples_3d(sphere, samples=lr, order=2_I4P, pose=pose, resident=.true.)
call integrate_loads(lr, fd, md, pressure=p, about=[0._R8P, 0._R8P, 0._R8P], device=.true.)  ! read on the device
print '(A,L1)', 'device = host, bit for bit: ', all(transfer([fd, md], 0_I8P, 6) == transfer([force, moment], 0_I8P, 6))
call load_samples_reset(lr)
call sphere%exit_device
!endregion loads

contains
   pure function clean(x) result(y)
   !< Rounded to 7 decimals, no negative zero (stable output).
   real(R8P), intent(in) :: x(:)       !< Values.
   real(R8P)             :: y(size(x)) !< Rounded.

   y = anint(x * 1.e7_R8P) / 1.e7_R8P + 0._R8P
   endfunction clean
endprogram ref_surface_loads_forces
