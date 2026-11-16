!run ref_surface_loads_cut ref_surface_loads_cut
!< Reference page "Surface loads", loads from the cut cells (`cut_loads_3d`, `cut_walls_gyration`): a marching-cubes
!< sphere of radius 0.5 moved to (0.3, 0.2, 0.1), cut by 24^3 cells with its wall terms, in the pressure p = 1 - z of a
!< fluid at rest under gravity -z: the force against the buoyancy (0, 0, V); the moment about the origin with the
!< centroid values only and with the gradient (exact); the same p on load samples of order 2 (`integrate_loads`); then
!< one value per fluid piece, at its centroid, as a cell-centred solver has it.
program ref_surface_loads_cut
use fossil, only : surface_stl_object, facet_object, extract_isosurface, cut_cells_3d_t, cut_cells_3d, cut_walls_3d_t, &
                   cut_walls_gyration, cut_loads_3d, load_samples_t, load_samples_3d, integrate_loads, rigid_pose_t, &
                   pose_from_axis_angle
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: sphere
type(cut_cells_3d_t)            :: cut
type(cut_walls_3d_t)            :: walls
type(load_samples_t)            :: ls
type(facet_object), allocatable :: facets(:)
real(R8P),          allocatable :: p(:), g(:,:), gyration(:,:), pp(:)
real(R8P)                       :: f(33, 33, 33), x(25), force(3), m0(3), m1(3), ffd(3), mfd(3), area, vol, xc(3), mb(3)
integer(I4P)                    :: i, j, k

do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33                         ! a sphere of radius 0.5 by marching cubes
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call sphere%adopt_facets(facets=facets)
call sphere%move(pose_from_axis_angle([1._R8P, 1._R8P, 0._R8P], 0.5_R8P, [0.3_R8P, 0.2_R8P, 0.1_R8P]))
call sphere%mass_properties(volume=vol, center=xc) ; mb = vol * [xc(2), -xc(1), 0._R8P]   ! V x_c x (0, 0, 1)
x = [(-0.53_R8P + 0.0625_R8P * i, i = 0, 24)]
!region cutloads
call cut_cells_3d(sphere, x + 0.1_R8P, x, x - 0.1_R8P, cut, walls=walls)   ! the cut and its wall terms
p = 1._R8P - walls%centroid(3, :)                                 ! the pressure at the walls' centroids
call cut_loads_3d(walls, force, m0, pressure=p, about=[0._R8P, 0._R8P, 0._R8P], area=area)
call cut_walls_gyration(walls, sphere, x + 0.1_R8P, x, x - 0.1_R8P, gyration)   ! for the exact moment
g = spread([0._R8P, 0._R8P, -1._R8P], 2, walls%n)                 ! the pressure gradient on every wall
call cut_loads_3d(walls, force, m1, pressure=p, gradient=g, gyration=gyration, about=[0._R8P, 0._R8P, 0._R8P])
!endregion cutloads
print '(A,I0,A,I0,A,F10.7,A,F10.7)', 'wall terms ', walls%n, ', fluid pieces ', cut%np, '; wetted area', area, ', V', vol
print '(A,3F11.7,A,ES8.1)', 'force ', clean(force), '   error to (0, 0, V)', maxval(abs(force - [0._R8P, 0._R8P, vol]))
print '(A,3F11.7,A,ES8.1)', 'moment, centroid values', clean(m0), '   error', maxval(abs(m0 - mb))
print '(A,3F11.7,A,L1)', 'moment, with the gradient', clean(m1), '   = V x_c x (0, 0, 1) to 1e-15: ', &
      maxval(abs(m1 - mb)) < 1.e-15_R8P
call load_samples_3d(sphere, samples=ls, order=2_I4P)             ! the same p on the facets' samples
call integrate_loads(ls, ffd, mfd, pressure=1._R8P - ls%points(3, :), about=[0._R8P, 0._R8P, 0._R8P])
print '(A,L1)', 'load samples (FD) = cut walls (FV) to 1e-15: ', &
      maxval(abs([ffd - force, mfd - m1])) < 1.e-15_R8P
pp = 1._R8P - cut%piece_centroid(3, :)                            ! one value per fluid piece, at its centroid
call cut_loads_3d(walls, force, m0, pressure_per_piece=pp)
print '(A,3F11.7,A,ES8.1)', 'per piece: force', clean(force), '   error', maxval(abs(force - [0._R8P, 0._R8P, vol]))

contains
   pure function clean(q) result(r)
   !< Rounded to 7 decimals, no negative zero (stable output).
   real(R8P), intent(in) :: q(:)       !< Values.
   real(R8P)             :: r(size(q)) !< Rounded.

   r = anint(q * 1.e7_R8P) / 1.e7_R8P + 0._R8P
   endfunction clean
endprogram ref_surface_loads_cut
