!run ref_fd_closure ref_fd_closure
!< Reference page "Finite-difference wall closure": the ghost points of the NACA 0012 section (the wing cut at z = 0.5)
!< on cell-centred lines of spacing chord/64 with reach 3 (the 7-point stencil of a sixth-order scheme), then of a
!< marching-cubes sphere on node lines with reach 2, the device path checked against the host one.
program ref_fd_closure
use fossil, only : surface_stl_object, curve_object, slice_surface, extract_isosurface, facet_object, ghost_points_t, &
                   ghost_points_2d, ghost_points_3d, GHOST_YM, GHOST_YP
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: wing, sphere
type(curve_object)              :: section
type(ghost_points_t)            :: gp, gd
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: h, xc(80), yc(16), xn(25), f(33, 33, 33)
integer(I4P)                    :: i, j, k, g, layer(3)

!region 2d
call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, section)            ! the airfoil in the plane z = 0.5, chord -0.25 .. 0.75
h = 1._R8P / 64._R8P
xc = [(-0.375_R8P + (i - 0.5_R8P) * h, i = 1, 80)]           ! the solver's cell centres
yc = [(-0.125_R8P + (j - 0.5_R8P) * h, j = 1, 16)]
call ghost_points_2d(section, xc, yc, [3_I4P, 3_I4P], gp)    ! reach 3: three ghost layers
print '(A,I0,A,I0,A,I0)', 'cell centres ', size(xc) * size(yc), ', ghost points ', gp%n, ', reach ', gp%reach(1)
print '(A,I0,A,I0)', '  seeing fluid below (-v) ', count(iand(gp%mask, GHOST_YM) /= 0), &
                     ', above (+v) ', count(iand(gp%mask, GHOST_YP) /= 0)
layer = 0
do g = 1, gp%n                                               ! the layer: steps to the nearest fluid point
   k = minval(gp%steps(:, g), mask=gp%steps(:, g) > 0) ; layer(k) = layer(k) + 1
enddo
print '(A,3(I0,:,", "))', '  layers 1, 2, 3: ', layer
g = 1
print '(A,2I4,A,2F9.5,A,B4.4,A,4I2)', 'ghost 1: ij', gp%ijk(:, g), '  point', gp%points(:, g), '  mask ', gp%mask(g), &
                                     '  steps', gp%steps(:, g)
!endregion 2d
!region 3d
do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33                     ! a sphere of radius 0.5 by marching cubes
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call sphere%adopt_facets(facets=facets)
xn = [(-0.6_R8P + 0.05_R8P * i, i = 0, 24)]                  ! node lines, the sphere's poles on them
call ghost_points_3d(sphere, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], gp)
call ghost_points_3d(sphere, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], gd, device=.true.)
print '(A,I0,A,I0)', 'sphere: points ', size(xn)**3, ', ghost points ', gp%n
print '(A,I0,A,I0)', '  first layer ', count([(any(gp%steps(:, g) == 1), g = 1, gp%n)]), &
                     ', second only ', count([(all(gp%steps(:, g) /= 1), g = 1, gp%n)])
print '(A,L1)', '  device path = host path: ', gd%n == gp%n .and. all(gd%ijk == gp%ijk) .and. &
                all(gd%mask == gp%mask) .and. all(gd%steps == gp%steps) .and. all(gd%points == gp%points)
!endregion 3d
endprogram ref_fd_closure
