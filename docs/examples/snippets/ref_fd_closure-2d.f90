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
