!run ck_mass_properties_2d ck_mass_properties_2d
program ck_mass_properties_2d
use fossil, only : curve_object, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(curve_object) :: tube
real(R8P)          :: area, center(2), sm(3), polar
integer(I4P)       :: status

! A hollow rectangular section 0.1 x 0.2 m, walls 0.01 m: the outer loop counter-clockwise, the hole clockwise.
call tube%adopt_loops(u=[-0.05_R8P,  0.05_R8P, 0.05_R8P, -0.05_R8P,   -0.04_R8P, -0.04_R8P, 0.04_R8P,  0.04_R8P], &
                      v=[-0.10_R8P, -0.10_R8P, 0.10_R8P,  0.10_R8P,   -0.09_R8P,  0.09_R8P, 0.09_R8P, -0.09_R8P], &
                      loop_sizes=[4_I4P, 4_I4P])
call tube%mass_properties(area=area, center=center, second_moments=sm, polar_moment=polar, status=status)
if (status /= STATUS_OK) error stop 'not a closed section'
print '(A,ES12.5,A)',  'area                 ', area, ' m^2'
print '(A,2F9.5)',     'centroid             ', center
print '(A,ES12.5,A)',  'int v^2 dA (bending about u) ', sm(2), ' m^4'
print '(A,ES12.5,A)',  'int u^2 dA (bending about v) ', sm(1), ' m^4'
print '(A,ES12.5,A)',  'polar moment         ', polar, ' m^4'
print '(A,ES12.5)',    'by hand: (b h^3 - b''h''^3)/12 = ', (0.1_R8P * 0.2_R8P ** 3 - 0.08_R8P * 0.18_R8P ** 3) / 12
endprogram ck_mass_properties_2d
