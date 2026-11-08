call mesh%mass_properties(mass=mass, center=centre, principal_moments=moments, principal_axes=axes, &
                          density=2700._R8P, status=status)                 ! aluminium, SI units: lengths in m
if (status /= STATUS_OK) error stop 'not a closed body'
print '(A,F9.2,A,3F8.4)', 'mass ', mass, ' kg, centre of mass', centre
do k = 1, 3
   print '(A,F9.3,A,3F8.4)', 'principal moment ', moments(k), ' about', merge(0._R8P, axes(:, k), abs(axes(:, k)) < 5.e-5_R8P)
enddo
call part%mass_properties(mass=mass_part, density=2700._R8P)
print '(A,ES9.2)', 'mass of the part before the remesh, relative difference: ', (mass - mass_part) / mass_part
