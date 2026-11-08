call mesh%mass_properties(mass=mass, center=centre, principal_moments=moments, density=2700._R8P, status=status)
if (status /= STATUS_OK) error stop 'not a closed body'
print '(A,F9.2,A,3F8.4)', 'mass ', mass, ', centre of mass', centre
print '(A,3F10.3)', 'principal moments', moments
