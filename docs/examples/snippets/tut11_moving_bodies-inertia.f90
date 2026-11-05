call surface%mass_properties(mass=mass, status=status)
print '(A,L1)', 'the placed dragon is a closed body: ', status == STATUS_OK               ! chapter 2: it is not
call surface%alpha_wrap(alpha=0.03_R8P, offset=0.003_R8P, wrapped=wrap, status=status)     ! a closed surrogate
call wrap%mass_properties(mass=mass, center=centre, inertia=inertia, principal_moments=moments, &
                          density=1000._R8P, status=status)                      ! a body of density 1000
if (status /= STATUS_OK) error stop 'the wrap is not closed'
print '(A,ES11.4,A,3F8.4)', 'mass ', mass, ', centre of mass', centre
print '(A,3ES11.4)', 'principal moments', moments
call wrap%move(p1)
call wrap%mass_properties(inertia=inertia1, density=1000._R8P)
print '(A,ES9.2)', 'moved wrap against inertia_to_world: ', maxval(abs(inertia1 - inertia_to_world(inertia, p1))) / &
                                                           maxval(abs(inertia))
