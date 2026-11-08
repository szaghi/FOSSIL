mesh = part
call mesh%isotropic_remesh(target_length=0.2_R8P, iterations=6_I4P, preserve_features=.true., status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
print '(A,I0,A,F7.4,A,L1)', 'remeshed: ', mesh%get_facets_number(), ' facets, longest edge ', &
                            mesh%largest_edge_len(), ', watertight ', mesh%is_watertight()
