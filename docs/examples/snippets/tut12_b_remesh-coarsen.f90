n0 = mesh%get_facets_number()
call mesh%isotropic_remesh(target_length=0.4_R8P, iterations=6_I4P, preserve_features=.true., status=status)
print '(A,I0,A,I0,A)', 'target 0.4 on the remeshed part: ', n0, ' -> ', mesh%get_facets_number(), &
                       ' facets (no collapse, issue #79)'
