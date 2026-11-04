print '(A,I0,A,I0)', 'facets      ', facets_before, ' -> ', surface%get_facets_number()
print '(A,I0)',      'degenerate  ', surface%get_degenerate_facets_removed()
print '(A,I0)',      'duplicate   ', surface%get_duplicate_facets_removed()
print '(A,I0)',      'non-manif.  ', surface%get_non_manifold_edges_number()
print '(A,F12.4)',   'volume      ', surface%get_volume()
print '(A,L1,A,L1,A,L1)', 'watertight  ', surface%is_watertight(), '   manifold ', surface%is_manifold(), &
                          '   volume ', surface%is_volume()
call surface%save_into_file(file_name='dragon-clean.stl')
print '(A)', 'saved dragon-clean.stl'
