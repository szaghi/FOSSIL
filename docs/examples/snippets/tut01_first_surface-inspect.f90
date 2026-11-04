bmin = surface%get_bmin() ; bmax = surface%get_bmax()
print '(A,I0)',       'facets      ', surface%get_facets_number()
print '(A,3F10.4)',   'box min     ', bmin%x, bmin%y, bmin%z
print '(A,3F10.4)',   'box max     ', bmax%x, bmax%y, bmax%z
print '(A,F12.4)',    'volume      ', surface%get_volume()
print '(A,L1,A,L1)',  'watertight  ', surface%is_watertight(), '   manifold ', surface%is_manifold()
