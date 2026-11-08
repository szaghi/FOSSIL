bmin%x = -15._R8P ; bmin%y = -5._R8P ; bmin%z =  0._R8P
bmax%x =   0._R8P ; bmax%y =  5._R8P ; bmax%z = 20._R8P
call head%load_from_file(file_name='dragon.stl', guess_format=.true., clip_min=bmin, clip_max=bmax)
print '(A,I0)', 'clipped load           facets ', head%get_facets_number()
bmin = head%get_bmin() ; bmax = head%get_bmax()
print '(A,F9.4,A,F9.4,A)', 'clipped load           x in [', bmin%x, ',', bmax%x, ']'
