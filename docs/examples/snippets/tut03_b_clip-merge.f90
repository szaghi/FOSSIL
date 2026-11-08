twin = half
call twin%mirror(normal=ey_R8P)                                     ! reflect through the plane y = 0
call twin%analyze
call half%merge_solids(other=twin)                                  ! merge_solids analyzes
print '(A,I0,A)', 'merged    ', half%get_facets_number(), ' facets (half + mirrored half)'
bmin = half%get_bmin() ; bmax = half%get_bmax()
print '(A,3F9.4)', 'box min   ', bmin%x, bmin%y, bmin%z
print '(A,3F9.4)', 'box max   ', bmax%x, bmax%y, bmax%z
