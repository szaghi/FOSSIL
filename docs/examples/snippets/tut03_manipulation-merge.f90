twin = half
call twin%mirror(normal=ey_R8P)                                     ! reflect through the plane y = 0
call twin%analyze
call half%merge_solids(other=twin)                                  ! merge_solids analyzes
print '(A,I0,A)', 'merged    ', half%get_facets_number(), ' facets (half + mirrored half)'
