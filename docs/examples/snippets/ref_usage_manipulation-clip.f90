surface = dragon
bmin%x = -15._R8P ; bmin%y = -5._R8P ; bmin%z =  0._R8P
bmax%x =   0._R8P ; bmax%y =  5._R8P ; bmax%z = 20._R8P
call surface%clip(bmin=bmin, bmax=bmax, remainder=remainder)  ! clip re-runs analyze on both
print '(A,I0,A,I0,A,I0)', 'clip           kept ', surface%get_facets_number(), ', remainder ', &
                          remainder%get_facets_number(), ', of ', dragon%get_facets_number()
