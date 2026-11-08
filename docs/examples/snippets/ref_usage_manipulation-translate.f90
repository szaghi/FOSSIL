surface = dragon
delta%x = 1._R8P ; delta%y = 0._R8P ; delta%z = 0._R8P
call surface%translate(delta=delta)                       ! by a 3-D vector
call surface%translate(x=2._R8P)                          ! by scalar components, any combination
call surface%translate(y=1.5_R8P)
call surface%translate(x=1._R8P, y=2._R8P, z=0.5_R8P)
call surface%analyze                                      ! box, volume, tree of the moved surface
