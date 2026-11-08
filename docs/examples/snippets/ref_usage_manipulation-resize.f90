surface = dragon
factor = 2._R8P * ex_R8P + 2._R8P * ey_R8P + 2._R8P * ez_R8P
call surface%resize(factor=factor)                        ! uniform, about the origin
call surface%analyze
