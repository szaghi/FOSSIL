surface = dragon
normal = ex_R8P                                           ! through the yz plane (normal = x axis)
call surface%mirror(normal=normal)
call surface%analyze
