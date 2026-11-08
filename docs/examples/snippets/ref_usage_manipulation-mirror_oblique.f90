surface = dragon
normal = ex_R8P + ey_R8P                                  ! through an oblique plane (normalised inside)
call surface%mirror(normal=normal)
call surface%analyze
