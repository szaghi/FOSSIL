surface = dragon
axis  = ex_R8P                                            ! about the x axis
angle = pi / 2._R8P                                       ! a quarter turn, in radians
call surface%rotate(axis=axis, angle=angle)               ! about the origin
call surface%analyze
