surface = dragon
call surface%rotate(axis=ex_R8P, angle=pi / 2._R8P, center=surface%get_centroid())  ! about the centroid
call surface%analyze
