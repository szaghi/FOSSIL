surface = dragon
call surface%resize(x=0.5_R8P, z=1.2_R8P)                 ! per axis; y keeps its size
call surface%analyze
