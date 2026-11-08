surface = dragon
matrix = rotation_matrix_R8P(axis=ex_R8P, angle=pi / 2._R8P)
call surface%rotate(matrix=matrix)                        ! the same quarter turn
call surface%analyze
