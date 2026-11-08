surface = dragon
matrix = mirror_matrix_R8P(normal=ex_R8P)
call surface%mirror(matrix=matrix)                        ! the same as normal=ex_R8P
call surface%analyze
