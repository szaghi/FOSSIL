do k = 1, 3
   call cube(k)%load_from_file(file_name='cube-01.stl', guess_format=.true.)
enddo
call cube(2)%translate(delta=vector_R8P(1.25_R8P, 0._R8P, 0._R8P)) ; call cube(2)%analyze   ! 0.25 from cube 1
call cube(3)%translate(delta=vector_R8P(0._R8P, 1.6_R8P, 0._R8P))  ; call cube(3)%analyze   ! 0.6 from cube 1
all = cube(1)
call all%merge_solids(other=cube(2)) ; call all%merge_solids(other=cube(3))   ! one surface, three bodies
call all%proximity_pairs(cutoff=0.62_R8P, pairs=pairs, distances=distances, status=status)
if (status /= STATUS_OK) error stop 'proximity_pairs failed'
print '(A,I0,A)', 'one surface, ', all%get_bodies_number(), ' bodies, cutoff 0.62:'
do i = 1, size(distances)
   print '(A,I0,A,I0,A,F7.4)', '  bodies ', pairs(1, i), ' and ', pairs(2, i), ': ', distances(i)
enddo
call proximity_pairs_of(surfaces=cube, cutoff=0.62_R8P, pairs=pairs, distances=distances)
print '(A)', 'three surfaces, cutoff 0.62:'
do i = 1, size(distances)
   print '(A,I0,A,I0,A,F7.4)', '  surfaces ', pairs(1, i), ' and ', pairs(2, i), ': ', distances(i)
enddo
