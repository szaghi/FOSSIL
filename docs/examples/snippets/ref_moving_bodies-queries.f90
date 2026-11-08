call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
i = 0
do k = 0, 2 ; do j = 0, 2 ; do l = 0, 2                                     ! 27 world points across the moved cube
   i = i + 1
   points(i) = vector_R8P(-0.1_R8P + 0.6_R8P * l, -0.1_R8P + 0.6_R8P * j, -0.1_R8P + 0.6_R8P * k)
enddo ; enddo ; enddo
call cube%distance_many(points=points, distances=d_posed, is_signed=.true., is_square_root=.true., pose=p1)
moved = cube
call moved%move(p1)                                                              ! refit, not rebuilt
call moved%distance_many(points=points, distances=d_moved, is_signed=.true., is_square_root=.true.)
rebuilt = moved
call rebuilt%analyze                                                             ! the tree rebuilt from scratch
call rebuilt%distance_many(points=points, distances=d_rebuilt, is_signed=.true., is_square_root=.true.)
print '(A,ES9.2,A,I0)', 'posed - moved: ', maxval(abs(d_posed - d_moved)), ', inside ', count(d_posed < 0._R8P)
print '(A,L1)', 'moved = rebuilt, bit for bit: ', all(d_moved == d_rebuilt)
