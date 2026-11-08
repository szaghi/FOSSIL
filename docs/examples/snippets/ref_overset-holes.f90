call wall_cutter_2d(blocks(1:1), wall)
call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.], distances=dist)
call hole_flags_2d(blocks, wall, holes_offset, offset=0.3_R8P, cuts=[.false., .true.])
print '(A,I0,A,I0,A,F7.4)', 'holes ', count(holes(2)%v), ', with offset 0.3 ', count(holes_offset(2)%v), &
                            ', nearest field point to the wall ', minval(dist(2)%v, mask=.not. holes(2)%v)
