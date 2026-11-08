call wall_cutter_2d(blocks(1:1), cutter)                     ! the cutter: the edges of the faces named 'wall'
print '(A,I0,A)', 'cutter: ', cutter%get_vertices_number(), ' vertices, one closed loop'
call hole_flags_2d(blocks, cutter, holes, cuts=[.false., .true.]) ! the background only: the O-grid's wall is the cutter
print '(A,I0,A,I0)', 'background holes ', count(holes(2)%v), ' of ', size(holes(2)%v)
call hole_map_2d(blocks, cutter, holes_map, [64, 64], cuts=[.false., .true.])  ! a Cartesian hole map instead
print '(A,L1)', 'hole map gives the same holes: ', all(holes_map(2)%v .eqv. holes(2)%v)
