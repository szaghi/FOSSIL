open(newunit=u, file='maybe_bad.stl', action='write')  ! one facet with a NaN coordinate
write(u, '(A)') 'solid bad', ' facet normal 0 0 1', '  outer loop', '   vertex 0 0 0', '   vertex 1 0 0', &
                '   vertex NaN 1 0', '  endloop', ' endfacet', 'endsolid bad'
close(u)
call surface%load_from_file(file_name='maybe_bad.stl', guess_format=.true., status=stat)
if (stat == STATUS_INVALID_INPUT) then
   write(*, '(A)') 'STL contains non-finite coordinates; aborting.'
   stop
endif
