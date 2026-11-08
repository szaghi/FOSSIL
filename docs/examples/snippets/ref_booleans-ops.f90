call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)) ; call b%analyze ! B overlaps A's corner
do k = 1, 4
   c = a                                                        ! boolean replaces self: work on a copy
   call c%boolean(other=b, op=OPS(k), status=status)
   if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
   print '(A,A,I4,A,F6.4,A,ES8.1,A,L1,A,L1)', NAMES(k), ': facets', c%get_facets_number(), ', volume ', &
         c%get_volume(), ', error ', abs(c%get_volume() - EXACT(k)), ', watertight ', c%is_watertight(), &
         ', manifold ', c%is_manifold()
enddo
