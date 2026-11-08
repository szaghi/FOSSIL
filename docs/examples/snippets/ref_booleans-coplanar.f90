call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(0.5_R8P, 0._R8P, 0._R8P)) ; call b%analyze   ! 4 faces coplanar with A's
c = a
call c%boolean(other=b, op=BOOL_DIFFERENCE, status=status)
print '(A,I3,A,F6.4,A,L1,A,L1)', 'coplanar difference: facets', c%get_facets_number(), ', volume ', c%get_volume(), &
      ', watertight ', c%is_watertight(), ', manifold ', c%is_manifold()
