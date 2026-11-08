call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(delta=vector_R8P(1._R8P + 1.e-13_R8P, 0._R8P, 0._R8P)) ; call b%analyze ! a gap of 1e-13
c = a
call c%boolean(other=b, op=BOOL_UNION, status=status)
print '(A,I0,A,F6.4)', 'side by side, no weld: shells ', c%get_shells_number(), ', volume ', c%get_volume()
c = a
call c%boolean(other=b, op=BOOL_UNION, status=status, weld=1.e-9_R8P)
print '(A,I0,A,F6.4,A,L1)', 'weld 1e-9:             shells ', c%get_shells_number(), ', volume ', c%get_volume(), &
      ', watertight ', c%is_watertight()
