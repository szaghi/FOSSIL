call part%boolean(bore, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
a_bore = 16 * 0.31_R8P**2 * sin(2 * PI / 32)
v_exact = v_exact - a_bore * 0.5_R8P                                ! the bore goes through the 0.5 of the base
call report('base + boss - bore', part)
print '(A,ES9.2)', '  volume - exact polygonal volume: ', part%get_volume() - v_exact
