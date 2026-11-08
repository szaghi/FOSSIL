call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
a_boss = 16 * 0.35_R8P**2 * sin(2 * PI / 32)                       ! area of the 32-gon of radius 0.35
v_exact = 1._R8P + a_boss * 0.4_R8P                                 ! the boss stands 0.4 above the base
call report('base + boss', part)
print '(A,ES9.2)', '  volume - exact polygonal volume: ', part%get_volume() - v_exact
