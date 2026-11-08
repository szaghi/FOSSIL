call split_box(a, 'low', 'high', vector_R8P(0._R8P, 0._R8P, 0._R8P))         ! z < 0.5 and z > 0.5 of a unit box
call split_box(b, 'b_low', 'b_high', vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))   ! the same, moved by (0.5, 0.5, 0.5)
call a%boolean(other=b, op=BOOL_UNION, status=status)                       ! each sub-triangle keeps its patch
if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
call areas('union', a)
