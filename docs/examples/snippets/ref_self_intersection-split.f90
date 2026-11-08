call split_shells(body, shell)                              ! one surface per shell
call shell(1)%boolean(other=shell(2), op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
call report('union:   ', shell(1))
