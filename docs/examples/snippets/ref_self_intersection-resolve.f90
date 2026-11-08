call body%resolve_self_intersections(status=status)         ! the boundary of the region wound at least once
if (status /= BOOL_STATUS_OK) error stop 'resolve_self_intersections failed'
call report('resolved: ', body)
