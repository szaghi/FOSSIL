call holed%alpha_wrap(alpha=-1._R8P, wrapped=wrapped, status=status)
print '(A,L1)', 'alpha < 0:   AWRAP_STATUS_BAD_INPUT  ', status == AWRAP_STATUS_BAD_INPUT
call holed%alpha_wrap(alpha=2._R8P, wrapped=wrapped, status=status)
print '(A,L1)', 'alpha > box: AWRAP_STATUS_DEGENERATE ', status == AWRAP_STATUS_DEGENERATE
