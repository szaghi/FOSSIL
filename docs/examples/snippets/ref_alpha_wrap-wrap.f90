call holed%alpha_wrap(alpha=0.6_R8P, offset=0.02_R8P, wrapped=wrapped, status=status)   ! ball wider than the hole
call report('alpha 0.60')
call holed%alpha_wrap(alpha=0.15_R8P, offset=0.033_R8P, wrapped=wrapped, status=status) ! ball through the hole
call report('alpha 0.15')
