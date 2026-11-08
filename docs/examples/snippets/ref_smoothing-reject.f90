call taubin%smooth(method=SMOOTH_METHOD_TAUBIN, lambda=0.5_R8P, mu=-0.4_R8P, status=status)   ! mu > -lambda
print '(A,L1)', 'mu = -0.4 > -lambda: SMOOTH_STATUS_BAD_INPUT ', status == SMOOTH_STATUS_BAD_INPUT
