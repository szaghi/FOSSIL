explicit = noisy
call explicit%smooth(method=SMOOTH_METHOD_EXPLICIT, lambda=0.5_R8P, iterations=10_I4P, status=status)
if (status /= SMOOTH_STATUS_OK) error stop 'explicit smoothing failed'
call report('explicit, 10    ', explicit)
taubin = noisy
call taubin%smooth(status=status)                              ! Taubin, 5 pairs of lambda 0.5, mu -0.53
if (status /= SMOOTH_STATUS_OK) error stop 'Taubin smoothing failed'
call report('Taubin, 5 pairs ', taubin)
