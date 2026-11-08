call a%isotropic_remesh(target_length=0.15_R8P, iterations=3_I4P, preserve_features=.true., status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
call areas('remeshed', a)                                                     ! patch boundaries stay where they are
