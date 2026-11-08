call field%compute(surface=surface, sign_algorithm=SIGN_FAST_WINDING_NUMBER, sweeps=sweeps, fallback_cells=fallback, &
                   status=status)
if (status /= STATUS_OK) error stop 'distance field failed'
call field%get_dims(ni=ni, nj=nj, nk=nk)
print '(A,I0,A,I0,A,I0)', 'grid         ', ni, ' x ', nj, ' x ', nk
print '(A,F8.4,A,I0,A)',  'band         ', field%get_band_width(), ' half-width, ', field%get_band_cells(), ' cells'
print '(A,I0,A,I0)',      'sweeps       ', sweeps, ', point-query signs ', fallback
print '(A,I0)',           'inside cells ', count(field%phi < 0._R8P)
