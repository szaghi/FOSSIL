call field%compute(surface=cube, sweeps=sweeps, fallback_cells=fallback, store_closest=.true., status=status)
if (status /= STATUS_OK) error stop 'distance field failed'
print '(A,I0,A,F7.3,A,F7.3)', 'grid ', n, '^3 on ', x(1), ' to ', x(n)
print '(A,I0,A,F6.3,A,I0,A,I0)', 'band cells ', field%get_band_cells(), ' (half-width ', field%get_band_width(), &
                                 '), point-query cells ', fallback, ', sweeps ', sweeps
