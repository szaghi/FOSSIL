call field%compute_band(curve=body, status=status)                ! exact distances within the band
if (status /= STATUS_OK) error stop 'compute_band failed'
signed = field                                                    ! a copy: each stage kept for the frames
call signed%compute_signs(curve=body, fallback_cells=fallback, status=status)  ! flood fill of the signs
final = signed
call final%sweep(sweeps=sweeps, status=status)                    ! fast sweeping to convergence
