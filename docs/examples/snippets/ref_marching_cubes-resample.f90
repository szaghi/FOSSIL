call body%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! the unit cube, 48 facets
do m = 0, 2
   n = 12 * 2**m                                                ! grid points along the longest side
   call body%resample_via_distance_field(resolution=n, surface_out=resampled, status=status)
   if (status /= MC_STATUS_OK) error stop 'resample failed'
   print '(A,I2,A,I5,A,L1,A,F7.5,A,ES8.2)', 'cube at ', n, ': ', resampled%get_facets_number(), ' facets, watertight ', &
         resampled%is_watertight(), ', volume ', resampled%get_volume(), ', deficit ', 1._R8P - resampled%get_volume()
enddo
