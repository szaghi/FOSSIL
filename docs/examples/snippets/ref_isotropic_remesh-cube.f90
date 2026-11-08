do k = 0, 1
   features = k == 1
   call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call cube%isotropic_remesh(target_length=0.2_R8P, iterations=5_I4P, preserve_features=features, status=status)
   if (status /= REM_STATUS_OK) error stop 'remesh failed'
   print '(A,L1,A,I4,A,F6.4,A,F6.4,A,F7.5,A,F7.5)', 'cube, preserve_features ', features, ': ', &
         cube%get_facets_number(), ' facets, edges ', cube%smallest_edge_len(), ' to ', cube%largest_edge_len(), &
         ', area ', cube%get_area(), ', volume ', cube%get_volume()
enddo
