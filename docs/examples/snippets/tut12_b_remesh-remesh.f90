do it = 1, 6
   if (it == 3 .or. it == 5) cycle
   mesh = part                                                              ! the part stays the reference
   call mesh%isotropic_remesh(target_length=0.2_R8P, iterations=it, preserve_features=.true., status=status)
   if (status /= REM_STATUS_OK) error stop 'remesh failed'
   print '(A,I0,A,I5,A,F7.4,A,F7.4,A,L1)', 'iterations ', it, ': ', mesh%get_facets_number(), ' facets, edges from ', &
         mesh%smallest_edge_len(), ' to', mesh%largest_edge_len(), ', watertight ', mesh%is_watertight()
enddo
