do k = 10, 14, 2                                             ! the first k - 1 cell layers of the O-grid only
   call thin(1)%destroy
   call thin(1)%initialize(u(:, 1:k), v(:, 1:k), name='thin')
   call detect_topology_2d(thin(1:1))
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
   call thin(1)%add_record(rec)
   rec%face = 4 ; rec%lo = [1, k, 1] ; rec%hi = [NI, k, 1] ; rec%name = 'overlap'
   call thin(1)%add_record(rec)
   thin(2) = blocks(2)
   call thin(1)%build_search
   call wall_cutter_2d(thin(1:1), wall)
   call hole_flags_2d(thin, wall, holes, cuts=[.false., .true.])
   call donor_candidates_2d(thin, donors)
   call classify_2d(thin, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   print '(A,I0,A,I0,A,I0,A,I0,A,I0,A)', 'O-grid of ', k - 1, ' layers: receptors ', count(masks(1)%iota < 0), &
         ', orphans ', count(masks(1)%orphan > 0), ' (stencil on a hole ', &
         count(masks(1)%orphan == ORPHAN_STENCIL_HOLE), ', receptor in the stencil ', &
         count(masks(1)%orphan == ORPHAN_RECEPTOR_DONOR), ')'
enddo
endblock
