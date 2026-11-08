call detect_topology_2d(ogrid)                            ! the seam: one TOPO_ABUT record per side
do r = 1, ogrid(1)%nrec
   print '(A,I0,A,L1,A,I0,A,2I3)', 'record ', r, ': abut ', ogrid(1)%rec(r)%kind == TOPO_ABUT, ', face ', &
         ogrid(1)%rec(r)%face, ', transform', ogrid(1)%rec(r)%transform(1:2)
enddo
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
call ogrid(1)%add_record(rec)                                ! face 3 (-j): the wall
rec%face = 4 ; rec%lo = [1, NJ, 1] ; rec%hi = [NI, NJ, 1] ; rec%name = 'overlap'
call ogrid(1)%add_record(rec)                                ! face 4 (+j): the boundary that receives
print '(A,I0,A)', 'records after add_record: ', ogrid(1)%nrec, ' (2 abut, wall, overlap)'
