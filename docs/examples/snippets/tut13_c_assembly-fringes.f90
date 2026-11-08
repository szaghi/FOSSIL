call fringe_levels_2d(blocks, masks)                         ! INNER 1, 2 around the holes; OUTER 1, 2 at the overlap
print '(A,4(1X,I0))', 'fringe INNER 1, INNER 2, OUTER 1, OUTER 2:', &
      sum([(count(masks(b)%fringe == FRINGE_INNER_1), b = 1, 2)]), sum([(count(masks(b)%fringe == FRINGE_INNER_2), b = 1, 2)]), &
      sum([(count(masks(b)%fringe == FRINGE_OUTER_1), b = 1, 2)]), sum([(count(masks(b)%fringe == FRINGE_OUTER_2), b = 1, 2)])
call overlap_width_2d(blocks, masks, d)                      ! Chesshire's width between the two fringes
print '(A,F7.4)', 'overlap width ', d(1, 2)
call build_connectivity(masks, donors, 3_I4P, 2_I4P, conn)   ! the receptor list a solver reads
call save_connectivity(conn, 'tut13-connectivity.txt')
open(newunit=un, file='tut13-connectivity.txt', action='read')
do i = 1, 5                                                  ! the header and the first receptor
   read(un, '(A)') line
   print '(A)', trim(line)
enddo
close(un)
