call donor_candidates_2d(blocks, donors)                     ! every vertex located in the other block
print '(A,I0,A,I0)', 'donor candidates: O-grid ', size(donors(1)%block), ', background ', size(donors(2)%block)
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes, iterations=iterations, &
                 circular=circular)                          ! 3 x 3 donor stencils, explicit interpolation
do b = 1, 2
   print '(2X,A,A,I0,A,I0,A,I0,A,I0)', blocks(b)%name(1:10), ': field ', count(masks(b)%iota == 1), ', holes ', &
         count(masks(b)%iota == 0), ', receptors ', count(masks(b)%iota < 0), ', orphans ', count(masks(b)%orphan > 0)
enddo
print '(A,I0,A,I0)', 'fixed point in ', iterations, ' iterations, circular pairs ', circular
