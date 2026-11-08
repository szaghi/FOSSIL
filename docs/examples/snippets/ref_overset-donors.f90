call overlap_matrix_2d(blocks, pair)
print '(A,4I2)', 'pair matrix ', pair
call donor_candidates_2d(blocks, donors, pair=pair, wall=wall)
print '(A,I0,A,I0)', 'candidates: O-grid vertices in the background ', size(donors(1)%block), &
                     ', background vertices in the O-grid ', size(donors(2)%block)
print '(A,2F8.4)', 'donor / receptor volume, O-grid vertices: min max ', minval(donors(1)%ratio), maxval(donors(1)%ratio)
