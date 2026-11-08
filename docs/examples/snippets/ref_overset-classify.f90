call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes_offset, iterations=iterations, &
                 circular=circular)
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'explicit: receptors ', count(masks(1)%iota < 0), ' + ', count(masks(2)%iota < 0), &
                                    ', orphans ', count(masks(1)%orphan > 0) + count(masks(2)%orphan > 0), &
                                    ', iterations ', iterations, ', circular ', circular
call classify_2d(blocks, donors, implicit, SELECT_IMPLICIT, 3_I4P, iterations=iterations)
print '(A,I0,A,I0,A,I0,A,I0)', 'implicit: receptors ', count(implicit(1)%iota < 0), ' + ', &
                               count(implicit(2)%iota < 0), ', holes ', count(implicit(2)%iota == 0), &
                               ', iterations ', iterations
call fringe_levels_2d(blocks, masks, levels=2_I4P)
print '(A,4I5)', 'fringe INNER 1, INNER 2 (background), OUTER 1, OUTER 2 (O-grid) ', &
                 count(masks(2)%fringe == FRINGE_INNER_1), count(masks(2)%fringe == FRINGE_INNER_2), &
                 count(masks(1)%fringe == FRINGE_OUTER_1), count(masks(1)%fringe == FRINGE_OUTER_2)
call overlap_width_2d(blocks, masks, width)
print '(A,F8.4)', 'overlap width ', width(1, 2)
