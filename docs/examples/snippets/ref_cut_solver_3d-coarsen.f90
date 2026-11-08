call cut_coarsen_3d(cut, links, g, g, g, 2_I4P, coarse, coarse_links, f2c, status)
gc = g(1:9:2)
call cut_cells_3d(cube, gc, gc, gc, direct)
print '(A,I0,A,I0,A,I0)', 'coarse: status ', status, ', pieces ', coarse%np, ', direct ', direct%np
print '(A,ES9.2,A,L1)', '  largest |alpha coarsened - alpha direct| ', maxval(abs(coarse%alpha - direct%alpha)), &
                        ', same piece counts ', all(coarse%pieces == direct%pieces)
