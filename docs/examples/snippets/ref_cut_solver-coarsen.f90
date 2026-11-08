call cut_coarsen_2d(cut, links, h, h, 2_I4P, coarse, coarse_links, fine_to_coarse, status)
print '(A,I0,A,I0,A,I0,A,F5.2)', 'coarse (ratio 2): status ', status, ', pieces ', coarse%np, ', links ', &
                                 size(coarse_links%piece), ', fraction of the cut column ', coarse%alpha(3, 3)
print '(A,I0,A,I0)', '  the fine piece lands in coarse piece ', fine_to_coarse(p), ', in coarse cell ', &
                     coarse%piece_cell(1, fine_to_coarse(p))
