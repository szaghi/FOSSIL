call cut_links_2d(cut, links)
call cut_coarsen_2d(cut, links, x, x, 2_I4P, coarse, coarse_links, fine_to_coarse, status)
call section%cut_cells(x=x(1:49:2), y=x(1:49:2), cut=direct)                      ! the coarse grid directly
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'coarse (ratio 2): status ', status, ', pieces ', coarse%np, ', direct ', direct%np, &
                               ', fine piece ', p, ' -> coarse piece ', fine_to_coarse(p)
print '(A,ES9.2)', '  largest |alpha coarsened - alpha direct| ', maxval(abs(coarse%alpha - direct%alpha))
