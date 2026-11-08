call holed%fill_holes(holes_filled=filled, holes_skipped=skipped, open_edges=open_edges)
print '(A,I0,A,I0,A,I0)', 'default max_hole_size 0.1: filled ', filled, ', skipped ', skipped, ', open edges ', open_edges
call holed%fill_holes(max_hole_size=0.5_R8P, holes_filled=filled, holes_skipped=skipped)
print '(A,I0,A,I0)', 'max_hole_size 0.5:         filled ', filled, ', skipped ', skipped
call holed%fill_holes(max_hole_size=1._R8P, holes_filled=filled, holes_skipped=skipped)
print '(A,I0,A,I0)', 'max_hole_size 1:           filled ', filled, ', skipped ', skipped
print '(A,I0,A,L1,A,F8.5)', 'result: ', holed%get_facets_number(), ' facets, watertight ', holed%is_watertight(), &
                            ', volume ', holed%get_volume()
