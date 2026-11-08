call merge_targets(links, cut%piece_cell, cut%piece_area, cut%piece_wet_normal, x, x, mt)
print '(A,I0,A,I0,A,I0)', 'merge: served ', count(mt%status /= MERGE_NONE), ', direct ', count(mt%status == MERGE_OK), &
                          ', through a chain ', count(mt%status == MERGE_CHAIN)
print '(A,2I3,A,I0)', '  smallest piece -> cell', cut%piece_cell(:, mt%final_target(p)), &
                      ', normal-directed targets ', mt%mix_start(p + 1) - mt%mix_start(p)
