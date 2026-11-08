call merge_targets(links, cut%piece_cell, cut%piece_area, cut%piece_wet_normal, h, h, mt)
if (mt%status(p) == MERGE_OK) print '(A,2I3,A,I0,A,F5.2)', '  merge target in cell', cut%piece_cell(:, mt%target(p)), &
   ', normal-directed targets ', mt%mix_start(p + 1) - mt%mix_start(p), ', weight ', mt%mix_weight(mt%mix_start(p))
