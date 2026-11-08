call srd_neighbourhoods(links, cut%piece_cell, cut%piece_volume, cut%piece_centroid, cut%piece_wet_normal, g, g, srd, z=g)
call merge_targets(links, cut%piece_cell, cut%piece_volume, cut%piece_wet_normal, g, g, mt, z=g)
print '(A,I0,A,I0,A,ES9.2)', 'SRD: small ', count(srd%status == SRD_COMPLETE), ', incomplete ', &
                             count(srd%status == SRD_INCOMPLETE), ', sum vhat - fluid volume ', &
                             sum(srd%vhat) - sum(cut%piece_volume)
print '(A,I0,A,I0)', 'merge: served ', count(mt%status /= MERGE_NONE), ', direct ', count(mt%status == MERGE_OK)
