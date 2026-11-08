call frd_neighbourhoods(links, cut%piece_cell, cut%piece_area, h, h, frd)
print '(A,I0,A,F6.3)', '  FRD set size ', frd%start(p + 1) - frd%start(p), ', own weight ', frd%weight(frd%start(p))
