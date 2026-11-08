call srd_neighbourhoods(links, cut%piece_cell, cut%piece_area, cut%piece_centroid, cut%piece_wet_normal, x, x, srd)
print '(A,I0,A,I0,A,I0)', 'SRD: regular ', count(srd%status == SRD_REGULAR), ', small ', &
                          count(srd%status == SRD_COMPLETE), ', incomplete ', count(srd%status == SRD_INCOMPLETE)
p = minloc(cut%piece_area, 1, mask=cut%piece_area > 0._R8P)          ! the smallest piece with area
print '(A,I0,A,2I3,A,ES9.2,A,I0,A)', '  smallest piece ', p, ' in cell', cut%piece_cell(:, p), ', area / h^2 ', &
                                     cut%piece_area(p) / h**2, ', neighbourhood of ', srd%start(p + 1) - srd%start(p), &
                                     ' pieces'
print '(A,F7.4,A,F7.4)', '  neighbourhood area / h^2 ', &
                         sum(cut%piece_area(srd%member(srd%start(p):srd%start(p + 1) - 1))) / h**2, ' >= tau = ', srd%tau
