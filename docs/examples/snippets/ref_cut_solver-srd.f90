h = [(0.1_R8P * i, i = 0, 10)] ! 10 x 10 cells of 0.1; the solid is x < 0.58: fluid fraction 0.2 in the cut column
call half%adopt_loops(u=[-1._R8P, 0.58_R8P, 0.58_R8P, -1._R8P], v=[-1._R8P, -1._R8P, 2._R8P, 2._R8P], loop_sizes=[4_I4P])
call half%cut_cells(x=h, y=h, cut=cut, moments_order=2_I4P)
call cut_links_2d(cut, links)                                                  ! 2-D: from the cut cells
call srd_neighbourhoods(links, cut%piece_cell, cut%piece_area, cut%piece_centroid, cut%piece_wet_normal, h, h, srd)
print '(A,I0,A,I0,A,I0)', '2-D: pieces ', cut%np, ', small ', count(srd%status == SRD_COMPLETE), ', regular ', &
                          count(srd%status == SRD_REGULAR)
p = minloc(cut%piece_area, 1, mask=cut%piece_area > 0._R8P .and. cut%piece_cell(2, :) == 5)
print '(A,2I3,A,F5.2,A,10I4)', '  piece in cell', cut%piece_cell(:, p), ', fraction ', cut%piece_area(p) / 0.01_R8P, &
                               ', SRD members', srd%member(srd%start(p):srd%start(p + 1) - 1)
print '(A,2I3,A,F6.3,A,F6.3)', '  member cells', cut%piece_cell(:, srd%member(srd%start(p) + 1)), &
                               ', eta_nb ', srd%eta_nb(p), ', eta_self ', srd%eta_self(p)
