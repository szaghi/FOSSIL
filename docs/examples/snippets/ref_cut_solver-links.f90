x = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)           ! [0, 1]^3 on the cell faces
call cut_cells_3d(cube, x, x, x, cut3, links=links3, edges=.true.)             ! 3-D: links with the cut cells
print '(A,I0,A,I0,A,I0,A,I0)', '3-D: pieces ', cut3%np, ' (', count(cut3%piece_volume == 0._R8P), &
                               ' of zero volume), links ', size(links3%piece), ', unmatched regions ', links3%unmatched
p = 1
do k = links3%start(p), links3%start(p + 1) - 1 ! piece p meets piece(k) over its face face(k), area area(k)
   print '(A,I0,A,I4,A,I0,A,F7.4)', '  piece ', p, ' -> ', links3%piece(k), ' over face ', links3%face(k), &
                                    ', area ', links3%area(k)
enddo
print '(A,I0,A,I0)', '  z edges with fluid ', count(cut3%edge_z > 0._R8P), ' of ', size(cut3%edge_z)
