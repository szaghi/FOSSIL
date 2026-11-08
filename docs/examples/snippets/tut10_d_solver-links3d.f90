call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call cut_cells_3d(surface=surface, x=x(1:49:2), y=x(1:49:2), z=x(1:49:2), cut=cut3, links=links3)
print '(A,I0,A,I0,A,I0)', '3-D: pieces ', cut3%np, ', links ', size(links3%piece), ', unmatched regions ', &
                          links3%unmatched
