call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%sanitize
strict = dragon
call dragon%fill_holes(holes_filled=filled, open_edges=open_edges, nonmanifold_edges=nonmanifold)
print '(A,I0,A,I0,A,I0)', 'dragon, default:              filled ', filled, ', open edges ', open_edges, &
                          ', non-manifold edges ', nonmanifold
call strict%fill_holes(holes_filled=filled, open_edges=open_edges, nonmanifold_edges=nonmanifold, manifold_only=.true.)
print '(A,I0,A,I0,A,I0)', 'dragon, manifold_only=.true.: filled ', filled, ', open edges ', open_edges, &
                          ', non-manifold edges ', nonmanifold
