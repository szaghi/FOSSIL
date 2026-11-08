print '(A,I0,A,I0)', 'pieces ', cut%np, ', cells with two or more ', count(cut%pieces > 1)
c = findloc(cut%pieces, 2)                                                         ! the first cell with two pieces
p = cut%first_piece(c(1), c(2), c(3))
print '(A,I0,A,I0,A,I0,A,F8.5,A,2F8.5)', 'cell (', c(1), ',', c(2), ',', c(3), ') alpha ', &
   cut%alpha(c(1), c(2), c(3)), ' = ', cut%piece_volume(p:p + 1) / h**3
print '(A,6F7.3)', '  faces of piece 1 (-x +x -y +y -z +z)', cut%piece_face(:, p)
print '(A,6F7.3)', '  faces of piece 2                    ', cut%piece_face(:, p + 1)
print '(A,6F7.3)', '  faces of the cell                   ', cut%face_x(c(1):c(1) + 1, c(2), c(3)), &
                   cut%face_y(c(1), c(2):c(2) + 1, c(3)), cut%face_z(c(1), c(2), c(3):c(3) + 1)
