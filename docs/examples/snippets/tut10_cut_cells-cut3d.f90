call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call cut_cells_3d(surface=surface, x=x(1:49:2), y=x(1:49:2), z=x(1:49:2), cut=cut3)    ! 24^3 cells, h = 0.05
print '(A,I0,A,I0,A,I0)', '3-D: cut cells ', count(cut3%alpha > 0._R8P .and. cut3%alpha < 1._R8P), &
                          ', flagged ', count(cut3%pieces < 0), ', cells with 2+ pieces ', count(cut3%pieces > 1)
c = findloc(cut3%pieces, 2)                                   ! the first cell with two fluid pieces
p = cut3%first_piece(c(1), c(2), c(3))
print '(A,I0,A,I0,A,I0,A,F8.5,A,2F8.5)', '     cell (', c(1), ',', c(2), ',', c(3), ') alpha ', &
   cut3%alpha(c(1), c(2), c(3)), '  pieces ', (cut3%piece_volume(p + k) / 0.05_R8P**3, k = 0, 1)
