print '(A,I0,A,I0)', 'pieces in the block ', cut%np, ', cells with two or more ', count(cut%pieces > 1)
do j = 1, 48 ; do i = 1, 48
   if (cut%pieces(i, j) < 2) cycle
   p = cut%first_piece(i, j)                          ! its pieces are p, p + 1, ...
   print '(A,I0,A,I0,A,F7.4)', '  cell (', i, ',', j, ') alpha ', cut%alpha(i, j)
   do k = p, p + cut%pieces(i, j) - 1
      print '(A,I0,A,F7.4,A,4F7.4)', '    piece ', k, ': area / h^2 ', cut%piece_area(k) / h**2, &
                                     ', apertures', cut%piece_aperture(:, k)
   enddo
enddo ; enddo
