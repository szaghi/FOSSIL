call grid(1)%geometry(volume=vol, face_i=fi, face_j=fj, face_k=fk)
gcl = 0._R8P ! the six outward face vectors of each cell sum to zero
do k = 1, 8
   do j = 1, 8
      do i = 1, 32
         gcl = max(gcl, maxval(abs(fi(:, i + 1, j, k) - fi(:, i, j, k) + fj(:, i, j + 1, k) - fj(:, i, j, k) + &
                                   fk(:, i, j, k + 1) - fk(:, i, j, k))))
      enddo
   enddo
enddo
print '(A,F12.9,A,ES9.2)', 'volume ', sum(vol), ', largest face-vector sum of a cell ', gcl
call grid(1)%validity(valid)
print '(A,I0,A,I0)', 'cells certified valid ', count(valid == CELL_VALID), ' of ', size(valid)
