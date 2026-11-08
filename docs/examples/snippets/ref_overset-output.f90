call build_connectivity(masks, donors, 3_I4P, 2_I4P, conn)
call save_connectivity(conn, 'cylinder.conn')
open(newunit=i, file='cylinder.conn', action='read')
do j = 1, 5
   read(i, '(A)') line ; print '(A)', line(1:min(len_trim(line), 100))
enddo
close(i)
