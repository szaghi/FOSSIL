do k = 1, 9 ! a wavy annulus 1 <= r <= 2, 32 x 8 cells around (clockwise) and out, 8 cells over 0 <= z <= 1
   do j = 1, 9
      do i = 1, 33
         t = -2 * PI * (i - 1) / 32 ; r = 1._R8P + (j - 1) / 8._R8P
         x(i, j, k) = r * cos(t) ; y(i, j, k) = r * sin(t)
         z(i, j, k) = (k - 1) / 8._R8P + 0.05_R8P * sin(2 * t) * sin(PI * (j - 1) / 8) * sin(PI * (k - 1) / 8)
      enddo
   enddo
enddo
x(33,:,:) = x(1,:,:) ; y(33,:,:) = y(1,:,:) ; z(33,:,:) = z(1,:,:)              ! the seam, exact
call grid(1)%initialize(x, y, z)
call detect_topology(grid, translations=reshape([0._R8P, 0._R8P, 1._R8P], [3, 1]))
do n = 1, grid(1)%nrec
   print '(A,I0,3A,I0,A,I0,A,3I3)', 'record ', n, ': ', KIND_NAME(grid(1)%rec(n)%kind), ' face ', &
         grid(1)%rec(n)%face, ' -> face ', grid(1)%rec(n)%donor_face, ', transform', grid(1)%rec(n)%transform
enddo
