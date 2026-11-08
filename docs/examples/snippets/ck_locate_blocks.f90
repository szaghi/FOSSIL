program ck_locate_blocks
use fossil, only : curvilinear_block_t, locate_in_blocks
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
type(curvilinear_block_t) :: blocks(2)
integer(I4P), allocatable :: block(:), cell(:,:), status(:)
real(R8P),    allocatable :: xi(:,:)
real(R8P)                 :: x(25, 13, 3), y(25, 13, 3), z(25, 13, 3), r, th, ph, g(9)
real(R8P)                 :: bx(9, 9, 9), by(9, 9, 9), bz(9, 9, 9)
integer(I4P)              :: i, j, k

do k = 1, 3 ; do j = 1, 13 ; do i = 1, 25                   ! block 1: a shell 1 <= r <= 1.5
   r = 1._R8P + 0.25_R8P * (k - 1) ; th = 2 * PI * (i - 1) / 24 ; ph = PI * (j - 1) / 12
   x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
enddo ; enddo ; enddo
x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:) ; x(:,[1,13],:) = 0._R8P ; y(:,[1,13],:) = 0._R8P
call blocks(1)%initialize(x, y, z)
g = [(-1.6_R8P + 0.4_R8P * (i - 1), i = 1, 9)]              ! block 2: a Cartesian box [-1.6, 1.6]^3
do k = 1, 9 ; do j = 1, 9
   bx(:, j, k) = g ; by(:, j, k) = g(j) ; bz(:, j, k) = g(k)
enddo ; enddo
call blocks(2)%initialize(bx, by, bz)
call blocks(1)%build_search ; call blocks(2)%build_search
! in the box only (inside the sphere), in both (the first block wins), in the box's corner only, in none
call locate_in_blocks(blocks, px=[0.3_R8P, 1.2_R8P, 1.5_R8P, 2._R8P], py=[0.2_R8P, 0.3_R8P, 1.5_R8P, 2._R8P], &
                      pz=[0.1_R8P, 0.2_R8P, 0._R8P, 2._R8P], block=block, cell=cell, xi=xi, status=status)
do i = 1, 4
   print '(A,I0,A,I0,A,3I3,A,3F7.4,A,I0)', 'point ', i, ': block ', block(i), ', cell', cell(:, i), &
                                         ', xi', xi(:, i), ', status ', status(i)
enddo
endprogram ck_locate_blocks
