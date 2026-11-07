program ck_block_locate
use fossil, only : curvilinear_block_2d_t, detect_topology_2d, CELL_VALID, LOC_INSIDE
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter         :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t) :: grid(1)
real(R8P)                    :: u(33, 9), v(33, 9), r, f
real(R8P),    allocatable    :: xi(:,:), w(:,:)
integer(I4P), allocatable    :: cell(:,:), st(:), nodes(:,:,:), sst(:), valid(:,:)
integer(I4P)                 :: i, j, k

do j = 1, 9                                                   ! an O-grid around the unit cylinder
   r = 1._R8P + 0.15_R8P * (j - 1)
   do i = 1, 33
      u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)                           ! the seam, written exactly
call grid(1)%initialize(u, v)
call detect_topology_2d(grid)                                 ! records: the seam on faces 1 and 2
print '(A,I0,A,2(I0,1X))', 'records ', grid(1)%nrec, ', faces ', grid(1)%rec(:)%face
call grid(1)%validity(valid)
print '(A,I0,A,I0)', 'valid cells ', count(valid == CELL_VALID), ' of ', size(valid)
call grid(1)%build_search
call grid(1)%locate_many([1.62_R8P], [0.03_R8P], cell, xi, st)  ! a point next to the seam
print '(A,L1,A,2I3,A,2F9.5)', 'inside ', st(1) == LOC_INSIDE, ', cell', cell(:, 1), ', xi', xi(:, 1)
call grid(1)%stencils(cell, xi, 4_I4P, nodes, w, sst)          ! width 4, order-matched xi
print '(A,16I3)', 'stencil i ', nodes(1, :, 1)
f = 0._R8P
do k = 1, 16
   f = f + w(k, 1) * exp(u(nodes(1, k, 1), nodes(2, k, 1))) * cos(v(nodes(1, k, 1), nodes(2, k, 1)))
enddo
print '(A,ES10.2)', 'interpolation error of exp(x) cos(y) ', abs(f - exp(1.62_R8P) * cos(0.03_R8P))
endprogram ck_block_locate
