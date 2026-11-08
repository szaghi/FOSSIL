program ck_topology
use fossil, only : curvilinear_block_2d_t, detect_topology_2d
use penf,   only : I4P, R8P
implicit none
character(9), parameter      :: KIND(4) = ['boundary ', 'abut     ', 'periodic ', 'singular ']
type(curvilinear_block_2d_t) :: blocks(2)
real(R8P)                    :: u(9, 9), v(9, 9)
integer(I4P)                 :: b, r, i, j

do j = 1, 9 ; do i = 1, 9                       ! block 1: [0, 1] x [0, 1], i along x
   u(i, j) = (i - 1) / 8._R8P ; v(i, j) = (j - 1) / 8._R8P
enddo ; enddo
call blocks(1)%initialize(u, v)
do j = 1, 9 ; do i = 1, 9                       ! block 2: [1, 2] x [0, 1], i along -y, j along x
   u(i, j) = 1._R8P + (j - 1) / 8._R8P ; v(i, j) = 1._R8P - (i - 1) / 8._R8P
enddo ; enddo
call blocks(2)%initialize(u, v)
call detect_topology_2d(blocks, translations=reshape([0._R8P, 1._R8P, 2._R8P, 0._R8P], [2, 2])) ! periods (0, 1), (2, 0)
do b = 1, 2
   do r = 1, blocks(b)%nrec
      associate(rec => blocks(b)%rec(r))
         print '(A,I0,A,I0,3A,I0,A,I0,A,2I3)', 'block ', b, ' face ', rec%face, ': ', KIND(rec%kind), &
               'block ', rec%donor_block, ' face ', rec%donor_face, ', transform', rec%transform(1:2)
      endassociate
   enddo
enddo
endprogram ck_topology
