!run figck-topology fig_ck_topology
!fig ck-topology
!< Figure of the cookbook recipe "Find how blocks connect": the two 8 x 8 blocks of ck_topology, seen from above, block
!< 2 indexed with i along -y. The records detect_topology_2d proposes are drawn on the faces: the abutting interface
!< in ember, the periodic faces (period (0, 1) and (2, 0)) in yellow; teal lines show each block's i and j directions.
program fig_ck_topology
use fossil,     only : curvilinear_block_2d_t, detect_topology_2d, TOPO_ABUT
use fossil_fig, only : fig_open, fig_segment, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(curvilinear_block_2d_t) :: blocks(2)
real(R8P)                    :: u(9, 9), v(9, 9)
type(vector_R8P)             :: a(200), c(200), o
integer(I4P)                 :: b, r, i, j, n, un

do j = 1, 9 ; do i = 1, 9
   u(i, j) = (i - 1) / 8._R8P ; v(i, j) = (j - 1) / 8._R8P
enddo ; enddo
call blocks(1)%initialize(u, v)
do j = 1, 9 ; do i = 1, 9
   u(i, j) = 1._R8P + (j - 1) / 8._R8P ; v(i, j) = 1._R8P - (i - 1) / 8._R8P
enddo ; enddo
call blocks(2)%initialize(u, v)
call detect_topology_2d(blocks, translations=reshape([0._R8P, 1._R8P, 2._R8P, 0._R8P], [2, 2]))
call fig_open(un, 'ck-topology', 'detect_topology_2d: abutting interface (ember), periodic faces (yellow)', &
              -90._R8P, 90._R8P, height=380_I4P)
do b = 1, 2
   n = 0
   do j = 1, 9 ; do i = 1, 8
      n = n + 1 ; a(n) = p(b, i, j) ; c(n) = p(b, i + 1, j)
      n = n + 1 ; a(n) = p(b, j, i) ; c(n) = p(b, j, i + 1)
   enddo ; enddo
   call fig_segments(un, 'ghost', a(1:n), c(1:n))
   do r = 1, blocks(b)%nrec                                    ! the record's face, its vertex range
      associate(rec => blocks(b)%rec(r))
         call fig_segment(un, merge('hole  ', 'accent', rec%kind == TOPO_ABUT), &
                          p(b, rec%lo(1), rec%lo(2), 1e-3_R8P), p(b, rec%hi(1), rec%hi(2), 1e-3_R8P))
      endassociate
   enddo
   o = p(b, 2, 2, 2e-3_R8P)                                    ! index directions
   call fig_segment(un, 'probe', o, p(b, 4, 2, 2e-3_R8P)) ; call fig_label(un, 'probe', p(b, 4, 2, 2e-3_R8P), 'i')
   call fig_segment(un, 'probe', o, p(b, 2, 4, 2e-3_R8P)) ; call fig_label(un, 'probe', p(b, 2, 4, 2e-3_R8P), 'j')
   call fig_label(un, 'accent', p(b, 5, 5, 2e-3_R8P), merge('block 1', 'block 2', b == 1))
enddo
call fig_close(un)
contains
   function p(b, i, j, dz) result(q)
   !< Vertex (i, j) of block b, lifted by dz (optional) towards the viewer.
   integer(I4P), intent(in)           :: b, i, j !< Block, vertex.
   real(R8P),    intent(in), optional :: dz      !< Lift.
   type(vector_R8P)                   :: q       !< Point.

   q = vector_R8P(blocks(b)%u(i, j), blocks(b)%v(i, j), 0._R8P)
   if (present(dz)) q%z = dz
   endfunction p
endprogram fig_ck_topology
