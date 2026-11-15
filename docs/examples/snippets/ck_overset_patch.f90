program ck_overset_patch
use fossil, only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, overset_patch_2d_t, overset_patch_2d, &
                   overlap_weights_t, overlap_weights_2d
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter         :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t) :: blocks(1)
type(topo_record_t)          :: rec
type(overset_patch_2d_t)     :: patch
type(overlap_weights_t)      :: w
real(R8P)                    :: u(33, 13), v(33, 13), x(41), r
integer(I4P)                 :: cells(2, 32), i, j, k, st

do j = 1, 13                                                  ! an O-grid around the unit cylinder, out to r = 2
   r = 1._R8P + (j - 1) / 12._R8P
   do i = 1, 33
      u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
call blocks(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)                                ! its outer boundary: the interface
x = [(-3._R8P + 0.15_R8P * (i - 1), i = 1, 41)]               ! a Cartesian background, edges along x and y
! A: the outer boundary cuts the background; fluxes cross the smallest faces
call overset_patch_2d(blocks, [1_I4P], x, x, patch, status=st)
print '(A,I0,A,I0,A,I0,A,I0)', 'status ', st, ': holes ', count(patch%hole), ', pieces ', patch%cut%np, &
                               ', smallest faces ', patch%ns
k = 1                                                         ! F . sf_area leaves the piece, enters the O-grid cell
print '(A,I0,A,2I3,A,2I3,A,2F8.4)', 'smallest face 1: piece ', patch%sf_piece(k), ' of cell', &
      patch%cut%piece_cell(:, patch%sf_piece(k)), ' -> O-grid cell', patch%sf_minor(2:3, k), ', sf_area', &
      patch%sf_area(:, k)
! B: the outer ring of O-grid cells receives from the background by intersection weights
cells(1, :) = [(i, i = 1, 32)] ; cells(2, :) = 12
call overlap_weights_2d(blocks(1), cells, x, x, w, status=st)
print '(A,I0,A,I0,A,I0,A,ES8.1)', 'status ', st, ': receptors ', w%nr, ', intersections ', w%nw, &
      ', max |sum w - 1| ', maxval([(abs(sum(w%weight(w%start(k):w%start(k + 1) - 1)) - 1._R8P), k = 1, w%nr)])
endprogram ck_overset_patch
