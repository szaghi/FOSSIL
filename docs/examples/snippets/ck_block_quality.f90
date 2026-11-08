program ck_block_quality
use fossil, only : curvilinear_block_2d_t, CELL_VALID, CELL_SINGULAR, CELL_INVALID
use penf,   only : I4P, R8P
implicit none
type(curvilinear_block_2d_t) :: grid
real(R8P)                    :: u(17, 9), v(17, 9)
real(R8P),    allocatable    :: skew(:,:), orth(:,:), aspect(:,:), minj(:,:)
integer(I4P), allocatable    :: valid(:,:)
integer(I4P)                 :: i, j

do j = 1, 9                                                   ! a 16 x 8 channel block over a bump
   do i = 1, 17
      u(i, j) = (i - 1) / 16._R8P
      v(i, j) = (j - 1) / 8._R8P * (0.5_R8P - 0.3_R8P * exp(-60._R8P * (u(i, j) - 0.5_R8P)**2)) + &
                0.3_R8P * exp(-60._R8P * (u(i, j) - 0.5_R8P)**2)
   enddo
enddo
u(5, 5) = u(5, 5) + 0.09_R8P                                  ! a vertex dragged past its neighbour: folded cells
call grid%initialize(u, v)
call grid%validity(valid)
print '(A,I0,A,I0,A,I0)', 'valid ', count(valid == CELL_VALID), ', singular ', count(valid == CELL_SINGULAR), &
                          ', invalid ', count(valid == CELL_INVALID)
call grid%quality(skewness=skew, orthogonality=orth, aspect=aspect, min_jacobian=minj)
print '(A,F6.3,A,F6.3,A,F6.2)', 'skewness max ', maxval(skew), ', orthogonality min ', minval(orth), &
                                ', aspect max ', maxval(aspect)
print '(A,F7.3)', 'scaled min J min ', minval(minj)
do j = 1, 8 ; do i = 1, 16
   if (valid(i, j) /= CELL_VALID) print '(A,2I3,A,F7.3,A,F6.3)', 'not valid: cell', i, j, ', min J ', minj(i, j), &
                                                                 ', skewness ', skew(i, j)
enddo ; enddo
endprogram ck_block_quality
