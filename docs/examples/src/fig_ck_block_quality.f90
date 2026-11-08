!run figck-block-quality fig_ck_block_quality
!fig ck-block-quality
!< Figure of the cookbook recipe "Check a block's validity and quality": the 16 x 8 channel block over a bump of
!< ck_block_quality, with one vertex dragged past its neighbour, seen from above. Each cell is coloured by its
!< equiangle skewness (0 square, 1 degenerate); the cells validity does not certify are outlined in ember.
program fig_ck_block_quality
use fossil,     only : curvilinear_block_2d_t, CELL_VALID
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(curvilinear_block_2d_t) :: grid
real(R8P)                    :: u(17, 9), v(17, 9)
real(R8P),    allocatable    :: skew(:,:)
integer(I4P), allocatable    :: valid(:,:)
type(vector_R8P)             :: c(4)
integer(I4P)                 :: i, j, un
character(100)               :: caption

do j = 1, 9
   do i = 1, 17
      u(i, j) = (i - 1) / 16._R8P
      v(i, j) = (j - 1) / 8._R8P * (0.5_R8P - 0.3_R8P * exp(-60._R8P * (u(i, j) - 0.5_R8P)**2)) + &
                0.3_R8P * exp(-60._R8P * (u(i, j) - 0.5_R8P)**2)
   enddo
enddo
u(5, 5) = u(5, 5) + 0.09_R8P
call grid%initialize(u, v)
call grid%validity(valid)
call grid%quality(skewness=skew)
write(caption, '(A,F5.3,A,I0,A)') 'cells coloured by skewness (max ', maxval(skew), '); ', &
                                  count(valid /= CELL_VALID), ' folded cells outlined'
call fig_open(un, 'ck-block-quality', trim(caption), -90._R8P, 90._R8P, height=360_I4P)
do j = 1, 8 ; do i = 1, 16
   call fig_quad(un, 'frac', p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1), value=skew(i, j))
enddo ; enddo
do j = 1, 8 ; do i = 1, 16
   c = [p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1)]
   call fig_segments(un, merge('ghost', 'hole ', valid(i, j) == CELL_VALID), c, [c(2), c(3), c(4), c(1)])
enddo ; enddo
call fig_points(un, 'accent', [p(5, 5)])
call fig_close(un)
contains
   function p(i, j) result(q)
   !< Vertex (i, j) as a point of the plane z = 0.
   integer(I4P), intent(in) :: i, j !< Vertex.
   type(vector_R8P)         :: q    !< Point.

   q = vector_R8P(u(i, j), v(i, j), 0._R8P)
   endfunction p
endprogram fig_ck_block_quality
