!< FOSSIL test: cell-centre blocks (issue #72, step 10).
program fossil_test_overset_cells
!< FOSSIL test: cell-centre blocks (issue #72, step 10).
!<
!< 1. the cell-centre block of a stretched Cartesian block: interior vertices at the cell centres, the boundary layer
!<    at the face centres, edge midpoints and corners; every vertex of the original block lies in it (it covers the
!<    block); 3-D the same on a sheared block.

use fossil, only : curvilinear_block_t, curvilinear_block_2d_t, cell_centre_block, cell_centre_block_2d, LOC_INSIDE
use penf,   only : I4P, R8P

implicit none

type(curvilinear_block_2d_t) :: b2, d2
type(curvilinear_block_t)    :: b3, d3
real(R8P)                    :: u(6, 5), v(6, 5), x(4, 3, 3), y(4, 3, 3), z(4, 3, 3), err
real(R8P),    allocatable    :: xi(:,:)
integer(I4P), allocatable    :: cell(:,:), st(:)
integer(I4P)                 :: i, j, k
logical                      :: all_ok

all_ok = .true.
do j = 1, 5
   do i = 1, 6
      u(i, j) = (i - 1) ** 1.5_R8P ; v(i, j) = 0.3_R8P * (j - 1)
   enddo
enddo
call b2%initialize(u, v) ; call cell_centre_block_2d(b2, d2)
err = max(abs(d2%u(3, 4) - 0.5_R8P * (u(2, 3) + u(3, 3))), abs(d2%v(3, 4) - 0.5_R8P * (v(2, 3) + v(2, 4))), &
          abs(d2%u(1, 3) - u(1, 2)), abs(d2%v(1, 3) - 0.5_R8P * (v(1, 2) + v(1, 3))), abs(d2%u(7, 6) - u(6, 5)))
call d2%build_search
call d2%locate_many(reshape(u, [30]), reshape(v, [30]), cell, xi, st)
call check(d2%ni == 7 .and. d2%nj == 6 .and. err <= 1e-15_R8P .and. all(st == LOC_INSIDE), &
           '2-D cell-centre block: centres, boundary layer, covers the block', err)
do k = 1, 3
   do j = 1, 3
      do i = 1, 4
         x(i, j, k) = 0.2_R8P * (i - 1) + 0.05_R8P * j ; y(i, j, k) = 0.3_R8P * (j - 1) ; z(i, j, k) = 0.25_R8P * (k - 1)
      enddo
   enddo
enddo
call b3%initialize(x, y, z) ; call cell_centre_block(b3, d3) ; call d3%build_search
call d3%locate_many(reshape(x, [36]), reshape(y, [36]), reshape(z, [36]), cell, xi, st)
err = abs(d3%x(2, 2, 2) - sum(x(1:2, 1:2, 1:2)) / 8)
call check(d3%ni == 5 .and. err <= 1e-15_R8P .and. all(st == LOC_INSIDE), '3-D cell-centre block covers the block', err)

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES11.3)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check
endprogram fossil_test_overset_cells
