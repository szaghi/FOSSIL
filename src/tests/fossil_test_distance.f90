!< FOSSIL, test distance computation.

program fossil_test_distance
!< FOSSIL, test distance computation.

use fossil, only : file_stl_object
use penf, only : I4P, R8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(file_stl_object)         :: file_stl            !< STL file.
type(vector_R8P), allocatable :: grid(:,:,:)         !< Grid.
real(R8P),        allocatable :: distance(:,:,:)     !< Distance of grid points to STL surface.
integer(I4P)                  :: ni, nj, nk          !< Grid dimensions.
integer(I4P)                  :: i, j, k             !< Counter.
real(R8P)                     :: Dx, Dy, Dz          !< Space steps.
type(vector_R8P)              :: emin                !< Minimum extent of the grid domain.
integer(I4P)                  :: file_unit           !< File unit.
logical                       :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call file_stl%initialize(file_name='src/tests/naca0012-binary.stl', is_ascii=.false.)
call file_stl%load_from_file

emin = -2._R8P * ex_R8P - 0.1_R8P * ey_R8P - 1._R8P * ez_R8P

ni = 32
nj = 32
nk = 32

Dx = 4._R8P  / ni ! x in [-2, 2]
Dy = 0.2_R8P / nj ! y in [-0.1, 0.1]
Dz = 4._R8P  / nk ! z in [-1, 3]

allocate(grid(1:ni, 1:nj, 1:nk))
allocate(distance(1:ni, 1:nj, 1:nk))

do k=1, nk
   do j=1, nj
      do i=1, ni
         grid(i, j, k) = emin + (i * Dx) * ex_R8P + (j * Dy) * ey_R8P + (k * Dz) * ez_R8P
      enddo
   enddo
enddo

print*, 'compute distances'
do k=1, nk
   do j=1, nj
      do i=1, ni
         distance(i, j, k) = file_stl%distance(point=grid(i, j, k))
      enddo
   enddo
enddo

print*, 'save output'
open(newunit=file_unit, file='fossil_test_distance.dat')
write(file_unit, '(A)')'VARIABLES = x y z distance'
write(file_unit, '(A)')'ZONE T="distance", I='//trim(str(ni))//', J='//trim(str(nj))//', K='//trim(str(nk))//''
do k=1, nk
   do j=1, nj
      do i=1, ni
         write(file_unit, '(A)') str(grid(i,j,k)%x)//' '//str(grid(i,j,k)%y)//' '//str(grid(i,j,k)%z)//' '//str(distance(i,j,k))
      enddo
   enddo
enddo
close(file_unit)

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
endprogram fossil_test_distance
