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

call file_stl%initialize(file_name='src/tests/naca0012-ascii.stl', is_ascii=.true.)
! call file_stl%initialize(file_name='src/tests/tubi-clean.stl', is_ascii=.false.)
call file_stl%load_from_file
call file_stl%compute_metrix

emin = -0.5_R8P * ex_R8P - 0.2_R8P * ey_R8P - 0.5_R8P * ez_R8P
! emin = -110._R8P * ex_R8P - 55._R8P * ey_R8P - 5._R8P * ez_R8P

ni = 64
nj = 64
nk = 64

Dx = 2._R8P  / ni ! x in [-0.5, 1.5]
Dy = 0.4_R8P / nj ! y in [-0.2, 0.2]
Dz = 2._R8P  / nk ! z in [-0.5, 1.5]
! Dx = 220._R8P / ni ! x in [-110, 110]
! Dy = 110._R8P / nj ! y in [-55, 55]
! Dz = 130._R8P / nk ! z in [-5, 125]

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
         ! distance(i, j, k) = file_stl%distance(point=grid(i, j, k), is_signed=.true., sign_algorithm='winding_number')
         ! distance(i, j, k) = file_stl%distance(point=grid(i, j, k), is_signed=.true., sign_algorithm='solid_angle')
         distance(i, j, k) = file_stl%distance(point=grid(i, j, k), is_signed=.true., sign_algorithm='ray_intersections')
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
