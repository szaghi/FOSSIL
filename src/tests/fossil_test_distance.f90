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
integer(I4P)                  :: file_unit           !< File unit.
logical                       :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call file_stl%initialize(file_name='src/tests/naca0012-binary.stl')
call file_stl%initialize(file_name='src/tests/dragon.stl')
call file_stl%load_from_file(guess_format=.true.)
call file_stl%compute_metrix
call file_stl%create_aabb_tree(refinement_levels=2)
call file_stl%aabb%save_geometry_tecplot_ascii(file_name='fossil_test_distance_aabb_tree.dat')
stop

associate(bmin=>file_stl%aabb%node(0)%bmin(), bmax=>file_stl%aabb%node(0)%bmax())
   ni = 64
   nj = 64
   nk = 64

   Dx = (bmax%x - bmin%x) / ni
   Dy = (bmax%y - bmin%y) / nj
   Dz = (bmax%z - bmin%z) / nk

   allocate(grid(-4:ni+5, -4:nj+5, -4:nk+5))
   allocate(distance(-4:ni+5, -4:nj+5, -4:nk+5))

   do k=-4, nk + 5
      do j=-4, nj + 5
         do i=-4, ni + 5
            grid(i, j, k) = bmin + (i * Dx) * ex_R8P + (j * Dy) * ey_R8P + (k * Dz) * ez_R8P
         enddo
      enddo
   enddo
endassociate

print*, 'compute distances'
do k=-4, nk + 5
   do j=-4, nj + 5
      do i=-4, ni + 5
         distance(i, j, k) = file_stl%distance(point=grid(i, j, k), is_signed=.true., sign_algorithm='ray_intersections')
      enddo
   enddo
enddo

print*, 'save output'
open(newunit=file_unit, file='fossil_test_distance.dat')
write(file_unit, '(A)')'VARIABLES = x y z distance'
write(file_unit, '(A)')'ZONE T="distance", I='//trim(str(ni+10))//', J='//trim(str(nj+10))//', K='//trim(str(nk+10))//''
do k=-4, nk + 5
   do j=-4, nj + 5
      do i=-4, ni + 5
         write(file_unit, '(A)') str(grid(i,j,k)%x)//' '//str(grid(i,j,k)%y)//' '//str(grid(i,j,k)%z)//' '//str(distance(i,j,k))
      enddo
   enddo
enddo
close(file_unit)

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
endprogram fossil_test_distance
