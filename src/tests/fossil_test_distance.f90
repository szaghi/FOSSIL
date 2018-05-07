!< FOSSIL, test distance computation.

program fossil_test_distance
!< FOSSIL, test distance computation.

use flap, only : command_line_interface
use fossil, only : file_stl_object
use penf, only : I4P, I8P, R8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(file_stl_object)         :: file_stl                !< STL file.
type(vector_R8P), allocatable :: grid(:,:,:)             !< Grid.
real(R8P),        allocatable :: distance(:,:,:)         !< Distance of grid points to STL surface.
character(999)                :: file_name_stl           !< Input STL file name.
integer(I4P)                  :: refinement_levels       !< AABB refinement levels used.
logical                       :: save_aabb_tree_geometry !< Sentinel to save AABB geometry.
logical                       :: save_aabb_tree_stl      !< Sentinel to save AABB stl.
logical                       :: test_brute_force        !< Sentinel to test also brute force.
character(999)                :: sign_algorithm          !< Algorithm used for "point in polyhedron" test.
integer(I4P)                  :: ni, nj, nk              !< Grid dimensions.
integer(I4P)                  :: i, j, k                 !< Counter.
real(R8P)                     :: Dx, Dy, Dz              !< Space steps.
integer(I4P)                  :: file_unit               !< File unit.
integer(I8P)                  :: timing(0:4)             !< Tic toc timing.
logical                       :: are_tests_passed(1)     !< Result of tests check.

are_tests_passed = .false.

call cli_parse
call file_stl%initialize(file_name=trim(adjustl(file_name_stl)))
call file_stl%load_from_file(guess_format=.true.)

call file_stl%build_connectivity
call file_stl%sanitize_normals

call file_stl%compute_metrix
call file_stl%create_aabb_tree(refinement_levels=refinement_levels)

are_tests_passed = int(file_stl%distance(point=0*ex_R8P), I4P) == 0_I4P

if (save_aabb_tree_geometry) call file_stl%aabb%save_geometry_tecplot_ascii(file_name='fossil_test_distance_aabb_tree.dat')
if (save_aabb_tree_stl) call file_stl%aabb%save_into_file_stl(base_file_name='fossil_test_distance_', is_ascii=.true.)

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

if (test_brute_force) then
   file_stl%aabb%is_initialized = .false.
   print*, 'compute distances brute force'
   call system_clock(timing(1))
   do k=-4, nk + 5
      do j=-4, nj + 5
         do i=-4, ni + 5
            distance(i, j, k) = file_stl%distance(point=grid(i, j, k), is_signed=.true., sign_algorithm=trim(sign_algorithm))
         enddo
      enddo
   enddo
   call system_clock(timing(2), timing(0))
   print*, 'brute force timing: ', real(timing(2) - timing(1))/ timing(0)

   print*, 'save output'
   open(newunit=file_unit, file='fossil_test_distance-brute.dat')
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
endif

file_stl%aabb%is_initialized = .true.
print*, 'compute distances AABB'
call system_clock(timing(3))
do k=-4, nk + 5
   do j=-4, nj + 5
      do i=-4, ni + 5
         distance(i, j, k) = file_stl%distance(point=grid(i, j, k), is_signed=.true., sign_algorithm=trim(sign_algorithm))
      enddo
   enddo
enddo
call system_clock(timing(4), timing(0))
print*, 'AABB timing: ', real(timing(4) - timing(3))/ timing(0)

print*, 'save output'
open(newunit=file_unit, file='fossil_test_distance-aabb.dat')
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
contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli   !< Test command line interface.
  integer(I4P)                 :: error !< Error trapping flag.

  call cli%init(progname='fossil_test_distance',                              &
                authors='S. Zaghi',                                           &
                help='Usage: ',                                               &
                examples=["fossil_test_distance --stl src/tests/dragon.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--stl',                      &
               help='STL (input) file name',        &
               required=.false.,                    &
               def='src/tests/naca0012-binary.stl', &
               act='store')

  call cli%add(switch='--ref_levels',         &
               help='AABB refinement levels', &
               required=.false.,              &
               def='2',                       &
               act='store')

  call cli%add(switch='--save_aabb_tree_geometry', &
               help='save AABB tree geometry',     &
               required=.false.,                   &
               def='.false.',                      &
               act='store_true')

  call cli%add(switch='--save_aabb_tree_stl',  &
               help='save AABB tree STL',      &
               required=.false.,               &
               def='.false.',                  &
               act='store_true')

  call cli%add(switch='--brute_force',         &
               help='test (also) brute force', &
               required=.false.,               &
               def='.false.',                  &
               act='store_true')

  call cli%add(switch='--sign_algorithm',                         &
               help='algorithm used to compute sign of distance', &
               required=.false.,                                  &
               def='ray_intersections',                           &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl',                     val=file_name_stl,           error=error) ; if (error/=0) stop
  call cli%get(switch='--ref_levels',              val=refinement_levels,       error=error) ; if (error/=0) stop
  call cli%get(switch='--save_aabb_tree_geometry', val=save_aabb_tree_geometry, error=error) ; if (error/=0) stop
  call cli%get(switch='--save_aabb_tree_stl',      val=save_aabb_tree_stl,      error=error) ; if (error/=0) stop
  call cli%get(switch='--brute_force',             val=test_brute_force,        error=error) ; if (error/=0) stop
  call cli%get(switch='--sign_algorithm',          val=sign_algorithm,          error=error) ; if (error/=0) stop
  endsubroutine cli_parse
endprogram fossil_test_distance
