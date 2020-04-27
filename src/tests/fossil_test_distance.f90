!< FOSSIL, test distance computation.

program fossil_test_distance
!< FOSSIL, test distance computation.

use flap, only : command_line_interface
use fossil, only : file_stl_object, surface_stl_object
use penf, only : I4P, I8P, R8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P
use fossil_aabb_tree_object, only : aabb_tree_object

implicit none

type(file_stl_object)         :: file_stl                !< STL file.
type(surface_stl_object)      :: surface_stl             !< STL surface.
type(vector_R8P), allocatable :: grid(:,:,:)             !< Grid.
real(R8P),        allocatable :: distance(:,:,:)         !< Distance of grid points to STL surface.
character(999)                :: file_name_stl           !< Input STL file name.
integer(I4P)                  :: refinement_levels       !< AABB refinement levels used.
logical                       :: save_aabb_tree_geometry !< Sentinel to save AABB geometry.
logical                       :: save_aabb_tree_stl      !< Sentinel to save AABB stl.
logical                       :: test_brute_force        !< Sentinel to test also brute force.
character(999)                :: sign_algorithm          !< Algorithm used for "point in polyhedron" test.
logical                       :: unsigned                !< Compute unsigned distance.
integer(I4P)                  :: ni, nj, nk, gi, gj, gk  !< Grid dimensions.
integer(I4P)                  :: i, j, k                 !< Counter.
real(R8P)                     :: Dx, Dy, Dz              !< Space steps.
integer(I4P)                  :: file_unit               !< File unit.
integer(I8P)                  :: timing(0:4)             !< Tic toc timing.
logical                       :: are_tests_passed(1)     !< Result of tests check.

are_tests_passed = .false.

call cli_parse
call file_stl%load_from_file(facet=surface_stl%facet, file_name=trim(adjustl(file_name_stl)), guess_format=.true.)
 call surface_stl%analize(aabb_refinement_levels=refinement_levels)
print '(A)', 'STL statistics before sanitization'
print '(A)', file_stl%statistics()
print '(A)', surface_stl%statistics()
call surface_stl%sanitize
call surface_stl%analize(aabb_refinement_levels=refinement_levels)
print '(A)', 'STL statistics after sanitization'
print '(A)', surface_stl%statistics()

are_tests_passed = int(surface_stl%distance(point=0*ex_R8P), I4P) == 0_I4P

if (save_aabb_tree_geometry) call surface_stl%aabb%save_geometry_tecplot_ascii(file_name='fossil_test_distance_aabb_tree.dat')
if (save_aabb_tree_stl) call file_stl%save_aabb_into_file(surface=surface_stl, base_file_name='fossil_test_distance_', &
                                                          is_ascii=.false.)

! stop

associate(bmin=>surface_stl%bmin, bmax=>surface_stl%bmax)
   Dx = (bmax%x - bmin%x) / ni
   Dy = (bmax%y - bmin%y) / nj
   Dz = (bmax%z - bmin%z) / nk

   allocate(grid(1-gi:ni+gi, 1-gj:nj+gj, 1-gk:nk+gk))
   allocate(distance(1-gi:ni+gi, 1-gj:nj+gj, 1-gk:nk+gk))

   do k=1 - gk, nk + gk
      do j=1 - gj, nj + gj
         do i=1 - gi, ni + gi
            grid(i, j, k) = bmin + (i * Dx) * ex_R8P + (j * Dy) * ey_R8P + (k * Dz) * ez_R8P
         enddo
      enddo
   enddo
endassociate

if (test_brute_force) then
   surface_stl%aabb%is_initialized = .false.
   print '(A)', 'compute distances brute force'
   call system_clock(timing(1))
   do k=1 - gk, nk + gk
      do j=1 - gj, nj + gj
         do i=1 - gi, ni + gi
            distance(i, j, k) = surface_stl%distance(point=grid(i, j, k), is_signed=.true., sign_algorithm=trim(sign_algorithm))
         enddo
      enddo
   enddo
   call system_clock(timing(2), timing(0))
   print '(A, F8.3)', 'brute force timing: ', real(timing(2) - timing(1))/ timing(0)

   print '(A)', 'save output'
   open(newunit=file_unit, file='fossil_test_distance-brute.dat')
   write(file_unit, '(A)')'VARIABLES = x y z distance'
   write(file_unit, '(A)')'ZONE T="distance", I='//trim(str(ni+10))//', J='//trim(str(nj+10))//', K='//trim(str(nk+10))//''
   do k=1 - gk, nk + gk
      do j=1 - gj, nj + gj
         do i=1 - gi, ni + gi
            write(file_unit, '(A)') str(grid(i,j,k)%x)//' '//str(grid(i,j,k)%y)//' '//str(grid(i,j,k)%z)//' '//str(distance(i,j,k))
         enddo
      enddo
   enddo
   close(file_unit)
endif

surface_stl%aabb%is_initialized = .true.
print '(A)', 'compute distances AABB'
call system_clock(timing(3))
do k=1 - gk, nk + gk
   do j=1 - gj, nj + gj
      do i=1 - gi, ni + gi
         distance(i, j, k) = surface_stl%distance(point=grid(i, j, k), is_signed=.not.unsigned, sign_algorithm=trim(sign_algorithm))
      enddo
   enddo
enddo
call system_clock(timing(4), timing(0))
print '(A, F8.3)', 'AABB timing: ', real(timing(4) - timing(3))/ timing(0)

print '(A)', 'save output'
open(newunit=file_unit, file='fossil_test_distance-aabb.dat')
write(file_unit, '(A)')'VARIABLES = x y z distance'
write(file_unit, '(A)')'ZONE T="distance", I='//trim(str(ni+10))//', J='//trim(str(nj+10))//', K='//trim(str(nk+10))//''
do k=1 - gk, nk + gk
   do j=1 - gj, nj + gj
      do i=1 - gi, ni + gi
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

  call cli%add(switch='--ni',                      &
               help='cells number in i direction', &
               required=.false.,                   &
               def='128',                          &
               act='store')

  call cli%add(switch='--nj',                      &
               help='cells number in j direction', &
               required=.false.,                   &
               def='128',                          &
               act='store')

  call cli%add(switch='--nk',                      &
               help='cells number in k direction', &
               required=.false.,                   &
               def='3',                            &
               act='store')

  call cli%add(switch='--gi',                            &
               help='ghost cells number in i direction', &
               required=.false.,                         &
               def='5',                                  &
               act='store')

  call cli%add(switch='--gj',                            &
               help='ghost cells number in j direction', &
               required=.false.,                         &
               def='5',                                  &
               act='store')

  call cli%add(switch='--gk',                            &
               help='ghost cells number in k direction', &
               required=.false.,                         &
               def='5',                                  &
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

  call cli%add(switch='--unsigned',              &
               help='compute unsigned distance', &
               required=.false.,                 &
               def='.false.',                    &
               act='store_true')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl',                     val=file_name_stl,           error=error) ; if (error/=0) stop
  call cli%get(switch='--ni',                      val=ni,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--nj',                      val=nj,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--nk',                      val=nk,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--gi',                      val=gi,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--gj',                      val=gj,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--gk',                      val=gk,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--ref_levels',              val=refinement_levels,       error=error) ; if (error/=0) stop
  call cli%get(switch='--save_aabb_tree_geometry', val=save_aabb_tree_geometry, error=error) ; if (error/=0) stop
  call cli%get(switch='--save_aabb_tree_stl',      val=save_aabb_tree_stl,      error=error) ; if (error/=0) stop
  call cli%get(switch='--brute_force',             val=test_brute_force,        error=error) ; if (error/=0) stop
  call cli%get(switch='--sign_algorithm',          val=sign_algorithm,          error=error) ; if (error/=0) stop
  call cli%get(switch='--unsigned',                val=unsigned,                error=error) ; if (error/=0) stop
  endsubroutine cli_parse
endprogram fossil_test_distance
