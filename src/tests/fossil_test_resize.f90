!< FOSSIL, test resize STL.

program fossil_test_resize
!< FOSSIL, test resize STL.

use flap, only : command_line_interface
use fossil, only : file_stl_object, surface_stl_object
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(file_stl_object)    :: file_stl            !< STL file.
type(surface_stl_object) :: surface             !< STL surface.
character(999)           :: file_name_stl       !< Input STL file name.
type(vector_R8P)         :: factor              !< Vectorial factor.
real(R8P)                :: x, y, z             !< Scalar factors.
logical                  :: are_tests_passed(4) !< Result of tests check.

are_tests_passed = .false.

call cli_parse
call file_stl%load_from_file(facet=surface%facet, file_name=trim(adjustl(file_name_stl)), guess_format=.true.)
call surface%analize
print '(A)', file_stl%statistics()
print '(A)', surface%statistics()

call surface%resize(factor=factor)
call file_stl%save_into_file(facet=surface%facet, file_name='fossil_test_resize-factor.stl')
are_tests_passed(1) = nint(surface%distance(point=2 * ex_R8P + 0 * ey_R8P + 0 * ez_R8P)) == 0
call surface%resize(factor=factor/4._R8P)
call surface%resize(x=x)
are_tests_passed(2) = nint(surface%distance(point=2 * ex_R8P + 0 * ey_R8P + 0 * ez_R8P)) == 0
call surface%resize(y=y)
are_tests_passed(3) = nint(surface%distance(point=2 * ex_R8P + 2 * ey_R8P + 0 * ez_R8P)) == 0
call surface%resize(z=z)
are_tests_passed(4) = nint(surface%distance(point=2 * ex_R8P + 2 * ey_R8P + 2 * ez_R8P)) == 0
call file_stl%save_into_file(facet=surface%facet, file_name='fossil_test_resize-xyz.stl')

call file_stl%load_from_file(facet=surface%facet, file_name=trim(adjustl(file_name_stl)), guess_format=.true.)

call surface%resize(factor=factor, respect_centroid=.true.)
call file_stl%save_into_file(facet=surface%facet, file_name='fossil_test_resize-factor-centroid.stl')

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli        !< Test command line interface.
  real(R8P)                    :: factor_(3) !< Vectorial factor.
  integer(I4P)                 :: error      !< Error trapping flag.

  call cli%init(progname='fossil_test_resize',                              &
                authors='S. Zaghi',                                         &
                help='Usage: ',                                             &
                examples=["fossil_test_resize --stl src/tests/dragon.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--stl',               &
               help='STL (input) file name', &
               required=.false.,             &
               def='src/tests/cube.stl',     &
               act='store')

  call cli%add(switch='--factor',       &
               help='vectorial factor', &
               required=.false.,        &
               nargs='+',               &
               def='2.0 2.0 2.0',       &
               act='store')

  call cli%add(switch='--x',     &
               help='factor x',  &
               required=.false., &
               def='2.0',        &
               act='store')

  call cli%add(switch='--y',     &
               help='factor y',  &
               required=.false., &
               def='2.0',        &
               act='store')

  call cli%add(switch='--z',     &
               help='factor z',  &
               required=.false., &
               def='2.0',        &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl',    val=file_name_stl, error=error) ; if (error/=0) stop
  call cli%get(switch='--factor', val=factor_,       error=error) ; if (error/=0) stop
  call cli%get(switch='--x',      val=x,             error=error) ; if (error/=0) stop
  call cli%get(switch='--y',      val=y,             error=error) ; if (error/=0) stop
  call cli%get(switch='--z',      val=z,             error=error) ; if (error/=0) stop
  factor%x = factor_(1)
  factor%y = factor_(2)
  factor%z = factor_(3)
  endsubroutine cli_parse
endprogram fossil_test_resize
