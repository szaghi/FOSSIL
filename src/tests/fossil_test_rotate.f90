!< FOSSIL, test rotate STL.

program fossil_test_rotate
!< FOSSIL, test rotate STL.

use flap, only : command_line_interface
use fossil, only : file_stl_object, surface_stl_object
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(file_stl_object)    :: file_stl            !< STL file.
type(surface_stl_object) :: surface             !< STL surface.
character(999)           :: file_name_stl       !< Input STL file name.
type(vector_R8P)         :: axis                !< Axis of rotation.
real(R8P)                :: angle               !< Angle of rotation.
logical                  :: are_tests_passed(2) !< Result of tests check.

are_tests_passed = .false.

call cli_parse
call file_stl%load_from_file(facet=surface%facet, file_name=trim(adjustl(file_name_stl)), guess_format=.true.)
call surface%analize

call surface%rotate(axis=axis, angle=angle)
call file_stl%save_into_file(facet=surface%facet, file_name='fossil_test_rotate.stl')
are_tests_passed(1) = nint(surface%distance(point=0 * ex_R8P - 1 * ey_R8P + 1 * ez_R8P)) == 0
are_tests_passed(2) = nint(surface%distance(point=1 * ex_R8P - 1 * ey_R8P + 0 * ez_R8P)) == 0

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli      !< Test command line interface.
  real(R8P)                    :: axis_(3) !< Axis of rotation.
  integer(I4P)                 :: error    !< Error trapping flag.

  call cli%init(progname='fossil_test_rotate',                              &
                authors='S. Zaghi',                                         &
                help='Usage: ',                                             &
                examples=["fossil_test_rotate --stl src/tests/dragon.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--stl',               &
               help='STL (input) file name', &
               required=.false.,             &
               def='src/tests/cube.stl',     &
               act='store')

  call cli%add(switch='--axis',         &
               help='axis of rotation', &
               required=.false.,        &
               nargs='+',               &
               def='1.0 0.0 0.0',       &
               act='store')

  call cli%add(switch='--angle',         &
               help='angle of rotation', &
               required=.false.,         &
               def='1.57079633',         &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl',   val=file_name_stl, error=error) ; if (error/=0) stop
  call cli%get(switch='--axis',  val=axis_,         error=error) ; if (error/=0) stop
  call cli%get(switch='--angle', val=angle,         error=error) ; if (error/=0) stop
  axis%x = axis_(1)
  axis%y = axis_(2)
  axis%z = axis_(3)
  endsubroutine cli_parse
endprogram fossil_test_rotate
