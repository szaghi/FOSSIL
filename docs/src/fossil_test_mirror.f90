!< FOSSIL, test mirror STL.

program fossil_test_mirror
!< FOSSIL, test mirror STL.

use flap, only : command_line_interface
use fossil, only : file_stl_object
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(file_stl_object) :: file_stl            !< STL file.
character(999)        :: file_name_stl       !< Input STL file name.
type(vector_R8P)      :: normal              !< Normal of mirroring plane.
logical               :: are_tests_passed(2) !< Result of tests check.

are_tests_passed = .false.

call cli_parse
call file_stl%initialize(file_name=trim(adjustl(file_name_stl)))
call file_stl%load_from_file(guess_format=.true.)
print '(A)', file_stl%statistics()

call file_stl%mirror(normal=normal)
call file_stl%save_into_file(file_name='fossil_test_mirror.stl')
are_tests_passed(1) = nint(file_stl%distance(point=0 * ex_R8P - 1 * ey_R8P + 1 * ez_R8P)) == 1
are_tests_passed(2) = nint(file_stl%distance(point=1 * ex_R8P - 1 * ey_R8P + 0 * ez_R8P)) == 1

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli        !< Test command line interface.
  real(R8P)                    :: normal_(3) !< Normal of mirroring plane.
  integer(I4P)                 :: error      !< Error trapping flag.

  call cli%init(progname='fossil_test_mirror',                              &
                authors='S. Zaghi',                                         &
                help='Usage: ',                                             &
                examples=["fossil_test_mirror --stl src/tests/dragon.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--stl',               &
               help='STL (input) file name', &
               required=.false.,             &
               def='src/tests/cube.stl',     &
               act='store')

  call cli%add(switch='--normal',                &
               help='normal of mirroring plane', &
               required=.false.,                 &
               nargs='+',                        &
               def='1.0 0.0 0.0',                &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl',    val=file_name_stl, error=error) ; if (error/=0) stop
  call cli%get(switch='--normal', val=normal_,       error=error) ; if (error/=0) stop
  normal%x = normal_(1)
  normal%y = normal_(2)
  normal%z = normal_(3)
  endsubroutine cli_parse
endprogram fossil_test_mirror
