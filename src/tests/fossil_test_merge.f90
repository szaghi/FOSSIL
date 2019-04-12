!< FOSSIL, test merge STL.

program fossil_test_merge
!< FOSSIL, test merge STL.

use flap, only : command_line_interface
use fossil, only : file_stl_object, surface_stl_object
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(file_stl_object)    :: file_stl            !< STL file handler.
character(999)           :: file_name_stl(2)    !< Input STL file names.
type(surface_stl_object) :: surface(2)          !< STL surface.
logical                  :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call cli_parse
call file_stl%load_from_file(facet=surface(1)%facet, file_name=trim(adjustl(file_name_stl(1))), guess_format=.true.)
call file_stl%load_from_file(facet=surface(2)%facet, file_name=trim(adjustl(file_name_stl(2))), guess_format=.true.)
! call surface%analize
call surface(1)%merge_solids(other=surface(2))
call file_stl%save_into_file(facet=surface(1)%facet, file_name='fossil_test_merge.stl')
are_tests_passed(1) = nint(surface(1)%bmax%x) > 0

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli   !< Test command line interface.
  integer(I4P)                 :: error !< Error trapping flag.

  call cli%init(progname='fossil_test_merge',                                             &
                authors='S. Zaghi',                                                       &
                help='Usage: ',                                                           &
                examples=["fossil_test_merge --stl dragon_part_1.stl dragon_part_2.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--stl',                                                  &
               help='STL (input) file name',                                    &
               required=.false.,                                                &
               nargs='2',                                                       &
               def='src/tests/dragon_part_1.stl src/tests/dragon_part_2.stl',   &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl', val=file_name_stl, error=error) ; if (error/=0) stop
  endsubroutine cli_parse
endprogram fossil_test_merge
