!< FOSSIL, test merge STL.

program fossil_test_merge
!< FOSSIL, test merge STL.

use flap, only : command_line_interface
use fossil, only : surface_stl_object
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

character(999)           :: file_name_stl(2)    !< Input STL file names.
type(surface_stl_object) :: surface(2)          !< STL surface.
logical                  :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call cli_parse
call surface(1)%load_from_file(file_name=trim(adjustl(file_name_stl(1))), guess_format=.true.)
call surface(2)%load_from_file(file_name=trim(adjustl(file_name_stl(2))), guess_format=.true.)
call surface(1)%merge_solids(other=surface(2))
call surface(1)%save_into_file(file_name='fossil_test_merge.stl')
associate (bmin_out => surface(1)%get_bmin())
   are_tests_passed(1) = nint(bmin_out%x) < -2
end associate

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli   !< Test command line interface.
  integer(I4P)                 :: error !< Error trapping flag.

  call cli%init(progname='fossil_test_merge',                                                     &
                authors='S. Zaghi',                                                               &
                help='Usage: ',                                                                   &
                examples=["fossil_test_merge --stl1 dragon_part_1.stl --stl2 dragon_part_2.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--stl1',                   &
               help='STL 1 (input) file name',    &
               required=.false.,                  &
               def='src/tests/dragon_part_1.stl', &
               act='store')

  call cli%add(switch='--stl2',                   &
               help='STL 2 (input) file name',    &
               required=.false.,                  &
               def='src/tests/dragon_part_2.stl', &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl1', val=file_name_stl(1), error=error) ; if (error/=0) stop
  call cli%get(switch='--stl2', val=file_name_stl(2), error=error) ; if (error/=0) stop
  endsubroutine cli_parse
endprogram fossil_test_merge
