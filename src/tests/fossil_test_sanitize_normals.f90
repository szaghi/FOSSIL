!< FOSSIL, test normals *sanitization*.
!<
!< The default input, cube-inconsistent.stl, is the unit cube with one facet wound the
!< wrong way: its signed volume is 1 - 2*(1/6) = 2/3 before the fix. sanitize_normals
!< must report at least one flipped facet and restore volume 1.

program fossil_test_sanitize_normals
!< FOSSIL, test normals *sanitization*.

use flap, only : command_line_interface
use fossil, only : surface_stl_object
use penf, only : I4P, R8P

implicit none

type(surface_stl_object) :: surface             !< STL surface.
character(999)           :: file_name_stl       !< Input STL file name.
logical                  :: are_tests_passed(3) !< Result of tests check.
integer(I4P)             :: flipped_count       !< Facets reversed by sanitize_normals.

are_tests_passed = .false.

call cli_parse
call surface%load_from_file(file_name=trim(adjustl(file_name_stl)), guess_format=.true.)

print*, 'volume before sanitize normals: ', surface%get_volume()
are_tests_passed(1) = abs(surface%get_volume() - 1._R8P) > 0.1_R8P   ! the defect is visible in the volume
call surface%sanitize_normals(flipped_count=flipped_count)
call surface%compute_volume
print*, 'volume after sanitize normals:  ', surface%get_volume()
print*, 'facets flipped:                 ', flipped_count
are_tests_passed(2) = abs(surface%get_volume() - 1._R8P) < 1.e-12_R8P
are_tests_passed(3) = flipped_count > 0_I4P

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. (all(are_tests_passed))) error stop 1
contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli   !< Test command line interface.
  integer(I4P)                 :: error !< Error trapping flag.

  call cli%init(progname='fossil_test_sanitize_normals',                              &
                authors='S. Zaghi',                                                   &
                help='Usage: ',                                                       &
                examples=["fossil_test_sanitize_normals --stl src/tests/dragon.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--stl',                        &
               help='STL (input) file name',          &
               required=.false.,                      &
               def='src/tests/cube-inconsistent.stl', &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl', val=file_name_stl, error=error) ; if (error/=0) stop
  endsubroutine cli_parse
endprogram fossil_test_sanitize_normals
