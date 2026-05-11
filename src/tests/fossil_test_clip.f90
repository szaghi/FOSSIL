!< FOSSIL, test clip STL.

program fossil_test_clip
!< FOSSIL, test clip STL.

use flap, only : command_line_interface
use fossil, only : surface_stl_object
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(surface_stl_object) :: surface             !< STL surface.
type(surface_stl_object) :: remainder           !< STL surface remainder after clip.
character(999)           :: file_name_stl       !< Input STL file name.
type(vector_R8P)         :: bmin, bmax          !< Bounding box extents.
logical                  :: are_tests_passed(2) !< Result of tests check.

are_tests_passed = .false.

call cli_parse
call surface%load_from_file(file_name=trim(adjustl(file_name_stl)), guess_format=.true.)

call surface%clip(bmin=bmin, bmax=bmax, remainder=remainder)
call surface%save_into_file(file_name='fossil_test_clip.stl')
associate (bmax_out => surface%get_bmax())
   are_tests_passed(1) = nint(bmax_out%x) <= 0
end associate
call remainder%save_into_file(file_name='fossil_test_clip_remainder.stl')

call surface%load_from_file(file_name=trim(adjustl(file_name_stl)), guess_format=.true., &
                            clip_min=bmin, clip_max=bmax)
associate (bmax_out => surface%get_bmax())
   are_tests_passed(2) = nint(bmax_out%x) <= 0
end associate

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli                !< Test command line interface.
  real(R8P)                    :: bmin_(3), bmax_(3) !< Bounding box extents.
  integer(I4P)                 :: error              !< Error trapping flag.

  call cli%init(progname='fossil_test_clip',                              &
                authors='S. Zaghi',                                       &
                help='Usage: ',                                           &
                examples=["fossil_test_clip --stl src/tests/dragon.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--stl',               &
               help='STL (input) file name', &
               required=.false.,             &
               def='src/tests/dragon.stl',   &
               act='store')

  call cli%add(switch='--bmin',                       &
               help='minimum extent of bounding box', &
               required=.false.,                      &
               nargs='+',                             &
               def='-15.0 -5.0 0.0',                  &
               act='store')

  call cli%add(switch='--bmax',                       &
               help='maximum extent of bounding box', &
               required=.false.,                      &
               nargs='+',                             &
               def='0.0 5.0 20.0',                    &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl',  val=file_name_stl, error=error) ; if (error/=0) stop
  call cli%get(switch='--bmin', val=bmin_,         error=error) ; if (error/=0) stop
  call cli%get(switch='--bmax', val=bmax_,         error=error) ; if (error/=0) stop
  bmin%x = bmin_(1)
  bmin%y = bmin_(2)
  bmin%z = bmin_(3)
  bmax%x = bmax_(1)
  bmax%y = bmax_(2)
  bmax%z = bmax_(3)
  endsubroutine cli_parse
endprogram fossil_test_clip
