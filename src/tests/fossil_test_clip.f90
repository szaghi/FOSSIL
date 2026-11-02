!< FOSSIL, test clip STL.
!<
!< `clip` keeps every facet that overlaps the box, whole (inclusion by overlap, not by
!< containment), so the kept surface may extend past the box. Checks:
!<   1. every kept facet's bounding box intersects the clip box (necessary for overlap);
!<   2. kept + remainder = original facet count;
!<   3. no remainder facet lies entirely inside the box;
!<   4. `load_from_file(clip_min=, clip_max=)` keeps the same number of facets as `clip`.

program fossil_test_clip
!< FOSSIL, test clip STL.

use flap, only : command_line_interface
use fossil, only : surface_stl_object, facet_object
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

type(surface_stl_object) :: surface             !< STL surface.
type(surface_stl_object) :: remainder           !< STL surface remainder after clip.
character(999)           :: file_name_stl       !< Input STL file name.
type(vector_R8P)         :: bmin, bmax          !< Bounding box extents.
logical                  :: are_tests_passed(4) !< Result of tests check.
integer(I4P)             :: n_original          !< Facets before clipping.
integer(I4P)             :: n_kept              !< Facets kept by clip.

are_tests_passed = .false.

call cli_parse
call surface%load_from_file(file_name=trim(adjustl(file_name_stl)), guess_format=.true.)
n_original = surface%get_facets_number()

call surface%clip(bmin=bmin, bmax=bmax, remainder=remainder)
call surface%save_into_file(file_name='fossil_test_clip.stl')
n_kept = surface%get_facets_number()
are_tests_passed(1) = n_kept > 0 .and. all_overlap_box(surface)
are_tests_passed(2) = n_kept + remainder%get_facets_number() == n_original
are_tests_passed(3) = none_inside_box(remainder)
call remainder%save_into_file(file_name='fossil_test_clip_remainder.stl')

call surface%load_from_file(file_name=trim(adjustl(file_name_stl)), guess_format=.true., &
                            clip_min=bmin, clip_max=bmax)
are_tests_passed(4) = surface%get_facets_number() == n_kept .and. all_overlap_box(surface)

print '(A,I0,A,I0,A,I0,A,4L2)', 'facets: original=', n_original, ' kept=', n_kept, &
      ' remainder=', remainder%get_facets_number(), '  checks:', are_tests_passed

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
contains
  function all_overlap_box(s) result(ok)
  !< True if every facet's bounding box intersects [bmin, bmax].
  type(surface_stl_object), intent(in) :: s  !< Surface.
  logical                              :: ok !< Check result.
  type(facet_object), pointer          :: f  !< Facet.
  integer(I4P)                         :: i  !< Counter.

  ok = .true.
  do i = 1, s%get_facets_number()
     f => s%facet_at(i)
     if (maxval(f%vertex(:)%x) < bmin%x .or. minval(f%vertex(:)%x) > bmax%x .or. &
         maxval(f%vertex(:)%y) < bmin%y .or. minval(f%vertex(:)%y) > bmax%y .or. &
         maxval(f%vertex(:)%z) < bmin%z .or. minval(f%vertex(:)%z) > bmax%z) ok = .false.
  enddo
  endfunction all_overlap_box

  function none_inside_box(s) result(ok)
  !< True if no facet has all three vertices inside [bmin, bmax] (such a facet overlaps the box).
  type(surface_stl_object), intent(in) :: s  !< Surface.
  logical                              :: ok !< Check result.
  type(facet_object), pointer          :: f  !< Facet.
  integer(I4P)                         :: i  !< Counter.

  ok = .true.
  do i = 1, s%get_facets_number()
     f => s%facet_at(i)
     if (all(f%vertex(:)%x >= bmin%x) .and. all(f%vertex(:)%x <= bmax%x) .and. &
         all(f%vertex(:)%y >= bmin%y) .and. all(f%vertex(:)%y <= bmax%y) .and. &
         all(f%vertex(:)%z >= bmin%z) .and. all(f%vertex(:)%z <= bmax%z)) ok = .false.
  enddo
  endfunction none_inside_box

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
