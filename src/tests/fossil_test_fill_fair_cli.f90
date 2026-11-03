!< FOSSIL test: refine/fair options through sanitize and fossilizer (issue #31, step 5).
program fossil_test_fill_fair_cli
!< FOSSIL test: refine/fair options through sanitize and fossilizer (issue #31, step 5).
!<
!< cube-02 with a 3 x 3 block of squares removed from the face z = 1 (tiny):
!<  - `sanitize(fill_holes=.true., refine_holes=.true.)` and `(..., fair_holes=.true.)` close
!<    the hole: watertight, volume 1 to 1e-12, more facets than the plain fill. The hole touches
!<    two edges of the cube: fairing must not bend the patch toward the side faces across those
!<    90-degree creases (the feature-aware ring), so the patch stays in its plane;
!<  - `fossilizer -i holed.stl --fill_holes --fair_holes --max_hole_size 1 --output out` writes a
!<    watertight surface with more facets than the plain fill.
!< Runs exe/fossilizer from the repository root, like every FOSSIL test.

use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P

implicit none

character(*), parameter         :: HOLED = 'fossil_test_fillfair_holed.stl'
character(*), parameter         :: OUT   = 'fossil_test_fillfair_out.stl'
type(surface_stl_object)        :: cube, s_plain, s_ref, s_fair, s_cli
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
integer(I4P)                    :: i, n, exitstat, u
logical                         :: ok

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
allocate(facets(cube%get_facets_number())) ; n = 0
do i = 1, cube%get_facets_number()
   f => cube%facet_at(i)
   if (f%centroid%z > 1._R8P - 1.e-9_R8P .and. f%centroid%x < 0.75_R8P .and. f%centroid%y < 0.75_R8P) cycle
   n = n + 1 ; facets(n) = f ; facets(n)%vertex_id = 0_I4P
enddo
facets = facets(1:n)
call s_plain%adopt_facets(facets=facets)
call s_plain%save_into_file(file_name=HOLED)
call s_ref%load_from_file(file_name=HOLED, guess_format=.true.)
call s_fair%load_from_file(file_name=HOLED, guess_format=.true.)
call s_plain%sanitize(fill_holes=.true., max_hole_size=1._R8P)
call s_ref%sanitize(fill_holes=.true., max_hole_size=1._R8P, refine_holes=.true.)
call s_fair%sanitize(fill_holes=.true., max_hole_size=1._R8P, fair_holes=.true.)
ok = s_plain%is_watertight() .and. s_ref%is_watertight() .and. s_fair%is_watertight() .and. &
     abs(s_ref%get_volume() - 1._R8P) < 1.e-12_R8P .and. abs(s_fair%get_volume() - 1._R8P) < 1.e-12_R8P .and. &
     s_ref%get_facets_number() > s_plain%get_facets_number() .and. s_fair%get_facets_number() > s_plain%get_facets_number()
print '(A,3I5,A,2F16.12)', 'sanitize facets plain/refine/fair: ', s_plain%get_facets_number(), s_ref%get_facets_number(), &
      s_fair%get_facets_number(), '  volumes refine/fair: ', s_ref%get_volume(), s_fair%get_volume()

call execute_command_line('exe/fossilizer -i '//HOLED//' --fill_holes --fair_holes --max_hole_size 1 --output '//OUT// &
                          ' > /dev/null', exitstat=exitstat)
call s_cli%load_from_file(file_name=OUT, guess_format=.true.)
print '(A,I0,A,I0,A,L1)', 'fossilizer --fair_holes: exit ', exitstat, '  facets ', s_cli%get_facets_number(), &
      '  watertight ', s_cli%is_watertight()
ok = ok .and. exitstat == 0 .and. s_cli%is_watertight() .and. s_cli%get_facets_number() > s_plain%get_facets_number()

open(newunit=u, file=HOLED, status='old') ; close(u, status='delete')
open(newunit=u, file=OUT, status='old') ; close(u, status='delete')
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1
endprogram fossil_test_fill_fair_cli
