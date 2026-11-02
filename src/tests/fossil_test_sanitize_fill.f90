!< FOSSIL, test `sanitize(fill_holes=...)` (issue #31, step 4).
!<
!<   a. cube minus 1 facet, default sanitize: the hole stays (behaviour unchanged).
!<   b. same with fill_holes=.true.: watertight, volume 1, centre inside.
!<   c. inside-out holed cube (reversed, analyzed) with filling: watertight, outward, volume 1.
!<   d. cube-02 with a 2x2 hole and 5 facets far from it reversed (mixed winding): filled after
!<      orientation, watertight, volume 1, centre inside.
!< Holes on the unit cube are large relative to the surface, so max_hole_size is passed.

program fossil_test_sanitize_fill

use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

logical                  :: ok(4)
type(surface_stl_object) :: s
real(R8P)                :: d

call holed(s, 'src/tests/cube.stl', 0.1_R8P, .false., 0)
call s%sanitize
ok(1) = (.not. s%is_watertight()) .and. s%get_facets_number() == 11
print '(A,I0,A,L1,A,L1)', 'a default sanitize: facets=', s%get_facets_number(), ' watertight=', s%is_watertight(), &
      ' ok=', ok(1)

call holed(s, 'src/tests/cube.stl', 0.1_R8P, .false., 0)
call s%sanitize(fill_holes=.true., max_hole_size=1._R8P)
ok(2) = report(s, 'b fill_holes=.true.        ')

call holed(s, 'src/tests/cube.stl', 0.1_R8P, .true., 0)
call s%sanitize(fill_holes=.true., max_hole_size=1._R8P)
ok(3) = report(s, 'c inside-out holed cube    ')

call holed(s, 'src/tests/cube-02.stl', 0.5_R8P, .false., 5)
call s%sanitize(fill_holes=.true., max_hole_size=0.5_R8P)
ok(4) = report(s, 'd mixed winding + 2x2 hole ')

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine holed(s, file_name, cell, invert, n_reverse)
   !< Remove the facets of the z = 0 face with centroid in [0, cell]^2 (cube.stl: facet 1);
   !< optionally reverse every facet, or the last n_reverse facets; adopt (analyze).
   type(surface_stl_object), intent(inout) :: s
   character(*),             intent(in)    :: file_name
   real(R8P),                intent(in)    :: cell
   logical,                  intent(in)    :: invert
   integer(I4P),             intent(in)    :: n_reverse
   type(surface_stl_object)                :: src
   type(facet_object), allocatable         :: keep(:)
   type(facet_object), pointer             :: fp
   integer(I4P)                            :: i, n

   call src%load_from_file(file_name=file_name, guess_format=.true.)
   call src%sanitize
   allocate(keep(src%get_facets_number())) ; n = 0
   do i = 1, src%get_facets_number()
      fp => src%facet_at(i)
      if (src%get_facets_number() == 12) then
         if (i == 1) cycle
      elseif (abs(fp%centroid%z) < 1.e-12_R8P .and. fp%centroid%x < cell .and. fp%centroid%y < cell) then
         cycle
      endif
      n = n + 1 ; keep(n) = fp
   enddo
   keep = keep(1:n)
   if (invert) call keep%reverse_normal
   do i = n - n_reverse + 1, n
      call keep(i)%reverse_normal
   enddo
   call s%adopt_facets(facets=keep)
   endsubroutine holed

   function report(s, label) result(passed)
   type(surface_stl_object), intent(in) :: s
   character(*),             intent(in) :: label
   logical                              :: passed
   real(R8P)                            :: d
   d = s%distance(point=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), is_signed=.true., is_square_root=.true.)
   passed = s%is_watertight() .and. abs(s%get_volume() - 1._R8P) < 1.e-12_R8P .and. abs(d + 0.5_R8P) < 1.e-12_R8P
   print '(A,A,I0,A,L1,A,F9.6,A,F8.4,A,L1)', label, ': facets=', s%get_facets_number(), ' watertight=', s%is_watertight(), &
         ' volume=', s%get_volume(), ' d(centre)=', d, ' ok=', passed
   endfunction report
endprogram fossil_test_sanitize_fill
