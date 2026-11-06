!< FOSSIL test: named patches through booleans and self-intersection repair (issue #66, step 4).
!<
!< A = the unit cube, B = the unit cube moved by (0.5, 0.5, 0.5); each has the patches `low` (the faces x, y, z = min)
!< and `high` (x, y, z = max), B's named `b_low` and `b_high`. Every output facet keeps its parent's patch, so the
!< area of each patch in the result is analytic: A's `high` faces have 0.25 each inside B, B's `b_low` faces 0.25 each
!< inside A. Expected areas [low, high, b_low, b_high]:
!<   union      [3, 2.25, 2.25, 3]      intersect  [0, 0.75, 0.75, 0]
!<   difference [3, 2.25, 0.75, 0]      symdiff    [3, 3, 3, 3]
!< Also: the name tables merge by name (the header's patch, low, high, b_low, b_high), and
!< `resolve_self_intersections` of A keeps its patches with their areas.
program fossil_test_patches_boolean

use fossil, only : surface_stl_object, facet_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF, &
                   BOOL_STATUS_OK
use penf,   only : I4P, R8P

implicit none

integer(I4P), parameter  :: OPS(4) = [BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF]
real(R8P),    parameter  :: EXPECTED(4, 4) = reshape([3._R8P, 2.25_R8P, 2.25_R8P, 3._R8P, &
                                                      0._R8P, 0.75_R8P, 0.75_R8P, 0._R8P, &
                                                      3._R8P, 2.25_R8P, 0.75_R8P, 0._R8P, &
                                                      3._R8P, 3._R8P, 3._R8P, 3._R8P], [4, 4])
type(surface_stl_object) :: a, b, r
integer(I4P)             :: k, st
real(R8P)                :: ar(4)
logical                  :: ok, ok_all

ok_all = .true.
call make(a, 0._R8P, 'low', 'high')
call make(b, 0.5_R8P, 'b_low', 'b_high')
do k = 1, 4
   r = a
   call r%boolean(other=b, op=OPS(k), status=st)
   ar = areas(r)
   ok = st == BOOL_STATUS_OK .and. all(abs(ar - EXPECTED(:, k)) < 1e-12_R8P) .and. r%get_patches_number() == 5 .and. &
        r%patch_name(4) == 'b_low' .and. r%patch_name(5) == 'b_high'
   print '(A,I0,A,4F10.6,A,L1)', 'op ', k, ': areas per patch ', ar, '  ok ', ok
   ok_all = ok_all .and. ok
enddo
r = a
call r%resolve_self_intersections(status=st)
ar = areas(r)
ok = st == BOOL_STATUS_OK .and. r%get_patches_number() == 3 .and. abs(ar(1) - 3._R8P) < 1e-12_R8P .and. &
     abs(ar(2) - 3._R8P) < 1e-12_R8P
print '(A,L1)', 'resolve_self_intersections keeps the patches: ', ok
ok_all = ok_all .and. ok
print '(A,L1)', 'Are all tests passed? ', ok_all

contains
   subroutine make(s, shift, nlow, nhigh)
   !< The unit cube moved by (shift, shift, shift), faces at the minimum in patch `nlow`, at the maximum in `nhigh`.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: shift
   character(*),             intent(in)    :: nlow, nhigh
   type(facet_object), pointer             :: f
   integer(I4P), allocatable               :: lo(:), hi(:)
   integer(I4P)                            :: i

   call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call s%translate(x=shift, y=shift, z=shift) ; call s%analyze
   allocate(lo(0), hi(0))
   do i = 1, s%get_facets_number()
      f => s%facet_at(i)
      if (f%normal%x + f%normal%y + f%normal%z < 0._R8P) then
         lo = [lo, i]
      else
         hi = [hi, i]
      endif
   enddo
   call s%set_patch(lo, nlow)
   call s%set_patch(hi, nhigh)
   endsubroutine make

   function areas(s) result(ar)
   !< Area per patch name [low, high, b_low, b_high] (the tables also hold the header's patch, left empty).
   type(surface_stl_object), intent(in) :: s
   real(R8P)                            :: ar(4)
   character(6), parameter              :: NM(4) = ['low   ', 'high  ', 'b_low ', 'b_high']
   type(facet_object), pointer          :: f
   integer(I4P)                         :: i, j

   ar = 0._R8P
   do i = 1, s%get_facets_number()
      f => s%facet_at(i)
      do j = 1, 4
         if (s%patch_name(f%patch) == trim(NM(j))) ar(j) = ar(j) + f%area()
      enddo
   enddo
   endfunction areas
endprogram fossil_test_patches_boolean
