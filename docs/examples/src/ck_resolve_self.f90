!run ck_resolve_self ck_resolve_self
program ck_resolve_self
use fossil, only : surface_stl_object, facet_object, intersection_pair_t, BOOL_UNION, BOOL_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)               :: soup, resolved, a, b
type(intersection_pair_t), allocatable :: pairs(:)
integer(I4P)                           :: status

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
b = a
call b%translate(delta=vector_R8P(0.3_R8P, 0.3_R8P, 0.3_R8P)) ; call b%analyze
soup = a
call soup%merge_solids(other=b) ; call soup%analyze          ! one surface, two crossing shells
resolved = soup
call resolved%resolve_self_intersections(status=status)
call resolved%find_self_intersections(pairs=pairs)
print '(A,I0,A,F7.4,A,L1,A,I0)', 'resolve_self_intersections: status ', status, ', volume ', resolved%get_volume(), &
                                 ', watertight ', resolved%is_watertight(), ', crossing pairs ', size(pairs)
if (.not. resolved%is_watertight()) then                      ! check the result: until issue #78 is fixed it is
   call shell(soup, 1_I4P, a) ; call shell(soup, 2_I4P, b)    ! the arrangement; unite the shells instead
   call a%boolean(other=b, op=BOOL_UNION, status=status)
   if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
   resolved = a
endif
call resolved%find_self_intersections(pairs=pairs)
print '(A,I0,A,F7.4,A,L1,A,I0)', 'result: facets ', resolved%get_facets_number(), ', volume ', resolved%get_volume(), &
                                 ', watertight ', resolved%is_watertight(), ', crossing pairs ', size(pairs)
contains
   subroutine shell(s, k, part)
   !< The facets of shell k of s, as a surface of their own.
   type(surface_stl_object), intent(in)  :: s
   integer(I4P),             intent(in)  :: k
   type(surface_stl_object), intent(out) :: part
   type(facet_object), allocatable       :: f(:)
   integer(I4P)                          :: i, n

   allocate(f(count([(s%facet_shell(i) == k, i = 1, s%get_facets_number())]))) ; n = 0
   do i = 1, s%get_facets_number()
      if (s%facet_shell(i) /= k) cycle
      n = n + 1 ; f(n) = s%facet_at(i)
   enddo
   call part%adopt_facets(facets=f)
   endsubroutine shell
endprogram ck_resolve_self
