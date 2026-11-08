!run ck_resolve_self ck_resolve_self
program ck_resolve_self
use fossil, only : surface_stl_object, facet_object, intersection_pair_t, BOOL_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)               :: cube
type(facet_object), allocatable        :: facets(:)
type(facet_object), pointer            :: f
type(intersection_pair_t), allocatable :: pairs(:)
integer(I4P)                           :: status, i, j

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
allocate(facets(cube%get_facets_number()))
do i = 1, size(facets)                                       ! push the corner (1, 1, 1) through the bottom face
   f => cube%facet_at(i) ; facets(i) = f
   do j = 1, 3
      if (facets(i)%vertex(j)%x == 1._R8P .and. facets(i)%vertex(j)%y == 1._R8P .and. facets(i)%vertex(j)%z == 1._R8P) &
         facets(i)%vertex(j) = vector_R8P(0.3_R8P, 0.35_R8P, -0.3_R8P)
   enddo
   call facets(i)%compute_normal
enddo
call cube%adopt_facets(facets=facets)                        ! one closed shell passing through itself
call cube%find_self_intersections(pairs=pairs)
print '(A,I0,A,L1,A,I0,A,F7.4)', 'folded:   ', cube%get_facets_number(), ' facets, watertight ', &
                                 cube%is_watertight(), ', crossing pairs ', size(pairs), ', volume ', cube%get_volume()
call cube%resolve_self_intersections(status=status)
if (status /= BOOL_STATUS_OK) then                          ! CDT_FAILED or NOT_CLOSED: the surface is unchanged
   print '(A,I0)', 'resolve_self_intersections failed, status ', status
   stop
endif
call cube%find_self_intersections(pairs=pairs)
print '(A,I0,A,L1,A,I0,A,F7.4)', 'resolved: ', cube%get_facets_number(), ' facets, watertight ', &
                                 cube%is_watertight(), ', crossing pairs ', size(pairs), ', volume ', cube%get_volume()
endprogram ck_resolve_self
