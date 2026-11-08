program ck_self_intersections
use fossil, only : surface_stl_object, intersection_pair_t, BOOL_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)               :: a, b, soup
type(intersection_pair_t), allocatable :: pairs(:)
integer(I4P)                           :: status

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
b = a
call b%translate(delta=vector_R8P(0.3_R8P, 0.35_R8P, 0.4_R8P)) ; call b%analyze
soup = a
call soup%merge_solids(other=b) ; call soup%analyze          ! two crossing cubes in one surface
call soup%find_self_intersections(pairs=pairs)
print '(A,I0,A,I0,A,F7.4)', 'merged: ', soup%get_facets_number(), ' facets, ', size(pairs), &
                            ' crossing pairs, volume ', soup%get_volume()
print '(A,2(3F7.3,2X))', 'first crossing segment ', pairs(1)%p%x, pairs(1)%p%y, pairs(1)%p%z, &
                         pairs(1)%q%x, pairs(1)%q%y, pairs(1)%q%z
call soup%resolve_self_intersections(status=status)          ! remove them: keep only the outer surface
if (status /= BOOL_STATUS_OK) error stop 'resolve_self_intersections failed'
call soup%find_self_intersections(pairs=pairs)
print '(A,I0,A,I0,A,F7.4,A,L1)', 'resolved: ', soup%get_facets_number(), ' facets, ', size(pairs), &
                                 ' crossing pairs, volume ', soup%get_volume(), ', watertight ', soup%is_watertight()
endprogram ck_self_intersections
